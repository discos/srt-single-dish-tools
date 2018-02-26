from __future__ import division, print_function
from srttools import CalibratorTable
from srttools.calibration import main_lcurve, _get_flux_quantity, main_cal
from srttools.read_config import read_config
from srttools.scan import list_scans
from srttools.simulate import sim_crossscans, _2d_gauss
from srttools.io import mkdir_p
from srttools.utils import compare_strings, HAS_MPL
import pytest
import logging
import glob

import os
import numpy as np
import subprocess as sp

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(x):
        return x


@pytest.fixture()
def logger():
    logger = logging.getLogger('Some.Logger')
    logger.setLevel(logging.DEBUG)

    return logger


np.random.seed(1241347)


def source_scan_func(x):
    return 52 * _2d_gauss(x, 0, sigma=2.5 / 60)


class TestCalibration(object):
    @classmethod
    def setup_class(klass):
        import os

        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, 'data')

        klass.config_file = \
            os.path.abspath(os.path.join(klass.datadir, "calibrators.ini"))

        klass.config_file_empty = \
            os.path.abspath(os.path.join(klass.datadir,
                                         "calibrators_nocal.ini"))

        klass.config = read_config(klass.config_file)
        klass.caldir = os.path.join(klass.datadir, 'sim', 'calibration')
        klass.caldir2 = os.path.join(klass.datadir, 'sim', 'calibration2')
        klass.caldir3 = os.path.join(klass.datadir, 'sim', 'calibration_bad')
        klass.crossdir = os.path.join(klass.datadir, 'sim', 'crossscans')
        if not os.path.exists(klass.caldir):
            print('Fake calibrators: DummyCal, 1 Jy.')
            mkdir_p(klass.caldir)
            sim_crossscans(5, klass.caldir)
        if not os.path.exists(klass.caldir2):
            print('Fake calibrators: DummyCal2, 1 Jy.')
            mkdir_p(klass.caldir2)
            sim_crossscans(5, klass.caldir2, srcname='DummyCal2')
        if not os.path.exists(klass.caldir3):
            print('Fake calibrators: DummyCal2, wrong flux 0.52 Jy.')
            mkdir_p(klass.caldir3)
            sim_crossscans(1, klass.caldir3, srcname='DummyCal2',
                           scan_func=source_scan_func)
        if not os.path.exists(klass.crossdir):
            print('Fake cross scans: DummySrc, 0.52 Jy.')
            mkdir_p(klass.crossdir)
            sim_crossscans(5, klass.crossdir, srcname='DummySrc',
                           scan_func=source_scan_func)

        klass.scan_list = \
            list_scans(klass.caldir, ['./']) + \
            list_scans(klass.caldir2, ['./']) + \
            list_scans(klass.crossdir, ['./'])

        klass.scan_list.sort()
        caltable = CalibratorTable()
        caltable.from_scans(klass.scan_list)
        caltable.update()

        klass.calfile = os.path.join(klass.curdir, 'test_calibration.hdf5')
        caltable.write(klass.calfile, overwrite=True)
        caltable.write(klass.calfile.replace('hdf5', 'csv'))

    def test_0_prepare(self):
        pass

    def test_script_is_installed(self):
        sp.check_call('SDTcal -h'.split(' '))

    def test_check_not_empty(self):
        caltable = CalibratorTable()
        assert not caltable.check_not_empty()

    def test_check_up_to_date_empty_return_false(self):
        caltable = CalibratorTable()
        assert not caltable.check_up_to_date()

    def test_calibrate_empty_return_none(self):
        caltable = CalibratorTable()
        assert caltable.calibrate() is None

    def test_update_empty_return_none(self):
        caltable = CalibratorTable()
        assert caltable.update() is None

    def test_get_fluxes_empty_return_none(self):
        caltable = CalibratorTable()
        assert caltable.get_fluxes() is None

    def test_check_class(self):
        caltable = CalibratorTable()
        caltable.from_scans(self.scan_list)
        with pytest.warns(UserWarning):
            caltable.check_up_to_date()

    def test_check_class_from_file(self):
        caltable = CalibratorTable.read(self.calfile, path='table')
        assert caltable.check_up_to_date()

    def test_Jy_over_counts_and_back(self):
        caltable = CalibratorTable.read(self.calfile, path='table')
        Jc, Jce = caltable.Jy_over_counts(channel='Feed0_LCP')
        Cj, Cje = caltable.counts_over_Jy(channel='Feed0_LCP')
        np.testing.assert_allclose(Jc, 1 / Cj)

    def test_Jy_over_counts_rough_one_bad_value(self, logger, caplog):
        caltable = CalibratorTable.read(self.calfile, path='table')

        flux_quantity = _get_flux_quantity('Jy/beam')
        caltable[flux_quantity + "/Counts"][0] += \
            caltable[flux_quantity + "/Counts Err"][0] * 2000
        Jc, Jce = caltable.Jy_over_counts_rough(channel='Feed0_LCP',
                                                map_unit='Jy/beam')
        assert 'Outliers: ' in caplog.text
        Cj, Cje = caltable.counts_over_Jy(channel='Feed0_LCP')
        np.testing.assert_allclose(Jc, 1 / Cj)

    def test_bad_file_missing_key(self, logger, caplog):
        caltable = CalibratorTable()
        caltable.from_scans([os.path.join(self.config['datadir'],
                                          'calibrators', 'summary.fits')])
        assert "Missing key" in caplog.text

    def test_bad_file_generic_error(self, logger, caplog):
        caltable = CalibratorTable()

        caltable.from_scans([os.path.join(self.config['datadir'],
                                          'calibrators', 'bubu.fits')])
        assert "Error while processing" in caplog.text

    def test_calibration_counts(self):
        """Simple calibration from scans."""

        caltable = CalibratorTable.read(self.calfile, path='table')
        caltable = caltable[compare_strings(caltable['Source'], 'DummyCal')]
        caltable_0 = caltable[compare_strings(caltable['Chan'], 'Feed0_LCP')]
        assert np.all(
            np.abs(caltable_0['Counts'] - 100.) < 3 * caltable_0['Counts Err'])
        caltable_1 = caltable[compare_strings(caltable['Chan'], 'Feed0_RCP')]
        assert np.all(
            np.abs(caltable_1['Counts'] - 80.) < 3 * caltable_1['Counts Err'])

    def test_calibration_width(self):
        """Simple calibration from scans."""

        caltable = CalibratorTable.read(self.calfile, path='table')
        caltable0 = caltable[compare_strings(caltable['Chan'], 'Feed0_LCP')]
        assert np.all(
            np.abs(caltable0['Width'] - 2.5/60.) < 5 * caltable0['Width Err'])
        caltable1 = caltable[compare_strings(caltable['Chan'], 'Feed0_RCP')]
        assert np.all(
            np.abs(caltable1['Width'] - 2.5/60.) < 5 * caltable1['Width Err'])

        beam, beam_err = caltable.beam_width(channel='Feed0_LCP')
        assert np.all(beam - np.radians(2.5/60) < 3 * beam_err)

    @pytest.mark.skipif('not HAS_MPL')
    def test_calibration_plot_two_cols(self):
        """Simple calibration from scans."""

        caltable = CalibratorTable.read(self.calfile, path='table')
        caltable.plot_two_columns('RA', "Flux/Counts", xerrcol="RA err",
                                  yerrcol="Flux/Counts Err", test=True)

    @pytest.mark.skipif('not HAS_MPL')
    def test_calibration_show(self):
        """Simple calibration from scans."""

        caltable = CalibratorTable.read(self.calfile, path='table')

        caltable.show()

    def test_calibrated_crossscans(self):
        caltable = CalibratorTable.read(self.calfile, path='table')
        dummy_flux, dummy_flux_err = \
            caltable.calculate_src_flux(source='DummySrc', channel='Feed0_LCP')
        assert (dummy_flux[0] - 0.52) < dummy_flux_err[0] * 3

    def test_check_consistency_fails_with_bad_data(self):
        scan_list = \
            list_scans(self.caldir, ['./']) + \
            list_scans(self.caldir2, ['./']) + \
            list_scans(self.caldir3, ['./']) + \
            list_scans(self.crossdir, ['./'])

        scan_list.sort()
        caltable = CalibratorTable()
        caltable.from_scans(scan_list)
        caltable.update()
        res = caltable.check_consistency(channel='Feed0_LCP')
        assert not np.all(res)

    def test_check_consistency(self):
        caltable = CalibratorTable.read(self.calfile, path='table')
        res = caltable.check_consistency(channel='Feed0_LCP')
        assert np.all(res)
        res = caltable.check_consistency(channel='Feed0_RCP')
        assert np.all(res)
        res = caltable.check_consistency()
        assert np.all(res)

    @pytest.mark.skipif('not HAS_MPL')
    def test_sdtcal_with_calfile(self):
        if os.path.exists("calibration_summary.png"):
            os.unlink("calibration_summary.png")
        with pytest.raises(SystemExit):
            main_cal([self.calfile])
        assert os.path.exists("calibration_summary.png")
        os.unlink("calibration_summary.png")

    @pytest.mark.skipif('not HAS_MPL')
    def test_sdtcal_show_with_config(self):
        main_cal(('-c ' + self.config_file + ' --check --show').split(" "))
        assert os.path.exists(self.config_file.replace(".ini", "_cal.hdf5"))
        assert os.path.exists("calibration_summary.png")

    def test_sdtcal_with_sample_config(self):
        if os.path.exists('sample_config_file.ini'):
            os.unlink('sample_config_file.ini')
        with pytest.raises(SystemExit):
            main_cal(['--sample-config'])
        assert os.path.exists('sample_config_file.ini')

    def test_sdtcal_no_config(self):
        # ValueError("Please specify the config file!")
        with pytest.raises(ValueError) as excinfo:
            main_cal([])
            assert "Please specify the config file!" in str(excinfo)

    def test_sdtcal_no_config_dir(self):
        ValueError("No calibrators specified in config file")
        with pytest.raises(ValueError) as excinfo:
            main_cal(["-c", self.config_file_empty])
            assert "No calibrators specified in config file" in str(excinfo)

    def test_lcurve_with_single_source(self):
        main_lcurve([self.calfile, '-s', 'DummySrc'])
        assert os.path.exists('DummySrc.csv')
        os.unlink('DummySrc.csv')

    def test_lcurve_with_all_sources(self):
        main_lcurve(['-c', self.config_file])
        assert os.path.exists('DummySrc.csv')
        assert os.path.exists('DummyCal.csv')
        assert os.path.exists('DummyCal2.csv')

    @classmethod
    def teardown_class(klass):
        """Clean up the mess."""
        if HAS_MPL:
            os.unlink('calibration_summary.png')
        for d in klass.config['list_of_directories']:
            hfiles = \
                glob.glob(os.path.join(klass.config['datadir'], d, '*.hdf5'))
            for h in hfiles:
                os.unlink(h)

            dirs = \
                glob.glob(os.path.join(klass.config['datadir'], d,
                                       '*_scanfit'))
            for dirname in dirs:
                shutil.rmtree(dirname)
