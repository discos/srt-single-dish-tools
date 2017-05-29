from __future__ import division, print_function
from ..calibration import CalibratorTable
from ..read_config import read_config
from ..scan import list_scans
from ..simulate import save_scan
from ..io import mkdir_p
import pytest

import os
import glob
import shutil
import numpy as np
import numpy.random as ra

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(x):
        return x


def _2d_gauss(x, y, sigma=3 / 60.):
    """A Gaussian beam"""
    return np.exp(-(x ** 2 + y ** 2) / (2 * sigma**2))

np.random.seed(1241347)


def calibrator_scan_func(x):
    return 100 * _2d_gauss(x, 0, sigma=3 / 60)


def sim_calibrators(ncross, caldir):
    src_ra = 185
    src_dec = 75
    timedelta = 0
    speed = 2.  # arcmin/s
    dt = 0.04
    dtheta = speed * dt
    scan_values = np.arange(-2, 2, dtheta/60)
    zero_values = np.zeros_like(scan_values)

    for i in tqdm(range(ncross)):
        ras = src_ra + scan_values / np.cos(np.radians(src_dec))
        decs = src_dec + zero_values
        times = np.arange(scan_values.size) * dt + timedelta

        scan = calibrator_scan_func(scan_values) + \
            ra.normal(0, 0.2, scan_values.size)
        save_scan(times, ras, decs, {'Ch0': scan, 'Ch1': scan},
                  filename=os.path.join(caldir, '{}_Ra.fits'.format(i)),
                  src_ra=src_ra, src_dec=src_dec, srcname='DummyCal')
        timedelta = times[-1] + 1

        ras = src_ra + zero_values
        decs = src_dec + scan_values
        times = np.arange(scan_values.size) * dt + timedelta

        scan = calibrator_scan_func(scan_values) + \
            ra.normal(0, 0.2, scan_values.size)
        save_scan(times, ras, decs, {'Ch0': scan, 'Ch1': scan},
                  filename=os.path.join(caldir, '{}_Dec.fits'.format(i)),
                  src_ra=src_ra, src_dec=src_dec, srcname='DummyCal')
        timedelta = times[-1] + 1


class TestCalibration(object):
    @classmethod
    def setup_class(klass):
        import os
        global DEBUG_MODE
        DEBUG_MODE = True

        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, 'data')

        klass.config_file = \
            os.path.abspath(os.path.join(klass.datadir, "calibrators.ini"))

        klass.config = read_config(klass.config_file)
        klass.caldir = os.path.join(klass.datadir, 'sim', 'calibration')
        if not os.path.exists(klass.caldir):
            print('Fake calibrators: DummyCal, 1 Jy.')
            mkdir_p(klass.caldir)
            sim_calibrators(5, klass.caldir)

        klass.scan_list = \
            list_scans(klass.caldir, ['./'])

        klass.scan_list.sort()
        caltable = CalibratorTable()
        caltable.from_scans(klass.scan_list)
        caltable.update()
        klass.calfile = os.path.join(klass.curdir, 'test_calibration.hdf5')
        caltable.write(klass.calfile, overwrite=True)

    def test_0_prepare(self):
        pass

    def test_check_class(self):
        caltable = CalibratorTable()
        caltable.from_scans(self.scan_list)
        with pytest.warns(UserWarning):
            caltable.check_up_to_date()

    def test_check_class_from_file(self):
        caltable = CalibratorTable.read(self.calfile)
        assert caltable.check_up_to_date()

    def test_bad_file(self):
        caltable = CalibratorTable()
        with pytest.warns(UserWarning) as record:
            caltable.from_scans([os.path.join(self.config['datadir'],
                                              'calibrators', 'summary.fits')])
        assert "Missing key" in record[0].message.args[0]

        with pytest.warns(UserWarning) as record:
            caltable.from_scans([os.path.join(self.config['datadir'],
                                              'calibrators', 'bubu.fits')])
        assert "Error while processing" in record[0].message.args[0]

    def test_calibration_counts(self):
        """Simple calibration from scans."""

        caltable = CalibratorTable.read(self.calfile)
        assert np.all(
            np.abs(caltable['Counts'] - 100.) < 3 * caltable['Counts Err'])

    def test_calibration_width(self):
        """Simple calibration from scans."""

        caltable = CalibratorTable.read(self.calfile)
        assert np.all(
            np.abs(caltable['Width'] - 3/60.) < 3 * caltable['Width Err'])

        beam, beam_err = caltable.beam_width()
        assert np.all(beam - np.radians(3/60) < 3 * beam_err)

    def test_calibration_plot(self):
        """Simple calibration from scans."""

        caltable = CalibratorTable.read(self.calfile)
        caltable.show()

    @classmethod
    def teardown_class(klass):
        """Clean up the mess."""
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
