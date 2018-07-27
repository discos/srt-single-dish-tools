from srttools.convert import convert_to_complete_fitszilla, main_convert
from srttools.scan import Scan
import numpy as np
import os
import pytest
import subprocess as sp
import shutil
import glob
from astropy.io import fits
from srttools.io import locations, mkdir_p
from srttools.simulate import simulate_map, sim_position_switching, \
    DEFAULT_CAL_OFFSET

try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

CI = os.getenv('CI') is not None
CI_MPL = CI and HAS_MPL


class Test1_Scan(object):
    @classmethod
    def setup_class(klass):
        import os

        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, 'data')

        klass.fname = \
            os.path.abspath(
                os.path.join(klass.datadir,
                             'srt_data_tp_multif.fits'))

        klass.skydip = \
            os.path.abspath(
                os.path.join(klass.datadir,
                             'gauss_skydip'))
        klass.example = \
            os.path.abspath(
                os.path.join(klass.datadir,
                             'example_polar'))
        klass.onoff = \
            os.path.abspath(
                os.path.join(klass.datadir,
                             'onoff_xarcos'))
        klass.nodding = \
            os.path.abspath(
                os.path.join(klass.datadir,
                             'nodding_xarcos'))
        klass.outdir = os.path.join('sim')
        klass.emptydir = os.path.join('sim', 'test_sdfits')

        klass.pswdir = os.path.join('sim', 'test_psw')
        for d in [klass.emptydir, klass.pswdir]:
            mkdir_p(d)
        sim_position_switching(klass.pswdir, nbin=1024)
        simulate_map(width_ra=2, width_dec=2., outdir=klass.emptydir)

    def test_converter_basic(self):
        convert_to_complete_fitszilla(self.fname, 'converted')
        os.unlink('converted.fits')

    def test_installed(self):
        sp.check_call('SDTconvert -h'.split(' '))

    def test_conversion(self):
        convert_to_complete_fitszilla(self.fname, 'converted')
        scan0 = Scan(self.fname, norefilt=False)
        scan1 = Scan('converted.fits', norefilt=False)
        for col in ['ra', 'el', 'az', 'dec']:
            assert np.allclose(scan0[col], scan1[col])
        os.unlink('converted.fits')

    def test_conversion_same_name_fails(self):
        with pytest.raises(ValueError):
            convert_to_complete_fitszilla(self.fname, self.fname)

    def test_main(self):
        main_convert([self.fname, '-f', 'fitsmod'])
        assert os.path.exists(self.fname.replace('.fits',
                                                 '_fitsmod.fits'))
        os.unlink(self.fname.replace('.fits', '_fitsmod.fits'))

    def test_main_dir(self):
        main_convert([self.skydip, '-f', 'fitsmod'])
        newfile = os.path.join(self.skydip,
                               'skydip_mod_fitsmod.fits')
        assert os.path.exists(newfile)
        os.unlink(newfile)

    def test_main_garbage_format(self):
        with pytest.warns(UserWarning):
            main_convert([self.fname, '-f', 'weruoiq'])

        assert not os.path.exists(self.fname.replace('.fits',
                                                     '_weruoiq.fits'))

    def test_main_nondir_mbfits(self):
        with pytest.raises(ValueError) as excinfo:
            main_convert([self.fname, '-f', 'mbfits'])

        assert "Input for MBFITS conversion must be " in str(excinfo)

    @pytest.mark.skipif('CI_MPL')
    def test_main_mbfitsw(self):
        main_convert([self.skydip, '-f', 'mbfitsw', '--test'])
        newfiles = glob.glob(self.skydip + '*KKG*.fits')
        assert len(newfiles) > 0
        # test that a new conversion does not make this fail
        newdir = main_convert([self.skydip, '-f', 'mbfitsw', '--test'])[0]

        shutil.rmtree(newdir)
        with fits.open(newfiles[0]) as hdul:
            header = hdul['SCAN-MBFITS'].header
            assert header['SCANTYPE'] == 'SKYDIP'
        for fname in newfiles:
            os.unlink(fname)

    @pytest.mark.skipif('CI_MPL')
    def test_main_mbfitsw_polar(self):
        newdir = main_convert([self.example, '-f', 'mbfitsw', '--test'])[0]
        newfiles = glob.glob(self.example + '*CCB*.fits')
        assert len(newfiles) > 0
        shutil.rmtree(newdir)
        with fits.open(newfiles[0]) as hdul:
            header = hdul['SCAN-MBFITS'].header
            assert header['SCANTYPE'] == 'MAP'
        for fname in newfiles:
            os.unlink(fname)

    @pytest.mark.skipif('CI_MPL')
    def test_main_mbfits(self):
        newdir = main_convert([self.skydip, '-f', 'mbfits', '--test'])[0]
        assert os.path.exists(newdir)
        assert os.path.isdir(newdir)
        assert os.path.exists(os.path.join(newdir, 'GROUPING.fits'))
        scanfile = os.path.join(newdir, 'SCAN.fits')
        assert os.path.exists(scanfile)
        with fits.open(scanfile) as hdul:
            header = hdul[1].header
            assert header['SCANTYPE'] == 'SKYDIP'
        shutil.rmtree(newdir)

    def test_main_classfits_onoff(self):
        newdir = main_convert([self.onoff, '-f', 'classfits', '--test'])[0]
        assert os.path.exists(newdir)
        # test that a new conversion does not make this fail
        newdir = main_convert([self.onoff, '-f', 'classfits', '--test'])[0]
        assert os.path.isdir(newdir)

    def test_main_classfits_sim(self):
        newdir = main_convert([self.pswdir, '-f', 'classfits', '--test'])[0]
        assert os.path.exists(newdir)
        assert os.path.isdir(newdir)
        probe_all = os.path.join(newdir, 'test_psw_all_feed0.fits')
        probe_cal = os.path.join(newdir, 'test_psw_cal_feed0.fits')
        probe_psw = os.path.join(newdir, 'test_psw_onoff_feed0.fits')
        with fits.open(probe_all) as hdul:
            good = (hdul[1].data['SIGNAL'] == 0) & (hdul[1].data['CAL_IS_ON'] == 0)
            off_spec = hdul[1].data['SPECTRUM'][good][0]
            good = (hdul[1].data['SIGNAL'] == 1)
            on_spec = hdul[1].data['SPECTRUM'][good][0]
        with fits.open(probe_cal) as hdul:
            cal_spec = hdul[1].data['SPECTRUM'][0]
        with fits.open(probe_psw) as hdul:
            onoff_spec = hdul[1].data['SPECTRUM'][0]
        assert np.isclose(np.max(on_spec - off_spec), 100, atol=0.1)
        idx = np.argmax(on_spec)
        max_onoff = onoff_spec[idx]
        ref_off = off_spec[idx]
        ref_cal = cal_spec[idx]

        assert np.isclose(max_onoff * ref_off, 100, atol=0.1)
        assert np.isclose(ref_cal, 100 / DEFAULT_CAL_OFFSET, atol=0.1)

    def test_main_classfits_nodding(self):
        newdir = main_convert([self.nodding, '-f', 'classfits', '--test'])[0]
        assert os.path.exists(newdir)
        assert os.path.isdir(newdir)
        shutil.rmtree(self.nodding + '_classfits')

    def test_main_sdfits_skydip(self):
        newdir = main_convert([self.skydip, '-f', 'sdfits', '--test'])[0]
        assert os.path.exists(newdir)
        assert os.path.isdir(newdir)
        newfiles = glob.glob(os.path.join(newdir, '*.fits'))

        with fits.open(newfiles[0]) as hdul:
            header = hdul['SINGLE DISH'].header
            assert np.isclose(header['OBSGEO-X'], locations['srt'].x.value)
        shutil.rmtree(self.skydip + '_sdfits')

    def test_main_sdfits_nodding(self):
        newdir = main_convert([self.nodding, '-f', 'sdfits', '--test'])[0]
        assert os.path.exists(newdir)
        # test that it doesn't fail when the directory is already present
        newdir = main_convert([self.nodding, '-f', 'sdfits', '--test'])[0]
        assert os.path.isdir(newdir)

        newfiles = glob.glob(os.path.join(newdir, '*.fits'))

        with fits.open(newfiles[0]) as hdul:
            header = hdul['SINGLE DISH'].header
            assert np.isclose(header['OBSGEO-X'], locations['srt'].x.value)

        shutil.rmtree(self.nodding + '_sdfits')

    def test_main_sdfits_tpmap(self):
        newdir = main_convert([self.emptydir, '-f', 'sdfits', '--test'])[0]
        assert os.path.exists(newdir)
        # test that it doesn't fail when the directory is already present
        newdir = main_convert([self.emptydir, '-f', 'sdfits', '--test'])[0]
        assert os.path.isdir(newdir)

        newfiles = glob.glob(os.path.join(newdir, '*.fits'))

        with fits.open(newfiles[0]) as hdul:
            header = hdul['SINGLE DISH'].header
            assert np.isclose(header['OBSGEO-X'], locations['srt'].x.value)

        shutil.rmtree(self.emptydir + '_sdfits')
