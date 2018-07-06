from srttools.convert import convert_to_complete_fitszilla, main_convert
from srttools.scan import Scan
import numpy as np
import os
import pytest
import subprocess as sp
import shutil
import glob
from astropy.io import fits

try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


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

    @pytest.mark.skipif('HAS_MPL')
    def test_main_mbfitsw(self):
        main_convert([self.skydip, '-f', 'mbfitsw', '--test'])
        newfiles = glob.glob(self.skydip + '*KKG*.fits')
        assert len(newfiles) > 0
        shutil.rmtree(self.skydip + '_mbfitsw')
        with fits.open(newfiles[0]) as hdul:
            header = hdul['SCAN-MBFITS'].header
            assert header['SCANTYPE'] == 'SKYDIP'
        for fname in newfiles:
            os.unlink(fname)

    @pytest.mark.skipif('HAS_MPL')
    def test_main_mbfitsw_polar(self):
        main_convert([self.example, '-f', 'mbfitsw', '--test'])
        newfiles = glob.glob(self.example + '*CCB*.fits')
        assert len(newfiles) > 0
        shutil.rmtree(self.example + '_mbfitsw')
        with fits.open(newfiles[0]) as hdul:
            header = hdul['SCAN-MBFITS'].header
            assert header['SCANTYPE'] == 'MAP'
        for fname in newfiles:
            os.unlink(fname)

    @pytest.mark.skipif('HAS_MPL')
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
        shutil.rmtree(self.skydip + '_mbfits')

    def test_main_classfits_onoff(self):
        newdir = main_convert([self.onoff, '-f', 'classfits', '--test'])[0]
        assert os.path.exists(newdir)
        assert os.path.isdir(newdir)
        shutil.rmtree(self.onoff + '_classfits')

    def test_main_classfits_nodding(self):
        newdir = main_convert([self.nodding, '-f', 'classfits', '--test'])[0]
        assert os.path.exists(newdir)
        assert os.path.isdir(newdir)
        shutil.rmtree(self.nodding + '_classfits')
