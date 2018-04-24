from srttools.convert import convert_to_complete_fitszilla, main_convert
from srttools.scan import Scan
import numpy as np
import os
import pytest
import subprocess as sp
import shutil


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

    def test_main_mbfits(self):
        main_convert([self.skydip, '-f', 'mbfits', '--test'])
        newdir = self.skydip + '_mbfits'
        assert os.path.exists(newdir)
        assert os.path.isdir(newdir)
        assert os.path.exists(os.path.join(newdir, 'GROUPING.fits'))
        assert os.path.exists(os.path.join(newdir, 'SCAN.fits'))
        shutil.rmtree(self.skydip + '_mbfits')

