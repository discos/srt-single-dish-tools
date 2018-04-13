from srttools.convert import convert_to_complete_fitszilla, converter_main
from srttools.scan import Scan
import numpy as np
import os
import pytest


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

    def test_converter_basic(self):
        convert_to_complete_fitszilla(self.fname, 'converted')
        os.unlink('converted.fits')

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
        converter_main([self.fname, '-f', 'fitsmod'])
        assert os.path.exists(self.fname.replace('.fits',
                                                 '_fitsmod.fits'))
        os.unlink(self.fname.replace('.fits', '_fitsmod.fits'))

    def test_main_garbage_format(self):
        with pytest.warns(UserWarning):
            converter_main([self.fname, '-f', 'weruoiq'])

        assert not os.path.exists(self.fname.replace('.fits',
                                                     '_weruoiq.fits'))
