from __future__ import (absolute_import, division,
                        print_function)
import numpy as np
try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

from srttools.destripe import basket_weaving, destripe_wrapper
np.random.seed(450720239)


class TestScanSet(object):
    @classmethod
    def setup_class(klass):
        x = np.linspace(-1, 1, 101)
        y = np.linspace(-1, 1, 101)
        xx, yy = np.meshgrid(x, y)

        zz = np.exp(- (xx ** 2 + yy ** 2) / 0.01)

        klass.img = zz
        klass.img_hor = 0
        klass.img_ver = 0

        example_scan = np.cumsum(np.random.normal(0, 0.01, len(zz[0])))

        for i in range(6):
            offset = np.random.choice(example_scan)
            zz_hor = np.array(
                [z + offset +
                 np.cumsum(np.random.normal(0, 0.01, len(z))) for z in zz])
            offset = np.random.choice(example_scan)
            zz_ver = np.array(
                [z + offset +
                 np.cumsum(np.random.normal(0, 0.01, len(z))) for z in zz.T]).T

            if i != 0:  # One less integration over x
                klass.img_hor += zz_hor
            klass.img_ver += zz_ver

        klass.img_hor /= 5
        klass.img_ver /= 6
        klass.img_hor[30:34, :] = 0
        klass.expo_hor = np.ones_like(klass.img_hor) * 5
        klass.expo_ver = np.ones_like(klass.img_ver) * 6
        klass.expo_hor[30:34, :] = 0

    def test_no_expo(self):
        img_clean = basket_weaving(self.img_hor, self.img_ver)

        diff_img = (self.img - img_clean)
        mean = np.mean(diff_img)
        std = np.std(diff_img)

        assert np.all(np.abs(diff_img - mean) < std * 5)

    def test_wrapper(self):
        img_clean = basket_weaving(self.img_hor, self.img_ver)
        img_clean_w = destripe_wrapper(self.img_hor, self.img_ver)

        assert np.all(img_clean == img_clean_w)

    def test_expo(self):
        img_clean = basket_weaving(self.img_hor, self.img_ver,
                                   expo_hor=self.expo_hor,
                                   expo_ver=self.expo_ver)

        diff_img = (self.img - img_clean)
        mean = np.mean(diff_img)
        std = np.std(diff_img)

        assert np.all(np.abs(diff_img - mean) < std * 5)
