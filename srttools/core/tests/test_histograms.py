from __future__ import division, print_function
from srttools.core import histograms as hist
import unittest
import numpy as np

class Test_Hist(unittest.TestCase):
    @classmethod
    def setup_class(klass):
        import os
        global DEBUG_MODE
        DEBUG_MODE = True


    def test_hist_numbers(self):
        a = np.random.poisson(100, 10000)
        b = np.random.poisson(100, 10000)
        bins = np.linspace(50, 150, 51)
        hnum, xbnum, ybnum = np.histogram2d(a, b, bins=(bins, bins))
        hh, xbh, ybh = hist.histogram2d(a, b, bins=(bins, bins))
        np.testing.assert_equal(hnum, hh)



