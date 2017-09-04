from __future__ import division, print_function
from srttools import histograms as hist
import numpy as np

np.random.seed(1742956)


class TestHist(object):
    @classmethod
    def setup_class(klass):
        klass.N = 10000
        klass.a = np.random.poisson(100, klass.N)
        klass.b = np.random.poisson(100, klass.N)
        klass.bins = np.linspace(50, 150, 51)

    def test_hist_numbers(self):
        hnum, xbnum, ybnum = np.histogram2d(self.a,
                                            self.b,
                                            bins=(self.bins, self.bins))
        hh, xbh, ybh = hist.histogram2d(self.a,
                                        self.b,
                                        bins=(self.bins, self.bins))
        np.testing.assert_equal(hnum, hh)

    def test_hist_numbers_normed_and_weights(self):
        w = np.random.uniform(1, 0.01, self.N)
        hnum, xbnum, ybnum = np.histogram2d(self.a,
                                            self.b,
                                            bins=(self.bins, self.bins),
                                            weights=w, normed=True)
        hh, xbh, ybh = hist.histogram2d(self.a,
                                        self.b,
                                        bins=(self.bins, self.bins),
                                        weights=w, normed=True)
        np.testing.assert_equal(hnum, hh)
