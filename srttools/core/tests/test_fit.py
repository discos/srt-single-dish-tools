# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function)
from srttools.core.fit import fit_baseline_plus_bell, purge_outliers, align
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(1231636)


def _test_shape(x):
    return 100 * np.exp(-(x - 50) ** 2 / 3)


class TestFit(object):
    @classmethod
    def setup_class(cls):
        cls.series = np.random.normal(0, 0.1, 1000)
        cls.t = np.arange(0, len(cls.series)/10, 0.1)

    def test_outliers1(self):
        """Test that outlier detection works."""
        series = np.copy(self.series)
        series[10] = 2
        series2 = purge_outliers(series)
        assert np.all(series2[:10] == series[:10])
        assert np.all(series2[11:] == series[11:])
        np.testing.assert_almost_equal(series2[10],
                                       (series[9] + series[11]) / 2)

    def test_outliers2(self):
        """Test that outlier detection works."""
        series = np.copy(self.series)
        series[10] = -2
        series2 = purge_outliers(series)
        assert np.all(series2[:10] == series[:10])
        assert np.all(series2[11:] == series[11:])
        np.testing.assert_almost_equal(series2[10],
                                       (series[9] + series[11]) / 2)

    def test_outliers3(self):
        """Test that outlier detection works."""
        series = np.copy(self.series)
        series[10] = 20
        series[11] = 20
        series2 = purge_outliers(series)

        assert np.all(series2[:10] == series[:10])
        assert np.all(series2[12:] == series[12:])

        assert np.all((series2[10:12] > series[9]) &
                      (series2[10:12] < series[12]))

    def test_outliers_bell(self):
        """Test that outlier detection works."""
        series = np.copy(self.series) + _test_shape(self.t) / 10
        series[10] = 2
        series2 = purge_outliers(series)
        assert np.all(series2[:10] == series[:10])
        assert np.all(series2[11:] == series[11:])
        np.testing.assert_almost_equal(series2[10],
                                       (series[9] + series[11]) / 2)

    def test_outliers_bell_larger(self):
        """Test that outlier detection works."""
        series = np.copy(self.series) + _test_shape(self.t)
        series[10] = 2
        series2 = purge_outliers(series)
        assert np.all(series2[:10] == series[:10])
        assert np.all(series2[11:] == series[11:])
        np.testing.assert_almost_equal(series2[10],
                                       (series[9] + series[11]) / 2)

    def test_fit_baseline_plus_bell(self):
        """Test that the fit procedure works."""

        x = np.arange(0, len(self.series)) * 0.1
        y = np.copy(self.series) + _test_shape(x) + x * 6 + 20

        model, _ = fit_baseline_plus_bell(x, y, ye=10, kind='gauss')

        np.testing.assert_almost_equal(model.mean_1, 50., 1)
        np.testing.assert_almost_equal(model.slope_0, 6., 1)
        assert np.abs(model.intercept_0 - 20.) < 2

    def test_minimize_align(self):
        """Test that the minimization of the alignment works."""

        x1 = np.arange(0, 100, 0.1)
        y1 = np.random.poisson(100, len(x1)) + _test_shape(x1)
        x2 = np.arange(0.02, 100, 0.1)
        y2 = np.random.poisson(100, len(x2)) + _test_shape(x2)
        x3 = np.arange(0.053, 98.34, 0.1)
        y3 = np.random.poisson(100, len(x3)) + _test_shape(x3)

        xs = [x1, x2, x3]
        ys = [y1, y2, y3]

        qs = [0, -60, 60]
        ms = [0, 0.3, -0.8]

        for ix, x in enumerate(xs):
            ys[ix] = ys[ix] + qs[ix] + ms[ix] * xs[ix]

        qs, ms = align(xs, ys)

        np.testing.assert_allclose(qs, [-60, 60], atol=3)
        np.testing.assert_allclose(ms, [0.3, -0.8], atol=0.05)
