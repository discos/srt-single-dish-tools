# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function)
from srttools.core.fit import fit_baseline_plus_bell, purge_outliers
import numpy as np

np.random.seed(1231636)


def _test_shape(x):
    return 1000 * np.exp(-(x - 50) ** 2 / 3)


class TestFit(object):
    @classmethod
    def setup_class(cls):
        cls.series = np.random.normal(0, 0.1, 1000)

    def test_outliers1(self):
        """Test that outlier detection works."""
        series = np.copy(self.series)
        series[10] = 2
        series2 = purge_outliers(series)
        assert np.all(series2[:10] == series[:10])
        assert np.all(series2[11:] == series[11:])
        np.testing.assert_almost_equal(series2[10], (series[9] + series[11]) / 2)

    def test_outliers2(self):
        """Test that outlier detection works."""
        series = np.copy(self.series)
        series[10] = -2
        series2 = purge_outliers(series)
        assert np.all(series2[:10] == series[:10])
        assert np.all(series2[11:] == series[11:])
        np.testing.assert_almost_equal(series2[10], (series[9] + series[11]) / 2)

    def test_outliers3(self):
        """Test that outlier detection works."""
        series = np.copy(self.series)
        series[10] = 20
        series[11] = 20
        series2 = purge_outliers(series)

        assert np.all(series2 == series)

    def test_fit_baseline_plus_bell(self):
        """Test that the fit procedure works."""

        x = np.arange(0, len(self.series)) * 0.1
        y = np.copy(self.series) + _test_shape(x) + x * 6 + 20

        model, _ = fit_baseline_plus_bell(x, y, ye=10, kind='gauss')

        np.testing.assert_almost_equal(model.mean_1, 50., 1)
        np.testing.assert_almost_equal(model.slope_0, 6., 1)
        assert np.abs(model.intercept_0 - 20.) < 2
