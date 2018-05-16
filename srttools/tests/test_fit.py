# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function)
from srttools.fit import fit_baseline_plus_bell, purge_outliers, align
from srttools.fit import baseline_rough, ref_mad, ref_std, _rolling_window
from srttools.fit import linear_fit, offset_fit, detrend_spectroscopic_data

import numpy as np
import pytest

np.random.seed(1231636)


def _test_shape(x):
    return 100 * np.exp(-(x - 50) ** 2 / 3)


def _setup_spectra(nx, ny):
    spectrum = np.vstack([np.linspace(0 + i, 2 + i, nx)
                          for i in np.linspace(0., 4, ny)]) + 1
    spectrum_noisy = np.random.normal(spectrum, 0.00001)
    x = np.arange(spectrum.shape[0])
    return x, spectrum_noisy


list_of_par_pairs = [(2 * n + 1, 2 * m + 1)
                     for (n, m) in zip(np.random.randint(1, 50, 10),
                                       np.random.randint(1, 50, 10))]


class TestStuff(object):
    @pytest.mark.parametrize('nx,ny', list_of_par_pairs)
    def test_detrend_spectroscopic_data(self, nx, ny):
        x, spectrum = _setup_spectra(nx, ny)
        detr, _ = detrend_spectroscopic_data(x, spectrum,
                                             kind='rough')
        assert np.allclose(detr, 0., atol=1e-3)

    @pytest.mark.parametrize('nx,ny', list_of_par_pairs)
    def test_detrend_spectroscopic_data_als(self, nx, ny):
        x, spectrum = _setup_spectra(nx, ny)
        detr, _ = detrend_spectroscopic_data(x, spectrum, kind='als',
                                             outlier_purging=False)
        assert np.allclose(detr, 0., atol=1e-2)

    @pytest.mark.parametrize('nx,ny', list_of_par_pairs)
    def test_detrend_spectroscopic_data_garbage(self, nx, ny):
        x, spectrum = _setup_spectra(nx, ny)
        detr, _ = detrend_spectroscopic_data(x, spectrum,
                                             kind='blabla')
        assert np.all(detr == spectrum)


class TestFit(object):
    @classmethod
    def setup_class(cls):
        cls.series = np.random.normal(0, 0.1, 1000)
        cls.t = np.arange(0, len(cls.series)/10, 0.1)

    def test_outliers1(self):
        """Test that outlier detection works."""
        series = np.copy(self.series)
        series[10] = 2
        with pytest.warns(UserWarning) as record:
            series2 = purge_outliers(series)
        assert "Found 1 outliers" in record[0].message.args[0]
        assert np.all(series2[:10] == series[:10])
        assert np.all(series2[11:] == series[11:])
        np.testing.assert_almost_equal(series2[10],
                                       (series[9] + series[11]) / 2)

    def test_outliers2(self):
        """Test that outlier detection works."""
        series = np.copy(self.series)
        series[10] = -2
        with pytest.warns(UserWarning) as record:
            series2 = purge_outliers(series)
        assert "Found 1 outliers" in record[0].message.args[0]
        assert np.all(series2[:10] == series[:10])
        assert np.all(series2[11:] == series[11:])
        np.testing.assert_almost_equal(series2[10],
                                       (series[9] + series[11]) / 2)

    def test_outliers3(self):
        """Test that outlier detection works."""
        series = np.copy(self.series)
        series[10] = 20
        series[11] = 20
        with pytest.warns(UserWarning) as record:
            series2 = purge_outliers(series)
        assert "Found 2 outliers" in record[0].message.args[0]

        assert np.all(series2[:10] == series[:10])
        assert np.all(series2[12:] == series[12:])

        lower = np.min([series[9], series[12]])
        upper = np.max([series[9], series[12]])
        assert np.all((series2[10:12] > lower) &
                      (series2[10:12] < upper))

    def test_outliers_bell(self):
        """Test that outlier detection works."""
        series = np.copy(self.series) + _test_shape(self.t) / 10
        series[10] = 2
        with pytest.warns(UserWarning) as record:
            series2 = purge_outliers(series)
        assert "Found 1 outliers" in record[0].message.args[0]
        assert np.all(series2[:10] == series[:10])
        assert np.all(series2[11:] == series[11:])
        np.testing.assert_almost_equal(series2[10],
                                       (series[9] + series[11]) / 2)

    def test_outliers_bell_larger(self):
        """Test that outlier detection works."""
        series = np.copy(self.series) + _test_shape(self.t)
        series[10] = 2
        with pytest.warns(UserWarning) as record:
            series2 = purge_outliers(series)
        assert "Found 1 outliers" in record[0].message.args[0]
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

    def test_fit_baseline_plus_bell_lorentz(self):
        """Test that the fit procedure works."""

        x = np.arange(0, len(self.series)) * 0.1
        y = np.copy(self.series) + _test_shape(x) + x * 6 + 20

        model, _ = fit_baseline_plus_bell(x, y, ye=10, kind='lorentz')

        np.testing.assert_almost_equal(model.x_0_1, 50., 1)

    def test_fit_baseline_plus_bell_invalid(self):
        """Test that the fit procedure works."""

        x = np.arange(0, len(self.series)) * 0.1
        y = np.copy(self.series) + _test_shape(x) + x * 6 + 20
        with pytest.raises(ValueError) as excinfo:
            model, _ = fit_baseline_plus_bell(x, y, ye=10, kind='zxcdf')
        assert 'kind has to be one of: gauss, lorentz' in str(excinfo)

    def test_fit_baseline_rough_return_baseline(self):
        """Test that the fit procedure works."""

        x = np.arange(0, len(self.series)) * 0.1
        y = np.copy(self.series) + _test_shape(x) + x * 6 + 20
        a, b = baseline_rough(x, y, return_baseline=True)

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

    def test_rolling_window_invalid(self):
        with pytest.raises(Exception):
            _rolling_window(1, 5)

    def test_ref_std_small_array(self):
        array = np.zeros(10)
        window = 10
        a = ref_std(array, window)
        assert a == 0

    def test_ref_mad_small_array(self):
        array = np.zeros(10)
        window = 10
        a = ref_mad(array, window)
        assert a == 0

    def test_linear_fit_err_not_implemented(self):
        x = np.arange(10)
        y = np.zeros(10)
        with pytest.warns(UserWarning) as record:
            _, err = linear_fit(x, y, [0, 0], return_err=True)
            assert "return_err not implemented" in record[0].message.args[0]
        assert err is None

    def test_offset_fit_err_not_implemented(self):
        x = np.arange(10)
        y = np.zeros(10)
        with pytest.warns(UserWarning) as record:
            _, err = offset_fit(x, y, 0, return_err=True)
            assert "return_err not implemented" in record[0].message.args[0]
        assert err is None
