# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function)
from srttools.core.fit import fit_baseline_plus_bell, align, purge_outliers
from srttools.core.global_fit import fit_full_image, display_intermediate
import numpy as np
import matplotlib.pyplot as plt
import unittest
from astropy.table import Table
from ..imager import ScanSet
import os
import glob

np.random.seed(1231636)

def test_outliers1():
    """Test that outlier detection works."""
    series = np.random.normal(0, 0.1, 100)
    series[10] = 2
    series2 = purge_outliers(series)
    assert np.all(series2[:10] == series[:10])
    assert np.all(series2[11:] == series[11:])
    np.testing.assert_almost_equal(series2[10], (series[9] + series[11]) / 2)

def test_outliers2():
    """Test that outlier detection works."""
    series = np.random.normal(0, 0.1, 100)
    series[10] = -2
    series2 = purge_outliers(series)
    assert np.all(series2[:10] == series[:10])
    assert np.all(series2[11:] == series[11:])
    np.testing.assert_almost_equal(series2[10], (series[9] + series[11]) / 2)

def test_outliers3():
    """Test that outlier detection works."""
    series = np.random.normal(0, 0.1, 100)
    series[10] = 20
    series[11] = 20
    series2 = purge_outliers(series)

    assert np.all(series2 == series)

def _test_shape(x):
    return 1000 * np.exp(-(x - 50) ** 2 / 3)

def test_fit_baseline_plus_bell():
    """Test that the fit procedure works."""

    x = np.arange(0, 100, 0.1)
    y = np.random.normal(0, 1, len(x)) + _test_shape(x) + x * 6 + 20

    model, _ = fit_baseline_plus_bell(x, y, ye=10, kind='gauss')

    np.testing.assert_almost_equal(model.mean_1, 50., 1)
    np.testing.assert_almost_equal(model.slope_0, 6., 1)
    assert np.abs(model.intercept_0 - 20.) < 2

# def test_minimize_align():
#     """Test that the minimization of the alignment works."""
#
#     x1 = np.arange(0, 100, 0.1)
#     y1 = np.random.poisson(100, len(x1)) + _test_shape(x1)
#     x2 = np.arange(0.02, 100, 0.1)
#     y2 = np.random.poisson(100, len(x2)) + _test_shape(x2)
#     x3 = np.arange(0.053, 98.34, 0.1)
#     y3 = np.random.poisson(100, len(x3)) + _test_shape(x3)
#
#     xs = [x1, x2, x3]
#     ys = [y1, y2, y3]
#
#     qs = [20, -60, 60]
#     ms = [0, 0.3, -0.8]
#
#     plt.figure('Original')
#     for ix, x in enumerate(xs):
#         ys[ix] = ys[ix] + qs[ix] + ms[ix] * xs[ix]
#
#         plt.plot(xs[ix], ys[ix], label=ix)
#
#     plt.legend()
#
#     qs, ms = align(xs, ys)
#
#     plt.figure('After fit')
#     for ix, x in enumerate(xs):
#         if ix == 0:
#             plt.plot(xs[0], ys[0], label=ix)
#         else:
#             ys[ix] = ys[ix] - (qs[ix-1] + ms[ix-1] * xs[ix])
#
#             plt.plot(xs[ix], ys[ix], label=ix)
#
#     plt.legend()
#     plt.show()
