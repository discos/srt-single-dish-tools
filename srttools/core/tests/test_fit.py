# -*- coding: utf-8 -*-

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
from ..fit import fit_baseline_plus_bell, align, purge_outliers
import numpy as np
import matplotlib.pyplot as plt
import unittest
from astropy.table import Table
from ..imager import ScanSet
import os
import glob

def _test_shape(x):
    return 1000 * np.exp(-(x - 50) ** 2 / 3)


def test_outliers():
    """Test that outlier detection works."""
    series = np.random.normal(0, 0.1, 100)
    series[10] = 2
    series2 = purge_outliers(series)
    assert np.all(series2[:10] == series[:10])
    assert np.all(series2[11:] == series[11:])
    np.testing.assert_almost_equal(series2[10], (series[9] + series[11]) / 2)

# def test_fit_baseline_plus_bell():
#     """Test that the fit procedure works."""
#     import matplotlib.pyplot as plt
#
#     x = np.arange(0, 100, 0.1)
#     y = np.random.poisson(1000, len(x)) + _test_shape(x) + x * 6 + 20
#
#     model, _ = fit_baseline_plus_bell(x, y, ye=10, kind='gauss')
#
#     plt.figure('After fit')
#     plt.errorbar(x, y, yerr=10)
#
#     plt.plot(x, model(x))
#
#     print('Amplitude: {}'.format(model.amplitude_1.value))
#     print('Stddev: {}'.format(model.stddev_1.value))
#     print('Mean: {}'.format(model.mean_1.value))
#     plt.show()
#
#
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
