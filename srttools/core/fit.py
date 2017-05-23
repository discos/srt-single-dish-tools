"""Useful fitting functions."""
from __future__ import (absolute_import, division,
                        print_function)
from scipy.optimize import curve_fit
import numpy as np
import traceback
import warnings


def _rolling_window(a, window):
    """A smart rolling window.

    Found at http://www.rigtorp.se/2011/01/01/rolling-statistics-numpy.html
    """
    try:
        shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
        strides = a.strides + (a.strides[-1],)
        return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)
    except Exception:
        warnings.warn(traceback.traceback.format_exc())
        raise


def ref_std(array, window):
    """Minimum standard deviation along an array."""

    return np.min(np.std(_rolling_window(array, window), 1))


def linear_fun(x, q, m):
    """A linear function."""
    return m * x + q


def linear_fit(time, lc, start_pars, return_err=False):
    """A linear fit with any set of data. Return the parameters."""
    par, pcov = curve_fit(linear_fun, time, lc, start_pars,
                          maxfev=6000)
    if return_err:
        pass
    else:
        return par


def offset(x, off):
    """An offset."""
    return off


def offset_fit(time, lc, offset_start=0, return_err=False):
    """A linear fit with any set of data. Return the parameters."""
    par, pcov = curve_fit(offset, time, lc, [offset_start],
                          maxfev=6000)
    if return_err:
        pass
    else:
        return par[0]


def baseline_rough(time, lc, start_pars=None, return_baseline=False):
    """Rough function to subtract the baseline."""
    if start_pars is None:
        m0 = 0
        q0 = min(lc)
        start_pars = [q0, m0]

    nbin = len(time)
    #    bins = np.arange(nbin, dtype=int)
    lc = lc.copy()
    time = time.copy()

    total_trend = 0

    local_std = ref_std(lc, np.max([nbin // 20, 20]))

    for percentage in [0.8, 0.15]:
        time_to_fit = time[1:-1]
        lc_to_fit = lc[1:-1]

        sorted_els = np.argsort(lc_to_fit)
        # Select the lowest half elements
        good = sorted_els[: int(nbin * percentage)]
        #    good = np.logical_or(bins <= nbin / 4, bins >= nbin / 4 * 3)

        print(np.std(lc_to_fit[good]), local_std)
        if np.std(lc_to_fit[good]) < 2 * local_std:
            good = np.ones(len(lc_to_fit), dtype=bool)

        time_filt = time_to_fit[good]
        lc_filt = lc_to_fit[good]
        back_in_order = np.argsort(time_filt)
        lc_filt = lc_filt[back_in_order]
        time_filt = time_filt[back_in_order]

        par = linear_fit(time_filt, lc_filt, start_pars)

        lc -= linear_fun(time, *par)
        total_trend += linear_fun(time, *par)

    if return_baseline:
        return lc, total_trend
    else:
        return lc


def purge_outliers(y):
    y = y.copy()
    idxs = np.arange(len(y))

    min_diff = ref_std(y, np.max([len(y) // 15, 30]))
    diffs = np.diff(y)
    diffs = np.append([0], diffs)
    diffs_before = np.array(diffs)[:-1]
    diffs_after = np.array(diffs)[1:]
    sign_rule = np.sign(diffs_before) != np.sign(diffs_after)
    outliers = (np.abs(diffs_before) > 10 * min_diff) & \
               (np.abs(diffs_after) > 10 * min_diff) & sign_rule
    outlier_idxs = idxs[:-1][outliers]
    for i in outlier_idxs:
        y[i] = (y[i - 1] + y[i + 1]) / 2

    if len(outlier_idxs) > 0:
        warnings.warn("Found {} outliers".format(len(outlier_idxs)))

    return y


def _als(y, lam, p, niter=10):
    """Baseline Correction with Asymmetric Least Squares Smoothing.

    Modifications to the routine from Eilers & Boelens 2005
    https://www.researchgate.net/publication/
        228961729_Technical_Report_Baseline_Correction_with_
        Asymmetric_Least_Squares_Smoothing
    The Python translation is partly from
    http://stackoverflow.com/questions/29156532/
        python-baseline-correction-library
    """
    from scipy import sparse
    L = len(y)
    D = sparse.csc_matrix(np.diff(np.eye(L), 2))
    w = np.ones(L)
    for i in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = sparse.linalg.spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z


def baseline_als(x, y, lam=None, p=None, niter=10, return_baseline=False,
                 offset_correction=True,
                 outlier_purging=True):
    """Baseline Correction with Asymmetric Least Squares Smoothing."""

    if lam is None:
        lam = 1e9
    if p is None:
        p = 0.001

    if outlier_purging:
        y = purge_outliers(y)

    z = _als(y, lam, p, niter=niter)

    # _, z2 = baseline_rough(x, y - z, return_baseline=True)
    # z += z2

    ysub = y - z
    offset = 0
    if offset_correction:
        std = ref_std(ysub, np.max([len(y) // 20, 20]))

        good = np.abs(ysub) < 10 * std

        offset = offset_fit(x[good], ysub[good], 0)

    if return_baseline:
        return ysub - offset, z + offset
    else:
        return ysub - offset


def fit_baseline_plus_bell(x, y, ye=None, kind='gauss'):
    """Fit a function composed of a linear baseline plus a bell function.

    kind:     'gauss' or 'lorentz'
    """
    assert kind in ['gauss', 'lorentz'], \
        'kind has to be one of: gauss, lorentz'
    from astropy.modeling import models, fitting

    base = models.Linear1D(slope=0, intercept=np.min(y), name='Baseline')

    xrange = np.max(x) - np.min(x)
    yrange = np.max(y) - np.min(y)

    if kind == 'gauss':
        bell = models.Gaussian1D(mean=np.mean(x), stddev=xrange / 20,
                                 amplitude=yrange, name='Bell')
        bell.amplitude.bounds = (0, None)
        bell.mean.bounds = (None, None)
        bell.stddev.bounds = (0, None)
        # max_name = 'mean'
    elif kind == 'lorentz':
        bell = models.Lorentz1D(x_0=np.mean(x), fwhm=xrange / 20,
                                amplitude=yrange, name='Bell')
        bell.amplitude.bounds = (0, None)
        bell.x_0.bounds = (None, None)
        bell.fwhm.bounds = (0, None)
        # max_name = 'x_0'

    mod_init = base + bell

    fit = fitting.LevMarLSQFitter()

    mod_out = fit(mod_init, x, y)

    return mod_out, fit.fit_info


def minimize_align(xs, ys, params):
    """Calculate the total variance of a series of scans.

    This functions subtracts a linear function from each of the scans
    (excluding the first one) and calculates the total variance.
    """
    params = np.array(params).flatten()
    qs = params[:len(xs) - 1]
    ms = params[len(xs) - 1:]

    x = xs[0].copy()
    y = ys[0].copy()
    for i in range(1, len(xs)):
        x = np.append(x, xs[i])
        scaled_y = ys[i] - (xs[i] * ms[i - 1] + qs[i - 1])
        y = np.append(y, scaled_y)

    order = np.argsort(x)

    x = x[order]
    y = y[order]

    x_range = [np.min(x), np.max(x)]

    xints = np.linspace(x_range[0], x_range[1], len(x) / 20)

    values = np.array([np.var(y[(x >= xints[k]) & (x < xints[k+1])])
                      for k in range(len(xints[:-1]))])
    good = values == values
    value = np.mean(values[good])

    return value


def objective_function(params, args):
    """Put the parameters in the right order to use with scipy's minimize."""
    return minimize_align(args[0], args[1], params)


def align(xs, ys):
    """Given the first scan, it aligns all the others to that."""
    from scipy.optimize import minimize

    qs = np.zeros(len(xs) - 1)
    ms = np.zeros(len(xs) - 1)

    result = minimize(objective_function, [qs, ms], args=[xs, ys],
                      options={'disp': True})

    qs = result.x[:len(xs) - 1]
    ms = np.zeros(len(xs) - 1)

    result = minimize(objective_function, [qs, ms], args=[xs, ys],
                      options={'disp': True})

    qs = result.x[:len(xs) - 1]
    ms = result.x[len(xs) - 1:]

    return qs, ms
