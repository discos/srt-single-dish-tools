"""Useful fitting functions."""
from __future__ import (absolute_import, division,
                        print_function)
from scipy.optimize import curve_fit
from scipy.signal import medfilt
import numpy as np
import traceback
import warnings
import collections
import copy
from .utils import mad, HAS_MPL


__all__ = ["contiguous_regions", "ref_std", "ref_mad", "linear_fun",
           "linear_fit", "offset", "offset_fit", "baseline_rough",
           "purge_outliers", "baseline_als", "fit_baseline_plus_bell",
           "total_variance", "align"]


def contiguous_regions(condition):
    """Find contiguous True regions of the boolean array "condition".

    Return a 2D array where the first column is the start index of the region
    and the second column is the end index.

    Parameters
    ----------
    condition : boolean array

    Returns
    -------
    idx : [[i0_0, i0_1], [i1_0, i1_1], ...]
        A list of integer couples, with the start and end of each True blocks
        in the original array

    Notes
    -----
    From http://stackoverflow.com/questions/4494404/
        find-large-number-of-consecutive-values-fulfilling-
        condition-in-a-numpy-array
    """  # NOQA
    # Find the indicies of changes in "condition"
    diff = np.logical_xor(condition[1:], condition[:-1])
    idx, = diff.nonzero()
    # We need to start things after the change in "condition". Therefore,
    # we'll shift the index by 1 to the right.
    idx += 1
    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = np.r_[0, idx]
    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = np.r_[idx, condition.size]
    # Reshape the result into two columns
    idx.shape = (-1, 2)
    return idx


def _rolling_window(a, window):
    """A smart rolling window.

    Found at http://www.rigtorp.se/2011/01/01/rolling-statistics-numpy.html
    """
    try:
        shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
        strides = a.strides + (a.strides[-1],)
        return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)
    except Exception:
        warnings.warn(traceback.format_exc())
        raise


def ref_std(array, window):
    """Minimum standard deviation along an array.

    If a data series is noisy, it is difficult to determine the underlying
    standard deviation of the original series. Here, the standard deviation is
    calculated in a rolling window, and the minimum is saved, because it will
    likely be the interval with less noise.

    Parameters
    ----------
    array : ``numpy.array`` object or list
        Input data
    window : int or float
        Number of bins of the window

    Returns
    -------
    ref_std : float
        The reference Standard Deviation
    """

    return np.std(np.diff(array)) / np.sqrt(2)


def ref_mad(array, window):
    """Ref. Median Absolute Deviation of an array, rolling median-subtracted.

    If a data series is noisy, it is difficult to determine the underlying
    statistics of the original series. Here, the MAD is calculated in a rolling
    window, and the minimum is saved, because it will likely be the interval
    with less noise.

    Parameters
    ----------
    array : ``numpy.array`` object or list
        Input data
    window : int or float
        Number of bins of the window

    Returns
    -------
    ref_std : float
        The reference MAD
    """
    return mad(np.diff(array)) / np.sqrt(2)


def linear_fun(x, q, m):
    """A linear function.

    Parameters
    ----------
    x : float or array
        The independent variable
    m : float
        The slope
    q : float
        The intercept

    Returns
    -------
    y : float or array
        The dependent variable
    """
    return m * np.asarray(x, dtype=np.float) + q


def linear_fit(x, y, start_pars, return_err=False):
    """A linear fit with any set of data.

    Parameters
    ----------
    x : array-like
    y : array-like
    start_pars : [q0, m0], floats
        Intercept and slope of linear function

    Returns
    -------
    par : [q, m], floats
        Fitted intercept and slope of the linear function
    """
    par, _ = curve_fit(linear_fun, x, y, start_pars,
                       maxfev=6000)
    if return_err:
        warnings.warn("return_err not implemented yet in linear_fit")
        return par, None
    else:
        return par


def offset(x, off):
    """An offset."""
    return off


def offset_fit(x, y, offset_start=0, return_err=False):
    """Fit a constant offset to the data.

    Parameters
    ----------
    x : array-like
    y : array-like
    offset_start : float
        Constant offset, initial value

    Returns
    -------
    offset : float
        Fitted offset
    """
    par, _ = curve_fit(offset, x, y, [offset_start],
                       maxfev=6000)
    if return_err:
        warnings.warn("return_err not implemented yet in offset_fit")
        return par[0], None
    else:
        return par[0]


def baseline_rough(x, y, start_pars=None, return_baseline=False, mask=None):
    """Rough function to subtract the baseline.

    Parameters
    ----------
    x : array-like
        the sample time/number/position
    y : array-like
        the data series corresponding to x
    start_pars : [q0, m0], floats
        Intercept and slope of linear function

    Other Parameters
    ----------------
    return_baseline : bool
        return the baseline?
    mask : array of bools
        Mask indicating the good x and y data. True for good, False for bad

    Returns
    -------
    y_subtracted : array-like, same size as y
        The initial time series, subtracted from the trend
    baseline : array-like, same size as y
        Fitted baseline
    """
    N = len(y)
    if start_pars is None:
        if N > 40:
            m0 = (np.median(y[-20:]) - np.median(y[:20])) / \
                       (np.mean(x[-20:]) - np.mean(x[:20]))
        else:
            m0 = (y[-1] - y[0])/(x[-1] - x[0])

        q0 = min(y)
        start_pars = [q0, m0]

    lc = y.copy()
    time = x.copy()

    if mask is None:
        mask = np.ones(len(time), dtype=bool)
    total_trend = 0

    if N < 20:
        par = linear_fit(time, lc, start_pars)

        lc = lc - linear_fun(time, *par)
        total_trend = total_trend + linear_fun(time, *par)
    else:
        local_std = ref_std(lc, np.max([N // 20, 20]))

        for percentage in [0.8, 0.15]:
            time_to_fit = time[mask][1:-1]
            lc_to_fit = lc[mask][1:-1]
            if len(time_to_fit) < len(start_pars):
                break

            sorted_els = np.argsort(lc_to_fit)
            # Select the lowest half elements
            good = sorted_els[: int(N * percentage)]

            if np.std(lc_to_fit[good]) < 2 * local_std:
                good = np.ones(len(lc_to_fit), dtype=bool)

            time_filt = time_to_fit[good]
            lc_filt = lc_to_fit[good]
            if len(time_filt) < len(start_pars):
                break
            back_in_order = np.argsort(time_filt)
            lc_filt = lc_filt[back_in_order]
            time_filt = time_filt[back_in_order]

            par = linear_fit(time_filt, lc_filt, start_pars)

            lc = lc - linear_fun(time, *par)
            total_trend = total_trend + linear_fun(time, *par)

    if return_baseline:
        return lc, total_trend
    else:
        return lc


def outlier_from_median_filt(y, window_size, down=True, up=True):
    diffs = y - medfilt(y, window_size)
    min_diff = mad(diffs)

    outliers = np.zeros(len(y), dtype=bool)
    if down:
        outliers = np.logical_or(outliers, -diffs > 10 * min_diff)
    if up:
        outliers = np.logical_or(outliers, diffs > 10 * min_diff)

    return outliers


def purge_outliers(y, window_size=5, up=True, down=True, mask=None,
                   plot=False):
    """Remove obvious outliers.

    Attention: This is known to throw false positives on bona fide, very strong
    Gaussian peaks
    """
    # Needs to be odd
    window_size = window_size // 2 * 2 + 1

    if mask is None:
        mask = np.ones(len(y), dtype=bool)

    bad_mask = np.logical_not(mask)
    if not (up or down):
        return y

    ysave = y
    y = y.copy()

    win1 = outlier_from_median_filt(y, window_size)
    win2 = outlier_from_median_filt(y, window_size * 2 + 1)
    local_outliers = win1 & win2

    Noutliers = len(local_outliers[local_outliers])
    if Noutliers > 0:
        warnings.warn("Found {} outliers".format(
            Noutliers),
                      UserWarning)

    outliers = np.logical_or(local_outliers, bad_mask)
    if not np.any(outliers):
        return y

    bad = contiguous_regions(outliers)
    for b in bad:
        if b[0] == 0:
            y[b[0]] = y[b[1]]
        elif b[1] >= len(y):
            y[b[0]:] = y[b[0] - 1]
        else:
            previous = y[b[0] - 1]
            next_bin = y[b[1]]
            dx = b[1] - b[0]
            y[b[0]:b[1]] = \
                (next_bin - previous)/(dx + 1) * \
                np.arange(1, b[1] - b[0] + 1) + previous

    if plot and HAS_MPL:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        plt.plot(ysave)
        plt.plot(y, zorder=3)
        plt.plot(medfilt(ysave, window_size), zorder=6, lw=1)
        plt.savefig("Bubu_" + str(np.random.randint(0, 10000000)) + '.png')
        plt.close(fig)

    return y


def _als(y, lam, p, niter=30):
    """Baseline Correction with Asymmetric Least Squares Smoothing.

    Modifications to the routine from Eilers & Boelens 2005
    https://www.researchgate.net/publication/
        228961729_Technical_Report_Baseline_Correction_with_
        Asymmetric_Least_Squares_Smoothing
    The Python translation is partly from
    http://stackoverflow.com/questions/29156532/
        python-baseline-correction-library

    Parameters
    ----------
    y : array-like
        the data series corresponding to x
    lam : float
        the lambda parameter of the ALS method. This control how much the
        baseline can adapt to local changes. A higher value corresponds to a
        stiffer baseline
    p : float
        the asymmetry parameter of the ALS method. This controls the overall
        slope tollerated for the baseline. A higher value correspond to a
        higher possible slope

    Other parameters
    ----------------
    niter : int
        The number of iterations to perform

    Returns
    -------
    z : array-like, same size as y
        Fitted baseline.
    """
    from scipy import sparse
    L = len(y)
    D = sparse.csc_matrix(np.diff(np.eye(L), 2))
    w = np.ones(L)
    for _ in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = sparse.linalg.spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z


def baseline_als(x, y, lam=None, p=None, niter=40, return_baseline=False,
                 offset_correction=True, mask=None,
                 outlier_purging=True):
    """Baseline Correction with Asymmetric Least Squares Smoothing.

    Parameters
    ----------
    x : array-like
        the sample time/number/position
    y : array-like
        the data series corresponding to x
    lam : float
        the lambda parameter of the ALS method. This control how much the
        baseline can adapt to local changes. A higher value corresponds to a
        stiffer baseline
    p : float
        the asymmetry parameter of the ALS method. This controls the overall
        slope tollerated for the baseline. A higher value correspond to a
        higher possible slope

    Other Parameters
    ----------------
    niter : int
        The number of iterations to perform
    return_baseline : bool
        return the baseline?
    offset_correction : bool
        also correct for an offset to align with the running mean of the scan
    outlier_purging : bool
        Purge outliers before the fit?
    mask : array of bools
        Mask indicating the good x and y data. True for good, False for bad

    Returns
    -------
    y_subtracted : array-like, same size as y
        The initial time series, subtracted from the trend
    baseline : array-like, same size as y
        Fitted baseline. Only returned if return_baseline is True
    """

    if not isinstance(outlier_purging, collections.Iterable):
        outlier_purging = (outlier_purging, outlier_purging)
    if lam is None:
        lam = 1e11
    if p is None:
        p = 0.001

    N = len(y)
    if N > 40:
        med_start = np.median(y[:20])
        med_stop = np.median(y[-20:])
        approx_m = (med_stop - med_start) / (N - 20)
    else:
        approx_m = (y[-1] - y[0]) / (N - 1)

    approx_q = y[0]

    approx_baseline = approx_m * np.arange(N) + approx_q
    y = y - approx_baseline

    y_mod = purge_outliers(y, up=outlier_purging[0],
                           down=outlier_purging[1],
                           mask=mask)

    z = _als(y_mod, lam, p, niter=niter)

    offset = 0
    ysub = y_mod - z
    if offset_correction:
        std = ref_std(ysub, np.max([len(y) // 20, 20]))

        good = np.abs(ysub) < 10 * std

        if len(ysub[good]) < 20:
            good = np.ones(len(ysub), dtype=bool)

        offset = np.median(ysub[good])
        if np.isnan(offset):
            offset = 0

    if return_baseline:
        return y - z - offset, z + offset + approx_baseline
    else:
        return y - z - offset


def detrend_spectroscopic_data(x, spectrum, kind='als', outlier_purging=True):
    """Take the baseline off the spectroscopic data.

    Examples
    --------
    >>> spectrum = np.vstack([np.arange(0 + i, 2 + i, 1/3)
    ...                       for i in np.arange(0., 4, 1/16)])
    >>> x = np.arange(spectrum.shape[0])
    >>> detr, _ = detrend_spectroscopic_data(x, spectrum, kind='rough')
    >>> np.allclose(detr, 0, atol=1e-3)
    True
    """
    y = np.sum(spectrum, axis=1)
    if kind == 'als':
        y_sub, baseline = baseline_als(x, y, return_baseline=True,
                                       outlier_purging=outlier_purging)
    elif kind == 'rough':
        y_sub, baseline = baseline_rough(x, y, return_baseline=True)
    else:
        warnings.warn('Baseline kind unknown')
        return spectrum, np.ones_like(spectrum)

    if len(spectrum.shape) == 1:
        return y_sub, baseline

    shape = spectrum.shape
    tiled_baseline = np.tile(baseline, (shape[1], 1)).transpose()
    tiled_norm = np.tile(y, (shape[1], 1)).transpose()

    tiled_baseline = tiled_baseline / tiled_norm * spectrum

    return spectrum - tiled_baseline, tiled_baseline


def fit_baseline_plus_bell(x, y, ye=None, kind='gauss'):
    """Fit a function composed of a linear baseline plus a bell function.

    Parameters
    ----------
    x : array-like
        the sample time/number/position
    y : array-like
        the data series corresponding to x

    Other parameters
    ----------------
    ye : array-like
        the errors on the data series
    kind: str
        Can be 'gauss' or 'lorentz'

    Returns
    -------
    mod_out : ``Astropy.modeling.model`` object
        The fitted model
    fit_info : dict
        Fit info from the Astropy fitting routine.
    """
    if kind not in ['gauss', 'lorentz']:
        raise ValueError('kind has to be one of: gauss, lorentz')
    from astropy.modeling import models, fitting

    approx_m = (np.median(y[-20:]) - np.median(y[:20])) / \
               (np.mean(x[-20:]) - np.mean(x[:20]))

    base = models.Linear1D(slope=approx_m, intercept=np.min(y),
                           name='Baseline')

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


def total_variance(xs, ys, params):
    """Calculate the total variance of a series of scans.

    This functions subtracts a linear function from each of the scans
    (excluding the first one) and calculates the total variance.

    Parameters
    ----------
    xs : list of array-like [array1, array2, ...]
        list of arrays containing the x values of each scan
    ys : list of array-like [array1, array2, ...]
        list of arrays containing the y values of each scan
    params : list of array-like [[q0, m0], [q1, m1], ...]
        list of arrays containing the parameters [m, q] for each scan.

    Returns
    -------
    total_variance : float
        The total variance of the baseline-subtracted scans.
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


def _objective_function(params, args):
    """Put the parameters in the right order to use with scipy's minimize."""
    return total_variance(args[0], args[1], params)


def align(xs, ys):
    """Given the first scan, it aligns all the others to that.

    Parameters
    ----------
    xs : list of array-like [array1, array2, ...]
        list of arrays containing the x values of each scan
    ys : list of array-like [array1, array2, ...]
        list of arrays containing the y values of each scan

    Returns
    -------
    qs : array-like
        The list of intercepts maximising the alignment, one for each scan
    ms : array-like
        The list of slopes maximising the alignment, one  for each scan
    """
    from scipy.optimize import minimize

    qs = np.zeros(len(xs) - 1)
    ms = np.zeros(len(xs) - 1)

    result = minimize(_objective_function, [qs, ms], args=[xs, ys],
                      options={'disp': True})

    qs = result.x[:len(xs) - 1]
    ms = np.zeros(len(xs) - 1)

    result = minimize(_objective_function, [qs, ms], args=[xs, ys],
                      options={'disp': True})

    qs = result.x[:len(xs) - 1]
    ms = result.x[len(xs) - 1:]

    return qs, ms
