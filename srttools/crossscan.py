"""Flux and upper limits of a point source from cross scans.

Each subscan is treated separately:

1. a smooth baseline (cubic least-squares spline with fixed, widely spaced
   knots) is removed;
2. the amplitude of a Gaussian of known width, centred at the **nominal**
   source position, is measured (a *forced fit*: the amplitude may be
   negative, nothing else is free).

The flux is then the average of the subscan amplitudes, and its error comes
from the scatter of those amplitudes between subscans. This error is correct
whatever the noise (white, 1/f, baseline drifts), as long as the subscans are
independent, which is not true for the samples inside a subscan.

Two ways of removing the baseline are available:

- ``joint=False``: the spline is fitted only outside a window of half-width
  ``window`` around the source, and the amplitude is measured afterwards on
  the residuals;
- ``joint=True``: spline and Gaussian are fitted together on all samples
  (one linear least-squares fit). The knots must be spaced widely enough
  that the spline cannot mimic the source.

All coordinates are in the same unit (e.g. pixels) and relative to the
nominal source position, unless ``center`` is given.
"""

import numpy as np
from scipy.interpolate import BSpline

__all__ = [
    "binned_profile",
    "combine_amplitudes",
    "consistency",
    "fit_subscan",
    "subscan_amplitudes",
    "upper_limit",
]

SPLINE_ORDER = 3


def _knots(lo, hi, knot_spacing, center=0.0, window=0.0):
    """Interior knots, equally spaced, none inside the window around the source."""
    n_interior = max(int(round((hi - lo) / knot_spacing)) - 1, 0)
    interior = np.linspace(lo, hi, n_interior + 2)[1:-1]
    interior = interior[np.abs(interior - center) >= window]
    k = SPLINE_ORDER
    return np.concatenate([[lo] * (k + 1), interior, [hi] * (k + 1)])


def _spline_basis(x, knots):
    return BSpline.design_matrix(x, knots, SPLINE_ORDER).toarray()


def fit_subscan(
    x, y, beam_sigma, center=0.0, window=None, knot_spacing=None, joint=False, min_points=20
):
    """Forced-fit amplitude of a Gaussian over a spline baseline, in one subscan.

    Parameters
    ----------
    x : array
        Coordinate along the scan
    y : array
        Data. Non-finite values are ignored
    beam_sigma : float
        Standard deviation of the Gaussian beam, same unit as ``x``

    Other Parameters
    ----------------
    center : float
        Nominal source position
    window : float
        Half-width of the region around the source excluded from the
        baseline fit (``joint=False``) and free of knots. Default: 1.5 beam
        FWHM
    knot_spacing : float
        Approximate distance between the spline knots. Default: 3 beam FWHM
    joint : bool
        Fit baseline and Gaussian together on all samples
    min_points : int
        Minimum number of samples used in the fit beyond the number of free
        parameters. Below it, NaN is returned

    Returns
    -------
    amplitude : float
        Peak of the Gaussian (NaN if the fit was impossible)
    resid : array
        ``y`` minus the baseline, so that the source is still there (NaN if
        the fit was impossible)
    """
    fwhm = beam_sigma * 2 * np.sqrt(2 * np.log(2))
    if window is None:
        window = 1.5 * fwhm
    if knot_spacing is None:
        knot_spacing = 3 * fwhm
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    resid = np.full(y.shape, np.nan)
    good = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(good) < min_points:
        return np.nan, resid
    xg, yg = x[good], y[good]
    knots = _knots(xg.min(), xg.max(), knot_spacing, center=center, window=window)
    basis = _spline_basis(xg, knots)
    template = np.exp(-((xg - center) ** 2) / (2 * beam_sigma**2))
    if joint:
        design = np.column_stack([basis, template])
        if xg.size < design.shape[1] + min_points:
            return np.nan, resid
        coef = np.linalg.lstsq(design, yg, rcond=None)[0]
        resid[good] = yg - basis @ coef[:-1]
        return coef[-1], resid

    outside = np.abs(xg - center) > window
    if np.count_nonzero(outside) < basis.shape[1] + min_points:
        return np.nan, resid
    coef = np.linalg.lstsq(basis[outside], yg[outside], rcond=None)[0]
    resid[good] = yg - basis @ coef
    amplitude = np.sum(template * resid[good]) / np.sum(template**2)
    return amplitude, resid


def subscan_amplitudes(x, y, scan_id, beam_sigma, **kwargs):
    """Run `fit_subscan` on each subscan.

    Parameters
    ----------
    x, y : array
        Coordinate and data of all samples
    scan_id : array of int
        Subscan each sample belongs to
    beam_sigma : float
        Standard deviation of the Gaussian beam
    **kwargs
        Passed to `fit_subscan`

    Returns
    -------
    ids : array
        Unique values of ``scan_id``
    amplitudes : array
        Amplitude in each subscan (NaN where the fit was impossible)
    """
    ids = np.unique(scan_id)
    amps = np.array(
        [fit_subscan(x[scan_id == i], y[scan_id == i], beam_sigma, **kwargs)[0] for i in ids]
    )
    return ids, amps


def combine_amplitudes(amplitudes, statistic="mean"):
    """Combine subscan amplitudes, with the error from their scatter.

    Parameters
    ----------
    amplitudes : array
        One value per subscan; NaNs are ignored

    Other Parameters
    ----------------
    statistic : str
        ``"mean"``: mean, error = standard deviation / sqrt(N).
        ``"median"``: median, error = sqrt(pi / 2) times the robust standard
        deviation (1.4826 times the median absolute deviation) / sqrt(N)

    Returns
    -------
    flux, err : float
    """
    a = np.asarray(amplitudes, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return np.nan, np.nan
    if a.size < 2:
        return a[0], np.nan
    if statistic == "mean":
        return np.mean(a), np.std(a, ddof=1) / np.sqrt(a.size)
    if statistic == "median":
        med = np.median(a)
        sdev = 1.4826 * np.median(np.abs(a - med))
        return med, np.sqrt(np.pi / 2) * sdev / np.sqrt(a.size)
    raise ValueError(f"Unknown statistic: {statistic}")


def consistency(first, second):
    """Difference between two (value, error) estimates, in units of its error."""
    (f1, e1), (f2, e2) = first, second
    return (f1 - f2) / np.sqrt(e1**2 + e2**2)


def upper_limit(flux, err, nsigma=3.0, cal_rel_err=0.0):
    """Upper limit on the flux: max(flux, 0) + nsigma * total error.

    The total error is the statistical error and the calibration error
    (``cal_rel_err`` times max(flux, 0)) added in quadrature.
    """
    f = max(flux, 0.0)
    return f + nsigma * np.sqrt(err**2 + (cal_rel_err * f) ** 2)


def binned_profile(x, y, scan_id, beam_sigma, bin_width=1.0, center=0.0, edges=None, **kwargs):
    """Average baseline-subtracted profile, with errors from the subscan scatter.

    Each subscan is baseline-subtracted with `fit_subscan` and binned; the
    profile is the mean of the binned subscans and the error is their
    standard deviation / sqrt(N), bin by bin.

    Parameters
    ----------
    x, y : array
        Coordinate and data of all samples
    scan_id : array of int
        Subscan each sample belongs to
    beam_sigma : float
        Standard deviation of the Gaussian beam

    Other Parameters
    ----------------
    bin_width : float
        Width of the bins. One bin is centred on ``center``
    edges : array
        Bin edges, overriding ``bin_width`` (e.g. to put several profiles on
        the same grid)
    **kwargs
        Passed to `fit_subscan`

    Returns
    -------
    centers, profile, err : array
    """
    x = np.asarray(x, dtype=float)
    if edges is None:
        good = np.isfinite(x)
        kmin = np.round((x[good].min() - center) / bin_width)
        kmax = np.round((x[good].max() - center) / bin_width)
        edges = center + (np.arange(kmin, kmax + 2) - 0.5) * bin_width
    ids = np.unique(scan_id)
    binned = np.full((ids.size, edges.size - 1), np.nan)
    for i, sid in enumerate(ids):
        mask = scan_id == sid
        _, resid = fit_subscan(x[mask], y[mask], beam_sigma, center=center, **kwargs)
        ok = np.isfinite(resid)
        if not np.any(ok):
            continue
        num = np.histogram(x[mask][ok], bins=edges, weights=resid[ok])[0]
        den = np.histogram(x[mask][ok], bins=edges)[0]
        binned[i, den > 0] = num[den > 0] / den[den > 0]
    n = np.sum(np.isfinite(binned), axis=0)
    profile = np.nansum(binned, axis=0) / np.maximum(n, 1)
    err = np.sqrt(
        np.nansum((binned - profile) ** 2, axis=0) / np.maximum(n - 1, 1) / np.maximum(n, 1)
    )
    profile[n == 0] = np.nan
    err[n < 2] = np.nan
    centers = (edges[:-1] + edges[1:]) / 2
    return centers, profile, err
