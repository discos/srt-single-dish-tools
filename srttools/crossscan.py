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

import warnings

import numpy as np
from scipy.interpolate import BSpline

from astropy.table import Table

from .io import get_channel_feed

__all__ = [
    "arm_offsets",
    "binned_profile",
    "combine_amplitudes",
    "consistency",
    "fit_subscan",
    "measure_crosses",
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


ARM_NAMES = {True: "RA", False: "Dec"}


def arm_offsets(table, feed, ra0, dec0):
    """Offset of each sample from the source along its arm, in degrees.

    Parameters
    ----------
    table : `astropy.table.Table`
        With ``ra`` and ``dec`` columns (radians, one column per feed) and
        ``direction`` (True for the RA arm)
    feed : int
        Feed index
    ra0, dec0 : float
        Source position, radians

    Returns
    -------
    offset : array
        RA offset times cos(dec0) on the RA arm, Dec offset on the Dec arm
    arm : array of bool
        True on the RA arm
    """
    arm = np.asarray(table["direction"], dtype=bool)
    ra = np.asarray(table["ra"])[:, feed]
    dec = np.asarray(table["dec"])[:, feed]
    dra = (ra - ra0 + np.pi) % (2 * np.pi) - np.pi
    offset = np.where(arm, dra * np.cos(dec0), dec - dec0)
    return np.degrees(offset), arm


def _channel_average(values):
    """Average over the channel axis (axis 1), ignoring NaNs."""
    n = np.sum(np.isfinite(values), axis=1)
    total = np.nansum(values, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n > 0, total / n, np.nan)


def _profile_stats(binned):
    """Mean over subscans (axis 0) and its error, ignoring NaNs."""
    n = np.sum(np.isfinite(binned), axis=0)
    profile = np.nansum(binned, axis=0) / np.maximum(n, 1)
    err = np.sqrt(
        np.nansum((binned - profile) ** 2, axis=0) / np.maximum(n - 1, 1) / np.maximum(n, 1)
    )
    profile[n == 0] = np.nan
    err[n < 2] = np.nan
    return profile, err


def measure_crosses(
    table,
    ra0,
    dec0,
    beam_fwhm,
    channels,
    suffix="",
    rfi_nsigma=5.0,
    statistic="mean",
    nsigma=3.0,
    cal_rel_err=0.0,
    bin_width=None,
    **kwargs,
):
    """Flux, error and upper limit of a point source from cross scans.

    Parameters
    ----------
    table : `astropy.table.Table`
        Cross-scan data, as in a ScanSet: ``ra``, ``dec`` (radians),
        ``direction`` (True for the RA arm), ``Scan_id``, the data columns
        and, optionally, the ``<channel>-filt`` columns (False for samples
        to ignore)
    ra0, dec0 : float
        Position where the flux is measured (radians)
    beam_fwhm : float
        Beam FWHM, degrees
    channels : list of str
        Channels to use

    Other Parameters
    ----------------
    suffix : str
        Data are read from the columns ``<channel><suffix>``, e.g. ``"_cal"``
    rfi_nsigma : float
        Ignore subscans with anomalous noise (residual RFI), see
        `srttools.imager.find_outlier_subscans`. ``None`` disables it
    statistic : str
        How subscan amplitudes are combined, see `combine_amplitudes`
    nsigma : float
        Detection threshold and width of the upper limit, in standard
        deviations
    cal_rel_err : float
        Relative calibration error, added in quadrature in upper limits
    bin_width : float
        Width of the profile bins, degrees. Default: 1/4 beam FWHM
    **kwargs
        Passed to `fit_subscan` (``window``, ``knot_spacing`` in degrees,
        ``joint``)

    Returns
    -------
    result : dict
        ``summary``: table with one row per channel and arm, per channel,
        per arm and in total (``label`` "TOTAL"); columns ``flux``, ``err``,
        ``nsub`` (subscans used), ``ul`` (upper limit), ``detected``.
        ``subscans``: table of the amplitude of each subscan and channel.
        ``consistency``: difference between the arms and between the first
        two channels, in standard deviations.
        ``profiles``: (bin centers, profile, error) for each channel and
        arm (key ``"<channel> <arm>"``), for each arm and for ``"TOTAL"``.
        ``flagged``: subscans rejected for anomalous noise, per channel
    """
    from .imager import find_outlier_subscans

    beam_sigma = beam_fwhm / (2 * np.sqrt(2 * np.log(2)))
    if bin_width is None:
        bin_width = beam_fwhm / 4
    ids = np.unique(table["Scan_id"])
    scan_id = np.asarray(table["Scan_id"])
    nchan = len(channels)
    amps = np.full((ids.size, nchan), np.nan)
    sub_arm = np.zeros(ids.size, dtype=bool)

    offsets = {}
    max_offset = 0.0
    for ch in channels:
        offset, arm = arm_offsets(table, get_channel_feed(ch), ra0, dec0)
        offsets[ch] = offset, arm
        max_offset = max(max_offset, np.max(np.abs(offset)))
    nedge = int(np.ceil(max_offset / bin_width - 0.5))
    edges = (np.arange(-nedge, nedge + 2) - 0.5) * bin_width
    binned = np.full((ids.size, nchan, edges.size - 1), np.nan)

    flagged = {}
    for ich, ch in enumerate(channels):
        offset, arm = offsets[ch]
        values = np.asarray(table[ch + suffix], dtype=float)
        good = np.isfinite(values)
        if ch + "-filt" in table.colnames:
            good &= np.asarray(table[ch + "-filt"], dtype=bool)
        flagged[ch] = []
        if rfi_nsigma is not None:
            flagged[ch] = find_outlier_subscans(
                values[good],
                offset[good] / bin_width,
                scan_id[good],
                arm[good],
                nsigma=rfi_nsigma,
            )
            good &= ~np.isin(scan_id, flagged[ch])
            if len(flagged[ch]) > 0:
                warnings.warn(
                    f"RFI rejection: {ch}: ignored {len(flagged[ch])}/{ids.size} subscans "
                    f"with anomalous noise (Scan_id {flagged[ch]})."
                )
        for i, sid in enumerate(ids):
            mask = (scan_id == sid) & good
            if not np.any(mask):
                continue
            sub_arm[i] = arm[mask][0]
            amp, resid = fit_subscan(offset[mask], values[mask], beam_sigma, **kwargs)
            amps[i, ich] = amp
            ok = np.isfinite(resid)
            if not np.any(ok):
                continue
            num = np.histogram(offset[mask][ok], bins=edges, weights=resid[ok])[0]
            den = np.histogram(offset[mask][ok], bins=edges)[0]
            binned[i, ich, den > 0] = num[den > 0] / den[den > 0]

    groups = []
    for ich, ch in enumerate(channels):
        for a in [True, False]:
            groups.append((f"{ch} {ARM_NAMES[a]}", ch, ARM_NAMES[a], [ich], sub_arm == a))
        groups.append((ch, ch, "", [ich], np.ones(ids.size, dtype=bool)))
    for a in [True, False]:
        groups.append((ARM_NAMES[a], "", ARM_NAMES[a], list(range(nchan)), sub_arm == a))
    groups.append(("TOTAL", "", "", list(range(nchan)), np.ones(ids.size, dtype=bool)))

    rows = []
    estimates = {}
    profiles = {}
    for label, ch, armname, ichans, subs in groups:
        per_sub = _channel_average(amps[subs][:, ichans])
        flux, err = combine_amplitudes(per_sub, statistic=statistic)
        estimates[label] = flux, err
        rows.append(
            dict(
                label=label,
                channel=ch,
                arm=armname,
                flux=flux,
                err=err,
                nsub=int(np.sum(np.isfinite(per_sub))),
                ul=upper_limit(flux, err, nsigma=nsigma, cal_rel_err=cal_rel_err),
                detected=bool(flux > nsigma * err),
            )
        )
        if ch == "" or armname != "":
            prof = _channel_average(binned[subs][:, ichans])
            profiles[label] = ((edges[:-1] + edges[1:]) / 2, *_profile_stats(prof))

    consist = {"RA-Dec": consistency(estimates["RA"], estimates["Dec"])}
    if nchan >= 2:
        diff = amps[:, 0] - amps[:, 1]
        dflux, derr = combine_amplitudes(diff, statistic=statistic)
        consist[f"{channels[0]}-{channels[1]}"] = dflux / derr

    subscans = Table()
    subscans["Scan_id"] = ids
    subscans["arm"] = [ARM_NAMES[a] for a in sub_arm]
    for ich, ch in enumerate(channels):
        subscans[ch] = amps[:, ich]

    return dict(
        summary=Table(rows=rows),
        subscans=subscans,
        consistency=consist,
        profiles=profiles,
        flagged=flagged,
    )
