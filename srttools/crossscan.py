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
- ``joint=True`` (default): spline and Gaussian are fitted together on all
  samples (one linear least-squares fit). The knots must be spaced widely
  enough that the spline cannot mimic the source.

The defaults (joint fit, knots every 3 beam FWHM, none within 1.5 FWHM of
the source) were chosen with a Monte Carlo simulation (white and 1/f noise,
baseline drifts, RFI, pointing offsets): all the variants tested give correct
errors, and this one gives the smallest (about 30% smaller than fitting the
baseline outside a window of 1.5 FWHM). See ``docs/technical_details.rst``.

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
    "beam_fwhm_from_frequency",
    "binned_profile",
    "calibrate_crosses",
    "combine_amplitudes",
    "consistency",
    "fit_subscan",
    "main_cross",
    "measure_crosses",
    "plot_cross_results",
    "subscan_amplitudes",
    "upper_limit",
    "write_cross_results",
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
    x, y, beam_sigma, center=0.0, window=None, knot_spacing=None, joint=True, min_points=20
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
        Fit baseline and Gaussian together on all samples (default). If
        False, fit the baseline outside the window first
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
    standard deviation / sqrt(N), bin by bin. Bins covered by less than half
    of the subscans are NaN.

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
    centers = (edges[:-1] + edges[1:]) / 2
    return (centers, *_profile_stats(binned))


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


def _profile_stats(binned, min_frac=0.5):
    """Mean over subscans (axis 0) and its error, ignoring NaNs.

    Bins covered by less than ``min_frac`` times the subscans of the
    best-covered bin (typically at the ends of the scans) are set to NaN.
    """
    n = np.sum(np.isfinite(binned), axis=0)
    binned = np.where(n >= min_frac * n.max(), binned, np.nan)
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
        ``consistency``: difference between the arms, between the first
        two channels and (with the joint fit) between the joint and window
        fits (``"joint-window"``, total over all subscans), in standard
        deviations. The last is negative for a source that does not look
        like the beam at the nominal position (pointing offset, extended
        source), whose flux the joint fit underestimates.
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
    joint = kwargs.get("joint", True)
    amps_window = np.full((ids.size, nchan), np.nan)
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
                    f"with anomalous noise (Scan_id {flagged[ch]}). Use --rfi-nsigma 0 to disable."
                )
        for i, sid in enumerate(ids):
            mask = (scan_id == sid) & good
            if not np.any(mask):
                continue
            sub_arm[i] = arm[mask][0]
            amp, resid = fit_subscan(offset[mask], values[mask], beam_sigma, **kwargs)
            amps[i, ich] = amp
            if joint:
                window_kwargs = {**kwargs, "joint": False}
                amps_window[i, ich] = fit_subscan(
                    offset[mask], values[mask], beam_sigma, **window_kwargs
                )[0]
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
    if joint:
        diff = _channel_average(amps) - _channel_average(amps_window)
        dflux, derr = combine_amplitudes(diff, statistic=statistic)
        consist["joint-window"] = dflux / derr

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


def beam_fwhm_from_frequency(table, channel, diameter=64.0):
    """Beam FWHM in degrees, as 1.22 lambda / D at the center of the band.

    This is the beam SDTimage uses to choose the pixel size.
    """
    import astropy.constants as c
    import astropy.units as u

    meta = table[channel].meta
    mid_freq = meta["frequency"] + 0.5 * meta["bandwidth"]
    return np.degrees((1.22 * c.c / mid_freq / (diameter * u.m)).to("").value)


def calibrate_crosses(table, channels, calibration):
    """Calibrated data of each channel (Jy), with the relative calibration error.

    The conversion from counts to Jy depends on the elevation of each sample
    (see `srttools.calibration.CalibratorTable.Jy_over_counts`).

    Returns
    -------
    calibrated : dict
        Calibrated data for each channel
    rel_err : float
        Median relative error of the conversion, over all samples and channels
    """
    from .calibration import CalibratorTable

    caltable = CalibratorTable.read(calibration)
    caltable.update()
    calibrated = {}
    rel_errs = []
    for ch in channels:
        elevation = np.asarray(table["el"])[:, get_channel_feed(ch)]
        fc, fce = caltable.Jy_over_counts(channel=ch, map_unit="Jy/beam", elevation=elevation)
        if fc is None:
            raise ValueError(f"No valid calibration for channel {ch}")
        calibrated[ch] = np.asarray(table[ch], dtype=float) * fc
        rel_errs.append(np.asarray(fce / fc))
    return calibrated, float(np.median(np.concatenate([np.ravel(r) for r in rel_errs])))


def _fits_name(label):
    return label.replace(" ", "_")


def _shape_warning(z, nsigma):
    """Message for a joint-minus-window difference of ``z`` standard deviations, or None."""
    if abs(z) <= nsigma:
        return None
    if z < 0:
        return (
            f"The joint fit is lower than the window fit by {-z:.1f} standard deviations: the "
            "source is wider than the beam or offset from the nominal position. Its flux is "
            "probably underestimated; check the profiles."
        )
    return (
        f"The joint fit is higher than the window fit by {z:.1f} standard deviations: the "
        "beam template may be wider than the real beam (see --beam-fwhm), or there is "
        "structure near the source, which affects the window fit more. Check the profiles."
    )


def write_cross_results(result, outroot, header):
    """Write ``<outroot>_cross.fits`` and ``<outroot>_cross.ecsv``.

    The FITS file has the main results as keywords in the primary header, and
    the tables ``SUMMARY``, ``PROFILES`` and ``SUBSCANS``. The ECSV file has
    the summary table, with the same keywords in its metadata.
    """
    from astropy.io import fits

    summary = result["summary"]
    total = summary[summary["label"] == "TOTAL"][0]
    header = dict(header)
    header.update(
        FLUX=float(total["flux"]),
        FLUXERR=float(total["err"]),
        UPLIM=float(total["ul"]),
        DETECTED=bool(total["detected"]),
        NSUB=int(total["nsub"]),
    )
    names = {"RA-Dec": "CONS_ARM", "joint-window": "CONS_SHP"}
    for key, value in result["consistency"].items():
        header[names.get(key, "CONS_POL")] = float(value)
    for ch, bad in result["flagged"].items():
        header[f"RFI_{ch}"[:8]] = len(bad)

    primary = fits.PrimaryHDU()
    for key, value in header.items():
        primary.header[key] = value

    profiles = Table()
    labels = list(result["profiles"].keys())
    profiles["offset"] = result["profiles"][labels[0]][0]
    profiles["offset"].unit = "deg"
    for label in labels:
        _, prof, err = result["profiles"][label]
        profiles[_fits_name(label)] = prof
        profiles[_fits_name(label) + "_err"] = err

    hdus = [primary]
    for name, table in [
        ("SUMMARY", summary),
        ("PROFILES", profiles),
        ("SUBSCANS", result["subscans"]),
    ]:
        hdu = fits.table_to_hdu(table)
        hdu.name = name
        hdus.append(hdu)
    fits.HDUList(hdus).writeto(outroot + "_cross.fits", overwrite=True)

    out = Table(summary, copy=True)
    out.meta.update(header)
    out.meta["consistency"] = {k: float(v) for k, v in result["consistency"].items()}
    out.write(outroot + "_cross.ecsv", overwrite=True)


def plot_cross_results(result, outroot, beam_fwhm, unit="", title=""):
    """Plot the profiles, with the fitted beam, in ``<outroot>_cross.png``.

    One panel per channel and arm, one per arm with all channels, and the
    total in ``<outroot>_cross_total.png``.
    """
    import matplotlib.pyplot as plt

    factor, unit_label = (1e3, "mJy") if unit == "Jy" else (1.0, unit or "counts")
    summary = result["summary"]
    beam_sigma = beam_fwhm / (2 * np.sqrt(2 * np.log(2)))

    def panel(ax, label):
        centers, prof, err = result["profiles"][label]
        row = summary[summary["label"] == label][0]
        arcmin = centers * 60
        ax.errorbar(arcmin, prof * factor, err * factor, fmt="o", ms=2, lw=0.6, color="k")
        model = row["flux"] * np.exp(-(centers**2) / (2 * beam_sigma**2))
        ax.plot(arcmin, model * factor, color="C3", lw=1)
        ax.axhline(0, color="grey", lw=0.5)
        status = "detected" if row["detected"] else f"UL {row['ul'] * factor:.3g}"
        ax.set_title(
            f"{label}: {row['flux'] * factor:.3g} ± {row['err'] * factor:.2g} "
            f"{unit_label} ({status})",
            fontsize=7,
        )
        ax.tick_params(labelsize=7)

    with plt.rc_context({"font.size": 7}):
        channels = [ch for ch in dict.fromkeys(summary["channel"]) if ch != ""]
        rows = [*channels, ""]
        fig, axes = plt.subplots(
            len(rows),
            2,
            figsize=(7, 1.6 * len(rows) + 0.4),
            sharex=True,
            squeeze=False,
            layout="constrained",
        )
        for irow, ch in enumerate(rows):
            for icol, arm in enumerate(["RA", "Dec"]):
                panel(axes[irow, icol], f"{ch} {arm}".strip())
            axes[irow, 0].set_ylabel(f"Flux [{unit_label}]")
        for ax in axes[-1]:
            ax.set_xlabel("Offset along the arm [arcmin]")
        fig.suptitle(title, fontsize=7)
        fig.savefig(outroot + "_cross.png", dpi=200)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(3.5, 2.65), layout="constrained")
        panel(ax, "TOTAL")
        ax.set_xlabel("Offset along the arm [arcmin]")
        ax.set_ylabel(f"Flux [{unit_label}]")
        fig.suptitle(title, fontsize=7)
        fig.savefig(outroot + "_cross_total.png", dpi=200)
        plt.close(fig)


def main_cross(args=None):
    """Flux or upper limit of a point source from cross scans."""
    import argparse
    import os

    description = (
        "Measure the flux (or upper limit) of a point source from cross scans. "
        "In each subscan, a smooth baseline is removed and the amplitude of a "
        "Gaussian beam at the source position is measured; the flux is the "
        "average over subscans and its error comes from the scatter between "
        "subscans."
    )
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "files",
        nargs="*",
        help="Intermediate ScanSet files from SDTimage (*_dump.hdf5)",
        default=[],
    )
    parser.add_argument(
        "-c", "--config", type=str, default=None, help="Config file (instead of files)"
    )
    parser.add_argument(
        "--calibrate",
        type=str,
        default=None,
        help="Calibration file from SDTcal. Without it, the calibrated columns "
        "(<channel>_cal) of the input files are used if present, otherwise counts",
    )
    parser.add_argument(
        "--cal-rel-err",
        type=float,
        default=None,
        help="Relative calibration error, added in quadrature to upper limits. "
        "Default: from the calibration file if given, otherwise 0",
    )
    parser.add_argument(
        "--source-coords",
        type=float,
        nargs=2,
        default=None,
        metavar=("RA", "DEC"),
        help="Position (deg) where the flux is measured. Default: the target "
        "position in the data",
    )
    parser.add_argument(
        "--beam-fwhm",
        type=float,
        default=None,
        help="Beam FWHM in arcmin. Default: 1.22 lambda / D at the center of the band",
    )
    parser.add_argument(
        "--nsigma",
        type=float,
        default=3.0,
        help="Detection threshold and width of the upper limit, in standard deviations",
    )
    parser.add_argument(
        "--rfi-nsigma",
        type=float,
        default=5.0,
        help="Ignore subscans whose noise exceeds the median of all subscans by "
        "more than this many robust standard deviations (e.g. residual RFI). "
        "0 disables the rejection",
    )
    parser.add_argument(
        "--statistic",
        choices=["mean", "median"],
        default="mean",
        help="How the subscan amplitudes are combined",
    )
    parser.add_argument(
        "--chans",
        type=str,
        default=None,
        help="Comma-separated channels to use (default: all, except Stokes Q and U)",
    )
    parser.add_argument(
        "-o", "--outroot", type=str, default=None, help="Root of the output file names"
    )
    parser.add_argument("--noplot", action="store_true", default=False, help="Do not plot")
    args = parser.parse_args(args)

    from .imager import ScanSet

    if args.config is not None:
        table = ScanSet(args.config, norefilt=True, nosub=True, plot=False)
        default_root = args.config.replace(".ini", "")
    elif len(args.files) > 0:
        table = ScanSet(args.files[0] if len(args.files) == 1 else args.files)
        default_root = args.files[0].replace("_dump.hdf5", "").replace(".hdf5", "")
    else:
        raise ValueError("Specify input files or a config file")
    outroot = args.outroot if args.outroot is not None else default_root
    outdir = os.path.dirname(outroot)
    if outdir != "":
        os.makedirs(outdir, exist_ok=True)

    if args.chans is not None:
        channels = args.chans.split(",")
    else:
        channels = [ch for ch in table.chan_columns if "Q" not in ch and "U" not in ch]

    cal_rel_err = 0.0
    if args.calibrate is not None:
        calibrated, cal_rel_err = calibrate_crosses(table, channels, args.calibrate)
        for ch in channels:
            table[ch + "_sdtcross"] = calibrated[ch]
        suffix, unit = "_sdtcross", "Jy"
    elif all(ch + "_cal" in table.colnames for ch in channels):
        suffix, unit = "_cal", "Jy"
    else:
        warnings.warn("No calibration: the flux is in counts.")
        suffix, unit = "", "counts"
    if args.cal_rel_err is not None:
        cal_rel_err = args.cal_rel_err

    if args.source_coords is not None:
        ra0, dec0 = np.radians(args.source_coords)
    else:
        ra0, dec0 = (float(table.meta[k].to("rad").value) for k in ["RA", "Dec"])

    if args.beam_fwhm is not None:
        beam_fwhm = args.beam_fwhm / 60
    else:
        beam_fwhm = beam_fwhm_from_frequency(table, channels[0])

    result = measure_crosses(
        table,
        ra0,
        dec0,
        beam_fwhm,
        channels,
        suffix=suffix,
        rfi_nsigma=args.rfi_nsigma if args.rfi_nsigma > 0 else None,
        statistic=args.statistic,
        nsigma=args.nsigma,
        cal_rel_err=cal_rel_err,
    )
    for key, value in result["consistency"].items():
        if np.abs(value) <= args.nsigma:
            continue
        if key == "joint-window":
            warnings.warn(_shape_warning(value, args.nsigma))
        else:
            warnings.warn(f"{key}: the two estimates differ by {value:.1f} standard deviations")

    source = str(table.meta.get("SOURCE", ""))
    header = dict(
        SOURCE=source,
        RA=float(np.degrees(ra0)),
        DEC=float(np.degrees(dec0)),
        BEAMFWHM=float(beam_fwhm),
        BUNIT=unit,
        NSIGMA=float(args.nsigma),
        CALERR=float(cal_rel_err),
        RFINSIG=float(args.rfi_nsigma),
        COMBINE=args.statistic,
    )
    write_cross_results(result, outroot, header)
    if not args.noplot:
        plot_cross_results(result, outroot, beam_fwhm, unit=unit, title=source)

    summary = Table(result["summary"], copy=True)
    for col in ["flux", "err", "ul"]:
        summary[col].format = ".4g"
    summary.pprint_all()
    return result
