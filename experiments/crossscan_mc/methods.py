"""The three current cross-scan methods, applied to simulated observations.

Each ``method_*`` function takes the output of
``mc_sim.simulate_observation`` and returns a dict with:

- ``flux``: flux estimate at the nominal source position;
- ``err``: the error that the method gives for it;
- ``ul``: the 3-sigma upper limit the method gives (NaN if none);
- ``detected``: whether the method's own numbers indicate a detection
  (detection criteria are defined in ``NOTES.md``, as none of the methods
  has one);
- ``nfail``: number of failed fits (dropped data);
- ``ncrash``: (a) only, fits raising an exception, which would stop SDTcal.

(a) ``fit_baseline_plus_bell`` on each subscan, as in ``_treat_scan``
    (``SDTcal``), averaged as in ``CalibratorTable.calculate_src_flux``
    (``SDTlcurve``). The two channels are then averaged.
(b) Matteo's ``plot_all_cross_scan_results.py``: the computing part is copied
    verbatim; only file reading and plotting are replaced.
(c) ``ScanSet.calculate_avg_cross``, called directly on an in-memory table.
"""

import warnings

import numpy as np
from mc_sim import CHANS, CRPIX, DEC0, NPIX, PIX_DEG, RA0, SRC_PIX
from numba import njit
from scipy.interpolate import make_splrep

from astropy.table import Table
from srttools.calibration import CalibratorTable
from srttools.fit import fit_baseline_plus_bell
from srttools.imager import ScanSet


# ----------------------------------------------------------------- (a)
def _arm_coord_deg(sub):
    """Sky coordinate along the arm, as ``scantype`` returns it (deg)."""
    off = (sub["x"] - SRC_PIX) * PIX_DEG
    if sub["arm"] == "RA":
        return RA0 + off / np.cos(np.radians(DEC0))
    return DEC0 + off


def fit_table_a(subscans):
    """Fits of method (a) as a ``CalibratorTable`` (one observation, one source).

    Also returns the number of failed fits (no covariance) and of crashes.
    """
    rows = []
    nfail = ncrash = 0
    for ch in CHANS:
        for isub, sub in enumerate(subscans):
            x = _arm_coord_deg(sub)
            y = sub[ch]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                try:
                    model, fit_info = fit_baseline_plus_bell(x, y, kind="gauss")
                except Exception:
                    ncrash += 1
                    continue
            try:
                uncert = fit_info["param_cov"].diagonal() ** 0.5
            except Exception:
                nfail += 1
                continue
            index = model.param_names.index("amplitude_1")
            rows.append(
                {
                    "Dir": "obs",
                    "File": f"sub{isub}",
                    "Source": "Src",
                    "Chan": ch,
                    "Scan Type": sub["arm"],
                    "Counts": model.amplitude_1.value,
                    "Counts Err": uncert[index],
                    "Width": model.stddev_1.value,
                    "Resid Std": np.std(y - model(x)),
                }
            )
    return CalibratorTable(Table(rows=rows)), nfail, ncrash


def method_a_flagged(table, nfail, ncrash, outlier_nsigma=None):
    """(a) from a ``fit_table_a`` table, optionally rejecting outlier fits."""
    table = CalibratorTable(table, copy=True)
    if outlier_nsigma is not None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            table.flag_outliers(nsigma=outlier_nsigma)
    fluxes, errs = [], []
    for ch in CHANS:
        good = (table["Chan"] == ch) & ~table["Outlier"]
        counts = np.asarray(table["Counts"][good])
        counts_err = np.asarray(table["Counts Err"][good])
        fluxes.append(np.mean(counts))
        errs.append(np.sqrt(np.mean(counts_err**2)))
    flux = np.mean(fluxes)
    err = np.sqrt(np.sum(np.array(errs) ** 2)) / len(errs)
    return dict(
        flux=flux,
        err=err,
        ul=np.nan,
        detected=flux > 3 * err,
        nfail=nfail,
        ncrash=ncrash,
        flagged=sorted(
            (str(c), int(str(f)[3:]))
            for c, f in zip(table["Chan"][table["Outlier"]], table["File"][table["Outlier"]])
        ),
    )


def method_a(subscans):
    fluxes, errs = [], []
    nfail = ncrash = 0
    for ch in CHANS:
        counts, counts_err = [], []
        for sub in subscans:
            x = _arm_coord_deg(sub)
            y = sub[ch]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                try:
                    model, fit_info = fit_baseline_plus_bell(x, y, kind="gauss")
                except Exception:
                    # In _treat_scan this exception is not caught: SDTcal
                    # would stop. Here we count it and drop the subscan.
                    ncrash += 1
                    continue
            # As in _treat_scan: if the covariance is missing, the scan is dropped
            try:
                uncert = fit_info["param_cov"].diagonal() ** 0.5
            except Exception:
                nfail += 1
                continue
            index = model.param_names.index("amplitude_1")
            counts.append(model.amplitude_1.value)
            counts_err.append(uncert[index])
        counts = np.array(counts)
        counts_err = np.array(counts_err)
        # calculate_src_flux with a conversion factor of 1 and no error on it:
        # calculated_flux_err = (counts_err / counts) * counts
        fluxes.append(np.mean(counts))
        errs.append(np.sqrt(np.mean(counts_err**2)))
    flux = np.mean(fluxes)
    err = np.sqrt(np.sum(np.array(errs) ** 2)) / len(errs)
    return dict(flux=flux, err=err, ul=np.nan, detected=flux > 3 * err, nfail=nfail, ncrash=ncrash)


# ----------------------------------------------------------------- tables
def to_dump_table(subscans, cal=True):
    """Build a table like the ``_dump.hdf5`` written by SDTimage."""
    cols = {k: [] for k in ["time", "x", "y", "direction", "Scan_id"] + CHANS}
    for i, sub in enumerate(subscans):
        n = sub["x"].size
        isra = sub["arm"] == "RA"
        cols["time"].append(sub["time"])
        cols["x"].append(sub["x"] if isra else np.full(n, SRC_PIX))
        cols["y"].append(np.full(n, SRC_PIX) if isra else sub["x"])
        cols["direction"].append(np.full(n, isra))
        cols["Scan_id"].append(np.full(n, i))
        for ch in CHANS:
            cols[ch].append(sub[ch])
    t = Table()
    for k, v in cols.items():
        v = np.concatenate(v)
        if k in ["x", "y"]:
            v = v[:, np.newaxis]
        t[k] = v
    if cal:
        for ch in CHANS:
            t[ch + "_cal"] = t[ch]
    return t


# ----------------------------------------------------------------- (b)
# Verbatim from plot_all_cross_scan_results.py
def get_data_from_direction(table, direction, column):
    good = table["direction"] == direction
    if direction:
        x = table["x"][good]
    else:
        x = table["y"][good]
    if len(np.shape(x)) == 2:
        x = x[:, 0]
    flux = table[column][good]
    return x, flux


@njit
def find_trend_change(data, nmax=-1):
    n = data.size
    if nmax < 0:
        nmax = n
    trend_edges = []
    previous_sign = np.sign(data[1] - data[0])
    previous_data = data[0]
    nchanges = 0
    previous_change = -1
    for i, d in enumerate(data):
        if i == 0:
            continue

        local_sign = np.sign(d - previous_data)
        if i == previous_change + 1:
            previous_sign = local_sign
            continue
        previous_data = d

        if local_sign != previous_sign:
            trend_edges.append(i)
            previous_sign = local_sign
            previous_change = i
            nchanges += 1
            if nchanges > nmax:
                break

    return trend_edges


def script_b(full_table, extnames, crpix, skip_nan=False):
    """``load_and_plot_data`` without file I/O and plots.

    Returns the TOTAL detrended profile, its error, the bin centers, the
    printed upper limit, and the avoid region of the last subscan.

    ``skip_nan=True`` applies the patch in ``script_skip_nan.patch``: samples
    flagged by SDTimage (NaN in the ``_cal`` columns) are dropped.
    """
    xmin = min(full_table["x"].min(), full_table["y"].min())
    xmax = max(full_table["x"].max(), full_table["y"].max())

    total_flux_sub = 0.0
    total_flux_sq = 0.0
    total_expo = 0.0
    for extname in extnames:
        extname = extname.replace("X_", "")
        direction = extname.split("_")[-1]
        ch = extname.replace("_" + direction, "").replace("FEED", "Feed")
        direction = 0 if direction.lower() == "horizontal" else 1

        x, fluxes = get_data_from_direction(full_table, direction, ch + "_cal")
        x = np.asarray(x)
        fluxes = np.asarray(fluxes)
        if skip_nan:
            finite = np.isfinite(fluxes)
            x, fluxes = x[finite], fluxes[finite]

        x = x - crpix
        direction_change_edges = find_trend_change(x)
        histogram_xbins = np.arange(np.floor(xmin) - crpix, np.ceil(xmax) + 1 - crpix)
        ave_sub = np.zeros(histogram_xbins.size - 1)
        ave_sub_sq = np.zeros(histogram_xbins.size - 1)
        edges = np.concatenate(([0], direction_change_edges, [x.size]))
        expo = np.zeros(histogram_xbins.size - 1)

        for e0, e1 in zip(edges[:-1], edges[1:]):  # noqa: RUF007 (verbatim)
            local_data = fluxes[e0 + 1 : e1 - 1].flatten()
            local_x = x[e0 + 1 : e1 - 1].flatten()
            order = np.argsort(local_x)
            local_x = local_x[order]
            local_data = local_data[order]
            knots_position = np.concatenate(
                (
                    [-0.001, -0.001, -0.001],
                    np.linspace(0, 0.4, 3)[1:],
                    np.linspace(0.60, 1, 3)[:-1],
                    [1.001, 1.001, 1.001],
                )
            )

            def percentage_to_x(p):
                return local_x.min() + p * (local_x.max() - local_x.min())

            avoid_region = percentage_to_x(np.array([0.35, 0.65]))
            knots = percentage_to_x(np.array(knots_position))

            good_x = (local_x < avoid_region[0]) | (local_x > avoid_region[1])

            spline = make_splrep(local_x[good_x], local_data[good_x], k=2, t=knots, s=1)
            local_data_sub = local_data - spline(local_x)

            ave_sub += np.histogram(local_x, bins=histogram_xbins, weights=local_data_sub)[0]
            ave_sub_sq += np.histogram(local_x, bins=histogram_xbins, weights=local_data_sub**2)[0]
            expo += np.histogram(local_x, bins=histogram_xbins)[0]
        plotting_xbins = histogram_xbins[:-1] + 0.5
        total_flux_sub += ave_sub
        total_flux_sq += ave_sub_sq
        total_expo += expo

    def limit_to_zero(array):
        return np.where(array < 0, 0, array)

    total_expo[total_expo == 0] = 1
    total_flux_sub /= total_expo
    total_flux_sdev = np.sqrt(total_flux_sq / total_expo - limit_to_zero(total_flux_sub) ** 2)
    total_flux_mean_err = total_flux_sdev / np.sqrt(total_expo)
    good_x_values = (plotting_xbins < avoid_region[0]) | (plotting_xbins > avoid_region[1])
    ul = np.nanmedian(total_flux_mean_err[good_x_values]) * 3
    return plotting_xbins, total_flux_sub, total_flux_mean_err, ul, avoid_region


EXTNAMES = [f"X_{ch.upper()}_{d}" for ch in CHANS for d in ["HORIZONTAL", "VERTICAL"]]


def method_b(subscans, table=None, skip_nan=False):
    if table is None:
        table = to_dump_table(subscans)
    xb, prof, prof_err, ul, _ = script_b(table, EXTNAMES, CRPIX, skip_nan=skip_nan)
    # Bin labelled 0 ("0 is the center" in the plots), i.e. [-0.5, 0.5)
    # around CRPIX
    i0 = np.argmin(np.abs(xb))
    flux, err = prof[i0], prof_err[i0]
    # In the TOTAL plot, a source is seen when the profile goes above the
    # 3-sigma line
    return dict(flux=flux, err=err, ul=ul, detected=flux > ul, nfail=0)


# ----------------------------------------------------------------- (c)
def scanset_from_table(table):
    t = Table(table, copy=True)
    for ch in CHANS:
        if ch + "_cal" in t.colnames:
            t.remove_column(ch + "_cal")
        t[ch + "-filt"] = np.ones(len(t), dtype=bool)
    t["el"] = np.full((len(t), 1), np.radians(45))
    t["x"].meta["frame"] = "icrs"
    t.meta["npix"] = np.array([NPIX, NPIX])
    return ScanSet(t)


def method_c(subscans, table=None, return_scanset=False, rfi_nsigma=None):
    """``rfi_nsigma=None``: the method as it was before RFI rejection."""
    if table is None:
        table = to_dump_table(subscans, cal=False)
    scanset = scanset_from_table(table)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        crosses = scanset.calculate_avg_cross(rfi_nsigma=rfi_nsigma)
    # Bin containing the nominal source position (pixel coordinate SRC_PIX)
    isrc = int(np.floor(SRC_PIX))
    vals, sdevs = [], []
    for ch in CHANS:
        for direction in [0, 1]:
            vals.append(crosses[f"{ch}-{direction}"][isrc])
            sdevs.append(crosses[f"{ch}-{direction}-Sdev"])
    sdevs = np.array(sdevs)
    comb_err = np.sqrt(np.sum(sdevs**2, axis=0)) / len(vals)
    flux = np.mean(vals)
    err = comb_err[isrc]
    # 3 x Sdev away from the source, as the black "3 sigma u.l." line in
    # Matteo's plots (median over bins more than 6 pixels from the source)
    centers = np.arange(NPIX) + 0.5
    off = np.abs(centers - SRC_PIX) > 6
    ul = 3 * np.median(comb_err[off])
    res = dict(flux=flux, err=err, ul=ul, detected=flux > ul, nfail=0)
    if return_scanset:
        return res, scanset
    return res


def run_all(subscans):
    """Run the three methods on one observation.

    As in SDTimage, ``calculate_avg_cross`` removes spikes from the channel
    data in place, and the ``_cal`` columns read by the script are written
    afterwards: so (b) is fed the data as modified by (c).
    """
    res = {"a": method_a(subscans)}
    res.update(run_bc(subscans))
    return res


def run_bc(subscans, rfi_nsigma=None):
    """Methods (b) and (c) only, optionally with RFI rejection.

    Flagged samples are passed to (b) as NaN, as SDTimage writes them in the
    ``_cal`` columns once RFI rejection is on, and (b) is the script patched
    to skip NaNs (``method_b(..., skip_nan=True)``).
    """
    res = {}
    res["c"], scanset = method_c(subscans, return_scanset=True, rfi_nsigma=rfi_nsigma)
    table = to_dump_table(subscans)
    for ch in CHANS:
        cal = np.array(scanset[ch], dtype=float)
        cal[~np.asarray(scanset[ch + "-filt"])] = np.nan
        table[ch + "_cal"] = cal
    res["b"] = method_b(subscans, table=table, skip_nan=rfi_nsigma is not None)
    res["flagged"] = {
        ch: sorted(set(np.asarray(scanset["Scan_id"])[~np.asarray(scanset[ch + "-filt"])]))
        for ch in CHANS
    }
    return res
