"""Compare SDTcross (d) with the cross-scan script (b) on the real SRT data.

The data in ``/Volumes/Elements/SRT`` are only read. For each observation:

- ``original``: both methods on the data as they are;
- ``null``: sign-flip null tests. The source measured by (d) (flux times the
  beam, at the nominal position) is subtracted from the good samples, then
  each subscan (both channels) is multiplied by +1 or -1, with as many +1 as
  -1 in each arm and direction of motion: what is left of the source cancels,
  the real noise stays. Over many random sign patterns, the scatter of the
  measured flux is the real error of the method, to be compared with the error
  it reports;
- ``inject``: the same sign-flipped data plus a fake point source (beam FWHM
  from the band center, as SDTcross) of known flux at the nominal position:
  recovered flux and detection rate;
- ``blocks``: the subscans are split in ``NBLOCKS`` consecutive blocks in
  time; chi-square of the block fluxes around their weighted mean;
- channels (LCP vs RCP) and arms (RA vs Dec) measured separately.

Method (b) is the computing part of ``plot_all_cross_scan_results.py``
(``experiments/crossscan_mc/methods.py``, verbatim), on the ``_cal`` columns
of the dump files, as Matteo runs it (flagged samples were stored as 0 in
those dumps). Its flux is the TOTAL profile in the bin at CRPIX1 (the one
labelled 0 in its plots), its error is the error bar of that bin, and a
detection is a flux above its "3 sigma u.l." line.

Method (d) is ``srttools.crossscan.measure_crosses`` with the SDTcross
defaults. Its error is the reported one; a detection is flux > 3 err.

Usage (from this directory)::

    PYTHONPATH=../..:../crossscan_mc python validate_real.py [--npatterns 50] [--nproc 10]

Writes ``real_results.ecsv`` (one row per observation, method and test) and
the figures (``plot_real.py``).
"""

import argparse
import logging
import os
import sys
import warnings
from multiprocessing import Pool

import numpy as np

from astropy.io import fits
from astropy.table import Table

DATADIR = "/Volumes/Elements/SRT"
OBSERVATIONS = [
    "HS1630+2355_Obs0",
    "Mrk915_Obs0",
    "SDSSJ094715_Obs0",
    "SDSSJ094715_Obs1",
    "SDSSJ133654_Obs0",
    "SDSSJ133654_Obs1",
]
CHANS = ["Feed0_LCP", "Feed0_RCP"]
NBLOCKS = 8
INJECT_SNR = 5.0  # fake source, in units of the error of (d) on the original data


def quiet():
    logging.getLogger().setLevel(logging.ERROR)
    warnings.simplefilter("ignore")


def load(name):
    """Dump table, CRPIX1 of the script, source position (rad), beam FWHM (deg)."""
    from srttools.crossscan import beam_fwhm_from_frequency

    root = os.path.join(DATADIR, "CCB_SARDARA_" + name)
    table = Table.read(root + "_dump.hdf5")
    with fits.open(root + "_cal.fits") as hdul:
        crpix = next(
            h.header["CRPIX1"] for h in hdul if h.header.get("EXTNAME", "").startswith("X_")
        )
    ra0, dec0 = (float(table.meta[k].to("rad").value) for k in ["RA", "Dec"])
    return table, crpix, ra0, dec0, beam_fwhm_from_frequency(table, CHANS[0])


def extnames(chans=CHANS, arms=("RA", "Dec")):
    # In the script, "HORIZONTAL" selects direction == 0, i.e. the Dec arm
    names = {"RA": "VERTICAL", "Dec": "HORIZONTAL"}
    return [f"X_{ch.upper()}_{names[a]}" for ch in chans for a in arms]


def run_b(table, crpix, chans=CHANS, arms=("RA", "Dec")):
    from methods import script_b

    stdout = sys.stdout
    sys.stdout = open(os.devnull, "w")
    try:
        xb, prof, err, ul, _ = script_b(table, extnames(chans, arms), crpix)
    finally:
        sys.stdout = stdout
    i0 = np.argmin(np.abs(xb))
    return dict(flux=prof[i0], err=err[i0], ul=ul, detected=bool(prof[i0] > ul))


def run_d(table, ra0, dec0, beam, label="TOTAL"):
    from srttools.crossscan import measure_crosses

    res = measure_crosses(table, ra0, dec0, beam, CHANS, suffix="_cal")
    row = res["summary"][res["summary"]["label"] == label][0]
    return dict(
        flux=row["flux"],
        err=row["err"],
        ul=row["ul"],
        detected=bool(row["detected"]),
        result=res,
    )


def subscan_groups(table, ra0, dec0):
    """Scan_id of each subscan, with its arm and direction of motion."""
    from srttools.crossscan import arm_offsets

    offset, arm = arm_offsets(table, 0, ra0, dec0)
    sid = np.asarray(table["Scan_id"])
    ids = np.unique(sid)
    groups = {}
    for i in ids:
        mask = sid == i
        key = (bool(arm[mask][0]), bool(offset[mask][-1] > offset[mask][0]))
        groups.setdefault(key, []).append(i)
    return groups


def balanced_signs(groups, rng):
    signs = {}
    for ids in groups.values():
        s = np.ones(len(ids))
        s[: len(ids) // 2] = -1
        rng.shuffle(s)
        signs.update(zip(ids, s))
    return signs


def modified_table(table, signs, ra0, dec0, beam, source=0.0, inject=0.0):
    """Copy of ``table`` with the source subtracted, sign-flipped ``_cal`` columns, and a fake source added.

    ``source`` is only subtracted from good samples (``<channel>-filt``): the
    flagged ones are 0 in the dumps read by (b) and stay 0.
    """
    from srttools.crossscan import arm_offsets

    new = Table(table, copy=True)
    sid = np.asarray(table["Scan_id"])
    factor = np.array([signs[i] for i in sid])
    sigma = beam / (2 * np.sqrt(2 * np.log(2)))
    offset, _ = arm_offsets(table, 0, ra0, dec0)
    model = np.exp(-(offset**2) / (2 * sigma**2))
    for ch in CHANS:
        good = np.asarray(table[ch + "-filt"], dtype=bool)
        data = np.asarray(table[ch + "_cal"]) - source * model * good
        new[ch + "_cal"] = data * factor + inject * model
    return new


def one_pattern(args):
    name, ipattern, source, inject = args
    quiet()
    table, crpix, ra0, dec0, beam = load(name)
    rng = np.random.default_rng([OBSERVATIONS.index(name), ipattern])
    signs = balanced_signs(subscan_groups(table, ra0, dec0), rng)
    rows = []
    for test, amp in [("null", 0.0), ("inject", inject)]:
        new = modified_table(table, signs, ra0, dec0, beam, source=source, inject=amp)
        for method, res in [("b", run_b(new, crpix)), ("d", run_d(new, ra0, dec0, beam))]:
            rows.append(
                dict(
                    obs=name,
                    method=method,
                    test=test,
                    pattern=ipattern,
                    true=amp,
                    flux=res["flux"],
                    err=res["err"],
                    ul=res["ul"],
                    detected=res["detected"],
                )
            )
    return rows


def blocks(table, nblocks=NBLOCKS):
    """Tables with the subscans of ``nblocks`` consecutive blocks."""
    ids = np.unique(table["Scan_id"])
    sid = np.asarray(table["Scan_id"])
    return [table[np.isin(sid, chunk)] for chunk in np.array_split(ids, nblocks)]


def chi2(fluxes, errs):
    fluxes, errs = np.asarray(fluxes), np.asarray(errs)
    w = 1 / errs**2
    mean = np.sum(w * fluxes) / np.sum(w)
    return np.sum((fluxes - mean) ** 2 * w), len(fluxes) - 1


def one_observation(name):
    """Original data, channels, arms and time blocks, for both methods."""
    quiet()
    table, crpix, ra0, dec0, beam = load(name)
    rows = []

    def add(method, test, res, **extra):
        rows.append(
            dict(
                obs=name,
                method=method,
                test=test,
                pattern=-1,
                true=np.nan,
                flux=res["flux"],
                err=res["err"],
                ul=res["ul"],
                detected=res["detected"],
                **extra,
            )
        )

    d = run_d(table, ra0, dec0, beam)
    add("b", "original", run_b(table, crpix))
    add("d", "original", d)
    summary = d["result"]["summary"]
    for label in CHANS + ["RA", "Dec"]:
        row = summary[summary["label"] == label][0]
        add("d", label, dict(flux=row["flux"], err=row["err"], ul=row["ul"], detected=False))
    for ch in CHANS:
        add("b", ch, run_b(table, crpix, chans=[ch]))
    for arm in ["RA", "Dec"]:
        add("b", arm, run_b(table, crpix, arms=[arm]))
    for iblock, sub in enumerate(blocks(table)):
        add("b", f"block{iblock}", run_b(sub, crpix))
        add("d", f"block{iblock}", run_d(sub, ra0, dec0, beam))
    return rows, dict(obs=name, d_flux=d["flux"], d_err=d["err"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--npatterns", type=int, default=50)
    parser.add_argument("--nproc", type=int, default=os.cpu_count())
    parser.add_argument("-o", "--output", default="real_results.ecsv")
    args = parser.parse_args()
    quiet()

    rows = []
    d_info = {}
    with Pool(args.nproc) as pool:
        for obs_rows, info in pool.imap_unordered(one_observation, OBSERVATIONS):
            rows.extend(obs_rows)
            d_info[info["obs"]] = info
            print(info["obs"], "done", flush=True)
        tasks = [
            (name, i, d_info[name]["d_flux"], INJECT_SNR * d_info[name]["d_err"])
            for name in OBSERVATIONS
            for i in range(args.npatterns)
        ]
        for i, res in enumerate(pool.imap_unordered(one_pattern, tasks)):
            rows.extend(res)
            if (i + 1) % 50 == 0:
                print(f"{i + 1}/{len(tasks)}", flush=True)
    Table(rows=rows).write(args.output, overwrite=True)


if __name__ == "__main__":
    main()
