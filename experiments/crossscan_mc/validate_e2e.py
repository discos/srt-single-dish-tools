"""Check that the array-level Monte Carlo reproduces the real srttools pipeline.

One simulated observation is written to FITS files (``save_scan``), one file
per subscan, and loaded with srttools' own code:

1. ``Scan`` (loading + default ALS baseline subtraction): the channel data
   must match the baseline-subtracted arrays used in the Monte Carlo. Since
   (b) is verbatim code and (c) is the real ``calculate_avg_cross``, this
   validates their input.
2. ``CalibratorTable.from_scans`` (the SDTcal path): the fitted amplitudes
   and errors must match those of ``methods.method_a``.

Usage::

    PYTHONPATH=../..:. python validate_e2e.py OUTDIR
"""

import glob
import logging
import os
import sys
import warnings

import numpy as np
from mc_sim import CHANS, DEC0, PIX_DEG, RA0, SRC_PIX, ideal_error, simulate_observation

from srttools.calibration import CalibratorTable
from srttools.scan import Scan
from srttools.simulate import save_scan

logging.getLogger().setLevel(logging.ERROR)
warnings.simplefilter("ignore")


def write_observation(subscans, outdir):
    os.makedirs(outdir, exist_ok=True)
    for i, sub in enumerate(subscans):
        off = (sub["x"] - SRC_PIX) * PIX_DEG
        if sub["arm"] == "RA":
            ra = RA0 + off / np.cos(np.radians(DEC0))
            dec = np.full_like(ra, DEC0)
        else:
            dec = DEC0 + off
            ra = np.full_like(dec, RA0)
        save_scan(
            sub["time"],
            ra,
            dec,
            {"Ch0": sub[CHANS[0] + "_raw"], "Ch1": sub[CHANS[1] + "_raw"]},
            filename=os.path.join(outdir, f"{i:03d}_{sub['arm']}.fits"),
            src_ra=RA0,
            src_dec=DEC0,
            srcname="Fake",
            counts_to_K=1.0,
        )


def per_subscan_fits(subscans):
    """(a) per-subscan amplitudes, as in methods.method_a, without averaging."""
    from methods import _arm_coord_deg

    from srttools.fit import fit_baseline_plus_bell

    out = {ch: [] for ch in CHANS}
    for ch in CHANS:
        for sub in subscans:
            try:
                model, info = fit_baseline_plus_bell(_arm_coord_deg(sub), sub[ch], kind="gauss")
                err = info["param_cov"].diagonal()[model.param_names.index("amplitude_1")] ** 0.5
                out[ch].append((model.amplitude_1.value, err))
            except Exception:
                out[ch].append((np.nan, np.nan))
    return out


def main(outdir):
    sig = ideal_error()
    for case in [dict(noise="white", drift="flat"), dict(noise="1/f", drift="drift")]:
        print(f"--- {case}")
        subscans = simulate_observation(30 * sig, offset_beams=0.0, rng=7, **case)
        d = os.path.join(outdir, case["noise"].replace("/", "") + "_" + case["drift"])
        write_observation(subscans, d)
        files = sorted(glob.glob(os.path.join(d, "*.fits")))

        # 1. Loading and baseline subtraction
        maxdiff = 0
        for f, sub in zip(files, subscans):
            scan = Scan(f, norefilt=True, nosave=True)
            for ich, ch in enumerate(CHANS):
                real = np.asarray(scan[scan.chan_columns()[ich]]).ravel()
                maxdiff = max(maxdiff, np.max(np.abs(real - sub[ch])))
        print(f"Scan loading: max |difference| = {maxdiff / sig:.2e} ideal errors")

        # 2. The SDTcal fits
        table = CalibratorTable()
        table.from_scans(files)
        mine = per_subscan_fits(subscans)
        for ich, ch in enumerate(CHANS):
            chname = sorted(set(table["Chan"]))[ich]
            t = table[table["Chan"] == chname]
            order = np.argsort([os.path.basename(str(n)) for n in t["File"]])
            real = np.array(t["Counts"])[order]
            real_err = np.array(t["Counts Err"])[order]
            arr = np.array(mine[ch])
            good = np.isfinite(arr[:, 0])
            if good.sum() != len(real):
                print(f"{chname}: {len(real)} fits in SDTcal, {good.sum()} here: cannot compare")
                continue
            rel = np.abs(real_err / arr[good, 1] - 1)
            bad = rel > 0.01
            print(
                f"{chname}: {len(real)} fits; max |dAmp| = "
                f"{np.max(np.abs(real - arr[good, 0])) / sig:.2e} ideal errors; "
                f"errors differ by >1% in {bad.sum()} fits. Their errors: SDTcal "
                f"{np.round(real_err[bad] / sig, 1)}, here {np.round(arr[good, 1][bad] / sig, 1)}; "
                f"median error of all fits: {np.median(real_err) / sig:.1f}"
            )


if __name__ == "__main__":
    main(sys.argv[1])
