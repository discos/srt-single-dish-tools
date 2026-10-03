"""Stage B: choose the settings of the new method (d), ``srttools.crossscan``.

Method (d), for each observation:

1. subscans with anomalous noise are rejected, channel by channel, with
   ``find_outlier_subscans`` (``rfi_nsigma=5``), as SDTimage does;
2. each subscan and channel gets a forced Gaussian amplitude at the nominal
   source position over a spline baseline (``crossscan.fit_subscan``);
3. the two channels are averaged subscan by subscan (so that the noise they
   share is not counted twice), then the subscans are combined with
   ``crossscan.combine_amplitudes`` (mean or median), whose error comes from
   the scatter between subscans.

Variants: baseline fitted outside a window (``joint=0``) for windows of 1,
1.5, 2 beam FWHM, or together with the Gaussian (``joint=1``, with or without
a knot-free window of 1.5 FWHM); knot spacing 2, 3, 5 FWHM.

Method (b) with RFI rejection (as in ``run_rfi_threshold.py``, nsigma=5) is
run on the same observations, for comparison. Seeds are those of
``run_mc.py``, so the observations are also the same as in Stage A.

Usage::

    PYTHONPATH=../..:. python run_crossscan_mc.py [--ntrials 200] [-o mc_crossscan.ecsv]

Output: one row per trial and method. Fluxes in units of the ideal error.
"""

import argparse
import itertools
import logging
import os
import sys
import warnings
from multiprocessing import Pool

from run_mc import DRIFTS, NOISES, OFFSETS, SNRS

from astropy.table import Table

TEST_SNRS = [0, 3, 5, 10, 30]
WINDOWS = [1.0, 1.5, 2.0]  # beam FWHM
SPACINGS = [2.0, 3.0, 5.0]  # beam FWHM
RFI_NSIGMA = 5.0


def variants():
    """(name, joint, window, knot_spacing), window and spacing in beam FWHM."""
    out = []
    for w, k in itertools.product(WINDOWS, SPACINGS):
        out.append((f"d_win{w}_k{k}", False, w, k))
    for w, k in itertools.product([0.0, 1.5], SPACINGS):
        out.append((f"d_joint_win{w}_k{k}", True, w, k))
    return out


def cells():
    """(cell id, snr, noise, offset, drift, rfi), with the ids of run_mc.py."""
    grid = list(itertools.product(SNRS, NOISES, OFFSETS, DRIFTS))
    out = []
    for snr in TEST_SNRS:
        for noise, offset, drift in itertools.product(NOISES, OFFSETS, DRIFTS):
            out.append((grid.index((snr, noise, offset, drift)), snr, noise, offset, drift, False))
        out.append((1000 + SNRS.index(snr), snr, "1/f", 0.0, "drift", True))
    return out


def per_subscan_amplitudes(subscans, flagged, joint, window, knot_spacing):
    """Channel-averaged amplitude of each subscan (NaN if all channels flagged)."""
    import numpy as np
    from mc_sim import BEAM_FWHM_PIX, BEAM_SIGMA_PIX, CHANS, SRC_PIX

    from srttools.crossscan import fit_subscan

    amps = np.full((len(subscans), len(CHANS)), np.nan)
    for i, sub in enumerate(subscans):
        for j, ch in enumerate(CHANS):
            if i in flagged[ch]:
                continue
            amps[i, j] = fit_subscan(
                sub["x"],
                sub[ch],
                BEAM_SIGMA_PIX,
                center=SRC_PIX,
                window=window * BEAM_FWHM_PIX,
                knot_spacing=knot_spacing * BEAM_FWHM_PIX,
                joint=joint,
            )[0]
    n = np.sum(np.isfinite(amps), axis=1)
    return np.where(n > 0, np.nansum(amps, axis=1) / np.maximum(n, 1), np.nan)


def one_trial(args):
    icell, snr, noise, offset, drift, rfi, itrial = args
    logging.getLogger().setLevel(logging.ERROR)
    warnings.simplefilter("ignore")
    sys.stdout = open(os.devnull, "w")
    from mc_sim import ideal_error, simulate_observation
    from methods import run_bc

    from srttools.crossscan import combine_amplitudes

    sig = ideal_error()
    obs = simulate_observation(
        snr * sig, noise=noise, offset_beams=offset, drift=drift, rfi=rfi, rng=[icell, itrial]
    )
    common = dict(snr=snr, noise=noise, offset=offset, drift=drift, rfi=rfi, trial=itrial)
    res_bc = run_bc(obs, rfi_nsigma=RFI_NSIGMA)
    flagged = {ch: set(v) for ch, v in res_bc.pop("flagged").items()}
    b = res_bc["b"]
    rows = [
        dict(
            **common,
            method="b_rfi",
            joint=False,
            window=0.0,
            knot_spacing=0.0,
            statistic="",
            flux=b["flux"] / sig,
            err=b["err"] / sig,
            ul=b["ul"] / sig,
            detected=bool(b["detected"]),
        )
    ]
    for name, joint, window, spacing in variants():
        amps = per_subscan_amplitudes(obs, flagged, joint, window, spacing)
        for statistic in ["mean", "median"]:
            flux, err = combine_amplitudes(amps, statistic=statistic)
            rows.append(
                dict(
                    **common,
                    method=f"{name}_{statistic}",
                    joint=joint,
                    window=window,
                    knot_spacing=spacing,
                    statistic=statistic,
                    flux=flux / sig,
                    err=err / sig,
                    ul=(max(flux, 0) + 3 * err) / sig,
                    detected=bool(flux > 3 * err),
                )
            )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ntrials", type=int, default=200)
    parser.add_argument("--nproc", type=int, default=os.cpu_count())
    parser.add_argument("-o", "--output", default="mc_crossscan.ecsv")
    args = parser.parse_args()
    tasks = [(*c, itrial) for c in cells() for itrial in range(args.ntrials)]
    rows = []
    with Pool(args.nproc) as pool:
        for i, res in enumerate(pool.imap_unordered(one_trial, tasks, chunksize=2)):
            rows.extend(res)
            if (i + 1) % 500 == 0:
                print(f"{i + 1}/{len(tasks)}", flush=True)
                Table(rows=rows).write(args.output, overwrite=True)
    Table(rows=rows).write(args.output, overwrite=True)


if __name__ == "__main__":
    main()
