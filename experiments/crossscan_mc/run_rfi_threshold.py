"""Choose the threshold of the RFI subscan rejection (``rfi_nsigma``).

Runs methods (b) and (c) only (``methods.run_bc``) for several thresholds, on
the same observations as ``run_mc.py`` (same seeds), plus a bright calibrator
with a pointing offset (the case most likely to cause false rejections).

Usage (from this directory, with srttools from this worktree)::

    PYTHONPATH=../..:. python run_rfi_threshold.py [--ntrials 100] [-o mc_rfi_threshold.ecsv]

Output: one row per trial, method and threshold, with the columns of
``run_mc.py`` plus ``nsigma`` (-1 = no rejection), ``nflag`` (subscans
flagged, summed over the two channels), ``nrfi`` (subscans with RFI added),
``nrfi_missed`` (subscans with RFI that were not flagged).
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

THRESHOLDS = [None, 3, 4, 5, 7, 10]
TEST_SNRS = [0, 5, 30]
BRIGHT_SNR = 3000  # source peak ~40 times the per-sample noise


def cells():
    """(cell id, snr, noise, offset, drift, rfi), with the ids of run_mc.py."""
    grid = list(itertools.product(SNRS, NOISES, OFFSETS, DRIFTS))
    out = []
    for snr in TEST_SNRS:
        for noise, drift in [("white", "flat"), ("1/f", "drift")]:
            out.append((grid.index((snr, noise, 0.0, drift)), snr, noise, 0.0, drift, False))
        out.append((1000 + SNRS.index(snr), snr, "1/f", 0.0, "drift", True))
    out.append((2000, BRIGHT_SNR, "1/f", 0.3, "drift", False))
    return out


def one_trial(args):
    icell, snr, noise, offset, drift, rfi, itrial = args
    logging.getLogger().setLevel(logging.ERROR)
    warnings.simplefilter("ignore")
    sys.stdout = open(os.devnull, "w")
    from mc_sim import CHANS, ideal_error, simulate_observation
    from methods import run_bc

    sig = ideal_error()
    obs = simulate_observation(
        snr * sig, noise=noise, offset_beams=offset, drift=drift, rfi=rfi, rng=[icell, itrial]
    )
    truth = {ch: {i for i, s in enumerate(obs) if s[ch + "_has_rfi"]} for ch in CHANS}
    rows = []
    for nsigma in THRESHOLDS:
        res = run_bc(obs, rfi_nsigma=nsigma)
        flagged = res.pop("flagged")
        nflag = sum(len(v) for v in flagged.values())
        nrfi = sum(len(v) for v in truth.values())
        missed = sum(len(truth[ch] - set(flagged[ch])) for ch in CHANS)
        for method, r in res.items():
            rows.append(
                dict(
                    snr=snr,
                    noise=noise,
                    offset=offset,
                    drift=drift,
                    rfi=rfi,
                    trial=itrial,
                    method=method,
                    nsigma=-1.0 if nsigma is None else float(nsigma),
                    flux=r["flux"] / sig,
                    err=r["err"] / sig,
                    ul=r["ul"] / sig,
                    detected=bool(r["detected"]),
                    nfail=0,
                    ncrash=0,
                    nflag=nflag,
                    nrfi=nrfi,
                    nrfi_missed=missed,
                )
            )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ntrials", type=int, default=100)
    parser.add_argument("--nproc", type=int, default=os.cpu_count())
    parser.add_argument("-o", "--output", default="mc_rfi_threshold.ecsv")
    args = parser.parse_args()
    tasks = [(*c, itrial) for c in cells() for itrial in range(args.ntrials)]
    rows = []
    with Pool(args.nproc) as pool:
        for i, res in enumerate(pool.imap_unordered(one_trial, tasks, chunksize=2)):
            rows.extend(res)
            if (i + 1) % 100 == 0:
                print(f"{i + 1}/{len(tasks)}", flush=True)
    Table(rows=rows).write(args.output, overwrite=True)


if __name__ == "__main__":
    main()
