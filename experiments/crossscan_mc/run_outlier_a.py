"""Choose the threshold of the outlier rejection of method (a) (``outlier_nsigma``).

Fits each observation once (``methods.fit_table_a``), then applies
``CalibratorTable.flag_outliers`` with several thresholds
(``methods.method_a_flagged``). Same observations (seeds) as
``run_rfi_threshold.py``.

Usage (from this directory, with srttools from this worktree)::

    PYTHONPATH=../..:. python run_outlier_a.py [--ntrials 100] [-o mc_outlier_a.ecsv]
    PYTHONPATH=../..:. python run_outlier_a.py --bright-rfi -o mc_outlier_a_bright_rfi.ecsv

Output columns as ``run_rfi_threshold.py`` (``method`` is always "a").
"""

import argparse
import logging
import os
import sys
import warnings
from multiprocessing import Pool

from run_rfi_threshold import cells

from astropy.table import Table

THRESHOLDS = [None, 3, 5, 7, 10]


def one_trial(args):
    icell, snr, noise, offset, drift, rfi, itrial = args
    logging.getLogger().setLevel(logging.ERROR)
    warnings.simplefilter("ignore")
    sys.stdout = open(os.devnull, "w")
    from mc_sim import CHANS, ideal_error, simulate_observation
    from methods import fit_table_a, method_a_flagged

    sig = ideal_error()
    obs = simulate_observation(
        snr * sig, noise=noise, offset_beams=offset, drift=drift, rfi=rfi, rng=[icell, itrial]
    )
    truth = {(ch, i) for i, s in enumerate(obs) for ch in CHANS if s[ch + "_has_rfi"]}
    table, nfail, ncrash = fit_table_a(obs)
    fitted = {(str(c), int(str(f)[3:])) for c, f in zip(table["Chan"], table["File"])}
    not_fitted = {(ch, i) for i in range(len(obs)) for ch in CHANS} - fitted
    rows = []
    for nsigma in THRESHOLDS:
        r = method_a_flagged(table, nfail, ncrash, outlier_nsigma=nsigma)
        flagged = set(r["flagged"])
        rows.append(
            dict(
                snr=snr,
                noise=noise,
                offset=offset,
                drift=drift,
                rfi=rfi,
                trial=itrial,
                method="a",
                nsigma=-1.0 if nsigma is None else float(nsigma),
                flux=r["flux"] / sig,
                err=r["err"] / sig,
                ul=r["ul"],
                detected=bool(r["detected"]),
                nfail=nfail,
                ncrash=ncrash,
                nflag=len(flagged),
                nrfi=len(truth),
                # Subscans without a fit (failed) count as caught
                nrfi_missed=len(truth - flagged - not_fitted),
            )
        )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ntrials", type=int, default=100)
    parser.add_argument("--nproc", type=int, default=os.cpu_count())
    parser.add_argument("-o", "--output", default="mc_outlier_a.ecsv")
    parser.add_argument(
        "--bright-rfi",
        action="store_true",
        help="Only run a bright source (S/N 3000, 0.3 beam offset) with RFI (cell 2001)",
    )
    args = parser.parse_args()
    run_cells = cells()
    if args.bright_rfi:
        run_cells = [(2001, 3000, "1/f", 0.3, "drift", True)]
    tasks = [(*c, itrial) for c in run_cells for itrial in range(args.ntrials)]
    rows = []
    with Pool(args.nproc) as pool:
        for i, res in enumerate(pool.imap_unordered(one_trial, tasks, chunksize=2)):
            rows.extend(res)
            if (i + 1) % 100 == 0:
                print(f"{i + 1}/{len(tasks)}", flush=True)
    Table(rows=rows).write(args.output, overwrite=True)


if __name__ == "__main__":
    main()
