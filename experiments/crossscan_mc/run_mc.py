"""Run the Monte Carlo grid for the three current cross-scan methods.

Usage (from this directory, with srttools from this worktree)::

    PYTHONPATH=../..:. python run_mc.py [--ntrials 200] [--nproc 10] [-o mc_results.ecsv]

Each trial uses its own seed, derived from the grid cell and trial number, so
that any single trial can be reproduced.
"""

import argparse
import itertools
import logging
import os
import sys
import warnings
from multiprocessing import Pool

from astropy.table import Table

SNRS = [0, 1, 2, 3, 5, 10, 30, 100]
NOISES = ["white", "1/f"]
OFFSETS = [0.0, 0.3]
DRIFTS = ["flat", "drift"]


def one_trial(args):
    icell, snr, noise, offset, drift, rfi, itrial = args
    logging.getLogger().setLevel(logging.ERROR)
    warnings.simplefilter("ignore")
    # FITPACK prints warnings on stdout
    sys.stdout = open(os.devnull, "w")
    from mc_sim import ideal_error, simulate_observation
    from methods import run_all

    seed = [icell, itrial]
    sig = ideal_error()
    obs = simulate_observation(
        snr * sig, noise=noise, offset_beams=offset, drift=drift, rfi=rfi, rng=seed
    )
    rows = []
    for method, res in run_all(obs).items():
        rows.append(
            dict(
                snr=snr,
                noise=noise,
                offset=offset,
                drift=drift,
                rfi=rfi,
                trial=itrial,
                method=method,
                # All in units of the ideal error
                flux=res["flux"] / sig,
                err=res["err"] / sig,
                ul=res["ul"] / sig,
                detected=bool(res["detected"]),
                nfail=res["nfail"],
                ncrash=res.get("ncrash", 0),
            )
        )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ntrials", type=int, default=200)
    parser.add_argument("--nproc", type=int, default=os.cpu_count())
    parser.add_argument("-o", "--output", default="mc_results.ecsv")
    parser.add_argument(
        "--rfi",
        action="store_true",
        help="Only run the RFI cells (1/f noise, drift, no offset); seeds "
        "differ from the main grid",
    )
    args = parser.parse_args()

    if args.rfi:
        cells = [(snr, "1/f", 0.0, "drift", True) for snr in SNRS]
        first_cell = 1000
    else:
        cells = [(*c, False) for c in itertools.product(SNRS, NOISES, OFFSETS, DRIFTS)]
        first_cell = 0
    tasks = [
        (first_cell + icell, *cell, itrial)
        for icell, cell in enumerate(cells)
        for itrial in range(args.ntrials)
    ]
    rows = []
    with Pool(args.nproc) as pool:
        for i, res in enumerate(pool.imap_unordered(one_trial, tasks, chunksize=4)):
            rows.extend(res)
            if (i + 1) % 500 == 0:
                print(f"{i + 1}/{len(tasks)}", flush=True)
                Table(rows=rows).write(args.output, overwrite=True)
    Table(rows=rows).write(args.output, overwrite=True)


if __name__ == "__main__":
    main()
