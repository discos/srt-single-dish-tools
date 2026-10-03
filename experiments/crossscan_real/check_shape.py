"""Joint-vs-window shape check of SDTcross (``CONS_SHP``) on the real data.

For each observation: the check on the original data, and on sign-flip null
data (source subtracted, subscan signs flipped, as in ``validate_real.py``),
where it should be distributed as a standard normal.

Usage (from this directory)::

    PYTHONPATH=../..:../crossscan_mc python check_shape.py [--npatterns 20]

Writes ``shape_check.ecsv``.
"""

import argparse
from multiprocessing import Pool

import numpy as np
import validate_real as v

from astropy.table import Table


def one(args):
    name, npatterns = args
    v.quiet()
    table, _, ra0, dec0, beam = v.load(name)
    d = v.run_d(table, ra0, dec0, beam)
    rows = [
        dict(obs=name, test="original", pattern=-1, z=d["result"]["consistency"]["joint-window"])
    ]
    groups = v.subscan_groups(table, ra0, dec0)
    for i in range(npatterns):
        signs = v.balanced_signs(groups, np.random.default_rng([v.OBSERVATIONS.index(name), i]))
        new = v.modified_table(table, signs, ra0, dec0, beam, source=d["flux"])
        z = v.run_d(new, ra0, dec0, beam)["result"]["consistency"]["joint-window"]
        rows.append(dict(obs=name, test="null", pattern=i, z=z))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--npatterns", type=int, default=20)
    args = parser.parse_args()
    v.quiet()
    with Pool(len(v.OBSERVATIONS)) as pool:
        rows = [
            r for obs in pool.map(one, [(o, args.npatterns) for o in v.OBSERVATIONS]) for r in obs
        ]
    t = Table(rows=rows)
    t.write("shape_check.ecsv", overwrite=True)
    for obs in v.OBSERVATIONS:
        g = t[t["obs"] == obs]
        null = g[g["test"] == "null"]["z"]
        orig = g[g["test"] == "original"]["z"][0]
        print(
            f"{obs:18s} original {orig:+6.2f}   null: mean {np.mean(null):+.2f} "
            f"std {np.std(null):.2f} max|z| {np.max(np.abs(null)):.2f}"
        )


if __name__ == "__main__":
    main()
