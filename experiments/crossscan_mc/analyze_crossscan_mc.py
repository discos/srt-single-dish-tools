"""Summary and figures for ``run_crossscan_mc.py``.

Usage::

    python analyze_crossscan_mc.py mc_crossscan.ecsv [--best d_joint_win0.0_k3.0_mean]

Writes ``summary_crossscan.ecsv`` (all methods and cells),
``ranking_crossscan.ecsv`` (one row per method, worst cases over the grid)
and ``fig_crossscan.png`` (``--best`` compared with (b) with RFI rejection).
Fluxes in units of the ideal error.

For each method and grid cell:

- ``bias``: mean of (flux - true);
- ``scatter``: standard deviation of the flux;
- ``pull_std``: standard deviation of (flux - true) / err; 1 if the errors
  are right;
- ``cov1``: fraction of trials with the truth within the 1-sigma error bar
  (0.683 if right);
- ``detected``: fraction of trials with a detection;
- ``ul_cov``: fraction of trials with the upper limit above the truth;
- ``med_err``: median reported error.

The "true" flux is the peak of the source at the nominal position: with a
pointing offset the forced fit is expected to be low (by 22% on average for
an offset of 0.3 beams in RA, on one arm only, and the attenuation of the
other arm).
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from astropy.table import Table

plt.rcParams.update({"font.size": 7, "axes.labelsize": 7, "legend.fontsize": 6})

KEYS = ["method", "noise", "drift", "rfi", "offset", "snr"]


def summarize(t):
    rows = []
    for g in t.group_by(KEYS).groups:
        true = g["snr"][0]
        flux = np.asarray(g["flux"])
        err = np.asarray(g["err"])
        resid = flux - true
        rows.append(
            dict(
                **{k: g[k][0] for k in KEYS},
                n=len(g),
                bias=np.mean(resid),
                scatter=np.std(flux),
                pull_std=np.std(resid / err),
                cov1=np.mean(np.abs(resid) < err),
                detected=np.mean(g["detected"]),
                ul_cov=np.mean(np.asarray(g["ul"]) >= true),
                med_err=np.median(err),
            )
        )
    return Table(rows=rows)


def rank(s):
    """Worst-case numbers of each method over the cells without offset.

    ``max_abs_bias_bright`` is at the highest S/N of the grid (30).
    """
    rows = []
    for g in s[s["offset"] == 0].group_by("method").groups:
        null = g[g["snr"] == 0]
        bright = g[g["snr"] == np.max(g["snr"])]
        rows.append(
            dict(
                method=g["method"][0],
                max_abs_pull_dev=np.max(np.abs(g["pull_std"] - 1)),
                min_cov1=np.min(g["cov1"]),
                max_false_det=np.max(null["detected"]),
                max_abs_bias_bright=np.max(np.abs(bright["bias"])),
                mean_med_err=np.mean(g["med_err"]),
            )
        )
    out = Table(rows=rows)
    out.sort("mean_med_err")
    return out


CASES = [
    ("white", "flat", False, "white, flat", "C0"),
    ("white", "drift", False, "white, drift", "C2"),
    ("1/f", "flat", False, "1/f, flat", "C1"),
    ("1/f", "drift", False, "1/f, drift", "C4"),
    ("1/f", "drift", True, "1/f, drift, RFI", "C3"),
]


def plot(s, best, fname):
    fig, axes = plt.subplots(1, 4, figsize=(7, 2.6), layout="constrained")
    for noise, drift, rfi, label, color in CASES:
        sel = (s["noise"] == noise) & (s["drift"] == drift) & (s["rfi"] == rfi)
        sel &= s["offset"] == 0
        for method, ls, mk in [("b_rfi", ":", "s"), (best, "-", "o")]:
            g = s[sel & (s["method"] == method)]
            g.sort("snr")
            kw = dict(color=color, ls=ls, marker=mk, ms=2, lw=0.8)
            lab = label if method == best else None
            axes[0].plot(g["snr"], g["bias"], label=lab, **kw)
            axes[1].plot(g["snr"], g["pull_std"], **kw)
            axes[2].plot(g["snr"], 100 * g["cov1"], **kw)
            axes[3].plot(g["snr"], g["scatter"], **kw)
    axes[0].set_ylabel("Bias [ideal errors]")
    axes[1].set_ylabel("Real scatter / reported error")
    axes[1].set_yscale("log")
    axes[1].axhline(1, color="grey", lw=0.5)
    axes[2].set_ylabel("Truth within 1-sigma error bar [%]")
    axes[2].axhline(68.3, color="grey", lw=0.5)
    axes[3].set_ylabel("Real scatter of the flux [ideal errors]")
    axes[0].axhline(0, color="grey", lw=0.5)
    for ax in axes:
        ax.set_xlabel("True S/N")
    handles, labels = axes[0].get_legend_handles_labels()
    handles += [
        plt.Line2D([], [], color="k", ls="-", marker="o", ms=2, lw=0.8),
        plt.Line2D([], [], color="k", ls=":", marker="s", ms=2, lw=0.8),
    ]
    labels += ["new method (d)", "script (b), RFI rejected"]
    fig.legend(handles, labels, loc="outside lower center", ncol=7, frameon=False)
    fig.savefig(fname, dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fname")
    parser.add_argument("--best", default="d_joint_win0.0_k3.0_mean")
    args = parser.parse_args()
    s = summarize(Table.read(args.fname))
    s.write("summary_crossscan.ecsv", overwrite=True)
    r = rank(s)
    r.write("ranking_crossscan.ecsv", overwrite=True)
    for c in r.colnames[1:]:
        r[c].format = ".3f"
    r.pprint_all()
    plot(s, args.best, "fig_crossscan.png")
    sel = s[(s["method"] == args.best) | (s["method"] == "b_rfi")]
    for c in ["bias", "scatter", "pull_std", "med_err"]:
        sel[c].format = ".2f"
    for c in ["cov1", "detected", "ul_cov"]:
        sel[c].format = ".3f"
    sel.sort(["rfi", "noise", "drift", "offset", "snr", "method"])
    sel.pprint_all()


if __name__ == "__main__":
    main()
