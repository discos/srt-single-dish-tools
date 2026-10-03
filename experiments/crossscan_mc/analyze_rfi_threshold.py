"""Summary and figure for ``run_rfi_threshold.py``.

Usage::

    python analyze_rfi_threshold.py mc_rfi_threshold.ecsv

Writes ``summary_rfi_threshold.ecsv`` and ``fig_rfi_threshold.png``. Fluxes in
units of the ideal error; ``nsigma = -1`` means no rejection.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from astropy.table import Table

plt.rcParams.update({"font.size": 7, "axes.labelsize": 7, "legend.fontsize": 6})

NCHSUB = 160  # subscans x channels per observation
CASES = [
    ("white", "flat", False, 0.0, "white, flat", "C0"),
    ("1/f", "drift", False, 0.0, "1/f, drift", "C1"),
    ("1/f", "drift", True, 0.0, "1/f, drift, RFI", "C3"),
    ("1/f", "drift", False, 0.3, "bright source, 0.3 beam offset", "k"),
]


def summarize(t):
    rows = []
    for g in t.group_by(["method", "noise", "drift", "rfi", "offset", "snr", "nsigma"]).groups:
        true = g["snr"][0]
        resid = np.asarray(g["flux"]) - true
        rows.append(
            dict(
                **{k: g[k][0] for k in ["method", "noise", "drift", "rfi", "offset", "snr"]},
                nsigma=g["nsigma"][0],
                n=len(g),
                bias=np.mean(resid),
                median_bias=np.median(resid),
                scatter=np.std(resid),
                err_ratio=np.std(resid) / np.median(g["err"]),
                detected=np.mean(g["detected"]),
                false_flag_frac=np.mean(g["nflag"] - (g["nrfi"] - g["nrfi_missed"])) / NCHSUB,
                missed_frac=(
                    np.sum(g["nrfi_missed"]) / np.sum(g["nrfi"]) if np.sum(g["nrfi"]) else np.nan
                ),
            )
        )
    return Table(rows=rows)


def plot(s, fname):
    fig, axes = plt.subplots(1, 4, figsize=(7, 2.8), layout="constrained")
    for noise, drift, rfi, offset, label, color in CASES:
        sel = (s["noise"] == noise) & (s["drift"] == drift) & (s["rfi"] == rfi)
        sel &= s["offset"] == offset
        for snr, ls in [(0, ":"), (30, "-"), (3000, "-")]:
            g = s[sel & (s["snr"] == snr) & (s["method"] == "b") & (s["nsigma"] > 0)]
            if len(g) == 0:
                continue
            g.sort("nsigma")
            kw = dict(color=color, ls=ls, marker="o", ms=2, lw=0.8)
            lab = f"{label}, S/N {snr}"
            axes[0].plot(g["nsigma"], 100 * g["false_flag_frac"], label=lab, **kw)
            if rfi:
                axes[1].plot(g["nsigma"], 100 * g["missed_frac"], **kw)
            if snr == 0:
                axes[2].plot(g["nsigma"], 100 * g["detected"], **kw)
            if snr == 30:
                axes[3].plot(g["nsigma"], g["bias"], **kw)
                gc = s[sel & (s["snr"] == snr) & (s["method"] == "c") & (s["nsigma"] > 0)]
                gc.sort("nsigma")
                axes[3].plot(gc["nsigma"], gc["bias"], **{**kw, "ls": "--"})
    axes[0].set_ylabel("Good subscans flagged [%]")
    axes[1].set_ylabel("RFI subscans not flagged [%]")
    axes[2].set_ylabel("(b) false detections at S/N 0 [%]")
    axes[3].set_ylabel("Bias at S/N 30 [ideal errors]\n(solid: (b), dashed: (c))")
    axes[3].axhline(0, color="grey", lw=0.5, ls=":")
    for ax in axes:
        ax.set_xlabel("rfi_nsigma")
    fig.legend(
        *axes[0].get_legend_handles_labels(), loc="outside lower center", ncol=4, frameon=False
    )
    fig.savefig(fname, dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fname")
    args = parser.parse_args()
    s = summarize(Table.read(args.fname))
    s.write("summary_rfi_threshold.ecsv", overwrite=True)
    plot(s, "fig_rfi_threshold.png")
    for c in ["bias", "median_bias", "scatter", "err_ratio"]:
        s[c].format = ".2f"
    for c in ["detected", "false_flag_frac", "missed_frac"]:
        s[c].format = ".4f"
    s.sort(["method", "rfi", "noise", "offset", "snr", "nsigma"])
    s.pprint_all()


if __name__ == "__main__":
    main()
