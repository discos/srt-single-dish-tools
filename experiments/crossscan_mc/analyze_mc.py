"""Metrics and figures from the Monte Carlo results of ``run_mc.py``.

Usage::

    python analyze_mc.py mc_results_current.ecsv mc_results_rfi.ecsv [--tag current]

All fluxes are in units of the ideal error (``mc_sim.ideal_error``); the true
flux is the source peak, ``snr``. Writes ``fig_<metric>_<tag>.png`` and
``summary_<tag>.ecsv`` (one row per method and grid cell).
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from astropy.table import Table, vstack

plt.rcParams.update({"font.size": 7, "axes.labelsize": 7, "legend.fontsize": 6})

COLORS = {"a": "C0", "b": "C1", "c": "C2"}
LABELS = {
    "a": "(a) SDTlcurve Gaussian fit",
    "b": "(b) cross-scan script",
    "c": "(c) calculate_avg_cross",
}
PANELS = [("white", "flat"), ("white", "drift"), ("1/f", "flat"), ("1/f", "drift")]


def summarize(t):
    rows = []
    keys = ["method", "noise", "drift", "offset", "rfi", "snr"]
    for g in t.group_by(keys).groups:
        true = g["snr"][0]
        flux, err, ul = np.array(g["flux"]), np.array(g["err"]), np.array(g["ul"])
        resid = flux - true
        rows.append(
            dict(
                **{k: g[k][0] for k in keys},
                n=len(g),
                bias=np.mean(resid),
                bias_err=np.std(resid) / np.sqrt(len(g)),
                median_bias=np.median(resid),
                scatter=np.std(flux),
                median_err=np.median(err),
                # >1: errors too small; <1: too large
                err_ratio=np.std(flux) / np.median(err),
                cov1=np.mean(np.abs(resid) < err),
                cov2=np.mean(np.abs(resid) < 2 * err),
                detected=np.mean(g["detected"]),
                ul_cov=np.mean(true <= ul) if np.all(np.isfinite(ul)) else np.nan,
                median_ul=np.median(ul),
                nfail=np.mean(g["nfail"]),
                ncrash=np.mean(g["ncrash"]),
            )
        )
    return Table(rows=rows)


def panel_plot(s, column, ylabel, fname, methods="abc", yscale="linear", hlines=(), ylim=None):
    s = s[~s["rfi"]]
    fig, axes = plt.subplots(1, 4, figsize=(7, 2.3), sharey=True, layout="constrained")
    for ax, (noise, drift) in zip(axes, PANELS):
        for m in methods:
            for offset, ls in [(0.0, "-"), (0.3, "--")]:
                g = s[
                    (s["method"] == m)
                    & (s["noise"] == noise)
                    & (s["drift"] == drift)
                    & (s["offset"] == offset)
                ]
                g.sort("snr")
                x = np.array(g["snr"], dtype=float)
                x[x == 0] = 0.5  # show S/N=0 on the log axis
                ax.plot(
                    x,
                    g[column],
                    ls=ls,
                    marker="o",
                    ms=2,
                    lw=0.8,
                    color=COLORS[m],
                    label=LABELS[m] if offset == 0 else None,
                )
        for h in hlines:
            ax.axhline(h, color="grey", lw=0.5, ls=":")
        ax.set_xscale("log")
        if yscale == "log" and not np.any(np.asarray(s[column]) > 0):
            yscale = "linear"
        ax.set_yscale(yscale)
        ax.set_title(f"{noise} noise, {drift} baseline", fontsize=7)
        ax.set_xlabel("True S/N (S/N=0 plotted at 0.5)")
        if ylim is not None:
            ax.set_ylim(ylim)
    axes[0].set_ylabel(ylabel)
    axes[0].legend(loc="best", frameon=False)
    axes[-1].plot([], [], "k--", lw=0.8, label="0.3 beam offset")
    axes[-1].legend(loc="best", frameon=False)
    fig.savefig(fname, dpi=300)
    plt.close(fig)


RFI_METRICS = [
    ("bias", "Mean (measured - true)\n[ideal errors]", "symlog", [0]),
    ("err_ratio", "Scatter / reported error", "log", [1]),
    ("detected", "Fraction detected", "linear", [0.00135]),
    ("ul_cov", "Fraction with true flux $\\leq$ UL", "linear", [0.99865]),
]


def rfi_plot(s, fname):
    """1/f noise, drifting baseline, no offset: without and with RFI."""
    s = s[(s["noise"] == "1/f") & (s["drift"] == "drift") & (s["offset"] == 0)]
    fig, axes = plt.subplots(1, 4, figsize=(7, 2.3), layout="constrained")
    for ax, (column, ylabel, yscale, hlines) in zip(axes, RFI_METRICS):
        for m in "abc":
            for rfi, ls in [(False, "-"), (True, "--")]:
                g = s[(s["method"] == m) & (s["rfi"] == rfi)]
                g.sort("snr")
                x = np.array(g["snr"], dtype=float)
                x[x == 0] = 0.5
                ax.plot(
                    x,
                    g[column],
                    ls=ls,
                    marker="o",
                    ms=2,
                    lw=0.8,
                    color=COLORS[m],
                    label=LABELS[m] if not rfi else None,
                )
        for h in hlines:
            ax.axhline(h, color="grey", lw=0.5, ls=":")
        ax.set_xscale("log")
        ax.set_yscale(yscale)
        ax.set_ylabel(ylabel)
        ax.set_xlabel("True S/N")
    axes[0].legend(loc="best", frameon=False)
    axes[1].plot([], [], "k--", lw=0.8, label="with RFI")
    axes[1].plot([], [], "k-", lw=0.8, label="no RFI")
    axes[1].legend(loc="best", frameon=False)
    fig.suptitle(
        "1/f noise, drifting baseline, no pointing offset (S/N=0 plotted at 0.5)", fontsize=7
    )
    fig.savefig(fname, dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fnames", nargs="+")
    parser.add_argument("--tag", default="current")
    args = parser.parse_args()
    tables = []
    for fname in args.fnames:
        t = Table.read(fname)
        if "rfi" not in t.colnames:
            t["rfi"] = False
        tables.append(t)
    t = vstack(tables)
    s = summarize(t)
    s.write(f"summary_{args.tag}.ecsv", overwrite=True)
    tag = args.tag
    panel_plot(
        s,
        "bias",
        "Mean (measured - true)\n[ideal errors]",
        f"fig_bias_{tag}.png",
        yscale="symlog",
        hlines=[0],
    )
    panel_plot(
        s,
        "err_ratio",
        "Scatter / reported error",
        f"fig_errratio_{tag}.png",
        yscale="log",
        hlines=[1],
    )
    panel_plot(
        s, "cov1", "Fraction within 1 reported $\\sigma$", f"fig_cov1_{tag}.png", hlines=[0.683]
    )
    panel_plot(
        s,
        "detected",
        "Fraction detected",
        f"fig_detected_{tag}.png",
        hlines=[0.00135],
        ylim=(-0.02, 1.02),
    )
    panel_plot(
        s,
        "ul_cov",
        "Fraction with true flux $\\leq$ UL",
        f"fig_ulcov_{tag}.png",
        methods="bc",
        hlines=[0.99865],
        ylim=(-0.02, 1.02),
    )
    if np.any(s["rfi"]):
        rfi_plot(s, f"fig_rfi_{tag}.png")
    with np.printoptions(precision=3):
        sel = s[(s["offset"] == 0) & np.isin(s["snr"], [0, 5, 30])]
        sel.sort(["method", "rfi", "noise", "drift", "snr"])
        for c in ["bias", "median_bias", "scatter", "median_err", "err_ratio", "median_ul"]:
            sel[c].format = ".2f"
        for c in ["cov1", "cov2", "detected", "ul_cov", "nfail", "ncrash"]:
            sel[c].format = ".3f"
        sel.remove_columns(["bias_err", "n", "offset"])
        sel.pprint_all()


if __name__ == "__main__":
    main()
