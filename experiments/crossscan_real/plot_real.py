"""Statistics and figures from ``validate_real.py``.

Usage (from this directory)::

    PYTHONPATH=../..:../crossscan_mc python plot_real.py real_results.ecsv

Writes ``real_stats.ecsv`` (one row per observation and method),
``fig_real_profiles.png`` (TOTAL profiles of both methods, recomputed from
the data) and ``fig_real_stats.png``.

Statistics, per observation and method (all should be close to the value in
brackets if the errors are right):

- ``null_ratio``: scatter of the null-test fluxes over the sign patterns,
  divided by the median error reported on the null data [1];
- ``null_false``: fraction of null tests with a detection [~0.001];
- ``inject_det``: fraction of injected sources (5 times the error of (d))
  detected;
- ``inject_bias``: mean recovered minus injected flux, in units of the
  reported error [0];
- ``block_chi2``: chi-square per degree of freedom of the fluxes of the time
  blocks around their weighted mean [1];
- ``z_pol``, ``z_arm``: LCP - RCP and RA - Dec, in units of the error of the
  difference (errors added in quadrature) [|z| < 2 most of the time].
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
import validate_real as v

from astropy.table import Table

plt.rcParams.update({"font.size": 7, "axes.labelsize": 7, "legend.fontsize": 6})

SHORT = {
    "HS1630+2355_Obs0": "HS1630 0",
    "Mrk915_Obs0": "Mrk915 0",
    "SDSSJ094715_Obs0": "J0947 0",
    "SDSSJ094715_Obs1": "J0947 1",
    "SDSSJ133654_Obs0": "J1336 0",
    "SDSSJ133654_Obs1": "J1336 1",
}


def statistics(t):
    rows = []
    for obs in v.OBSERVATIONS:
        for method in ["b", "d"]:
            sel = (t["obs"] == obs) & (t["method"] == method)

            def get(test):
                return t[sel & (t["test"] == test)]

            null, inject = get("null"), get("inject")
            blocks = t[sel & np.char.startswith(np.asarray(t["test"], dtype=str), "block")]
            c2, dof = v.chi2(blocks["flux"], blocks["err"])

            def z(a, b):
                a, b = get(a)[0], get(b)[0]
                return (a["flux"] - b["flux"]) / np.hypot(a["err"], b["err"])

            orig = get("original")[0]
            rows.append(
                dict(
                    obs=obs,
                    method=method,
                    flux=orig["flux"],
                    err=orig["err"],
                    ul=orig["ul"],
                    detected=orig["detected"],
                    null_ratio=np.std(null["flux"]) / np.median(null["err"]),
                    null_false=np.mean(null["detected"]),
                    inject_det=np.mean(inject["detected"]),
                    inject_bias=np.mean((inject["flux"] - inject["true"]) / inject["err"]),
                    block_chi2=c2 / dof,
                    z_pol=z(*v.CHANS),
                    z_arm=z("RA", "Dec"),
                )
            )
    return Table(rows=rows)


def plot_stats(s, fname):
    fig, axes = plt.subplots(1, 4, figsize=(7, 2.4), layout="constrained")
    x = np.arange(len(v.OBSERVATIONS))
    for method, color, marker, label in [
        ("b", "C0", "s", "script (b)"),
        ("d", "C3", "o", "SDTcross (d)"),
    ]:
        g = s[s["method"] == method]
        kw = dict(color=color, marker=marker, ls="", ms=4, label=label)
        axes[0].plot(x, g["null_ratio"], **kw)
        axes[1].plot(x, g["block_chi2"], **kw)
        axes[2].plot(x - 0.1, np.abs(g["z_pol"]), **kw)
        axes[2].plot(x + 0.1, np.abs(g["z_arm"]), **{**kw, "mfc": "none", "label": None})
        axes[3].plot(x, 100 * g["null_false"], **kw)
    axes[0].set_ylabel("Null scatter / reported error")
    axes[1].set_ylabel(r"Time blocks: $\chi^2$/d.o.f.")
    axes[2].set_ylabel("|LCP-RCP| (filled), |RA-Dec| (empty) [sigma]")
    axes[3].set_ylabel("Null tests with a detection [%]")
    for ax, ref in zip(axes, [1, 1, 2, 0]):
        ax.axhline(ref, color="grey", lw=0.5)
        ax.set_xticks(x, [SHORT[o] for o in v.OBSERVATIONS], rotation=60)
    for ax in axes[:2]:
        ax.set_yscale("log")
    axes[0].legend(frameon=False)
    fig.savefig(fname, dpi=300)
    plt.close(fig)


def plot_profiles(fname):
    """TOTAL profiles of (b) and (d) for each observation."""
    from methods import script_b

    from astropy.io import fits

    nobs = len(v.OBSERVATIONS)
    fig, axes = plt.subplots(nobs, 2, figsize=(7, 1.45 * nobs), layout="constrained")
    for i, obs in enumerate(v.OBSERVATIONS):
        table, crpix, ra0, dec0, beam = v.load(obs)
        root = f"{v.DATADIR}/CCB_SARDARA_{obs}"
        pix_arcmin = abs(fits.getheader(root + "_cal.fits", 1)["CDELT2"]) * 60
        xb, prof, err, ul, _ = script_b(table, v.extnames(), crpix)
        ax = axes[i, 0]
        ax.errorbar(xb * pix_arcmin, prof * 1e3, err * 1e3, fmt="o", ms=1.5, lw=0.5, color="k")
        ax.axhline(ul * 1e3, color="C0", lw=0.8, ls="--")
        ax.axhline(0, color="grey", lw=0.5)
        i0 = np.argmin(np.abs(xb))
        ax.set_title(
            f"{obs} (b): {prof[i0] * 1e3:.2f} ± {err[i0] * 1e3:.2f} mJy, "
            f"3σ line {ul * 1e3:.2f}",
            fontsize=6,
        )

        d = v.run_d(table, ra0, dec0, beam)
        centers, dprof, derr = d["result"]["profiles"]["TOTAL"]
        sigma = beam / (2 * np.sqrt(2 * np.log(2)))
        ax = axes[i, 1]
        ax.errorbar(centers * 60, dprof * 1e3, derr * 1e3, fmt="o", ms=1.5, lw=0.5, color="k")
        model = d["flux"] * np.exp(-(centers**2) / (2 * sigma**2))
        ax.plot(centers * 60, model * 1e3, color="C3", lw=0.8)
        ax.axhline(0, color="grey", lw=0.5)
        status = "detected" if d["detected"] else f"UL {d['ul'] * 1e3:.2f}"
        ax.set_title(
            f"{obs} (d): {d['flux'] * 1e3:.2f} ± {d['err'] * 1e3:.2f} mJy ({status})",
            fontsize=6,
        )
        for ax in axes[i]:
            ax.set_xlim(-16, 16)
            ax.tick_params(labelsize=6)
        axes[i, 0].set_ylabel("mJy")
    for ax in axes[-1]:
        ax.set_xlabel("Offset along the arm [arcmin]")
    fig.savefig(fname, dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fname")
    args = parser.parse_args()
    v.quiet()
    s = statistics(Table.read(args.fname))
    s.write("real_stats.ecsv", overwrite=True)
    plot_stats(s, "fig_real_stats.png")
    plot_profiles("fig_real_profiles.png")
    for c in ["flux", "err", "ul"]:
        s[c] = s[c] * 1e3
        s[c].format = ".2f"
    for c in [
        "null_ratio",
        "null_false",
        "inject_det",
        "inject_bias",
        "block_chi2",
        "z_pol",
        "z_arm",
    ]:
        s[c].format = ".2f"
    s.pprint_all()


if __name__ == "__main__":
    main()
