"""Measure the noise properties of real cross-scan data, to set up the simulations.

For each observation (``*_dump.hdf5`` written by SDTimage), and each subscan
(``Scan_id``) whose samples are all good (``<ch>-filt``):

- sample time, samples per subscan, number of subscans, scan length in beams;
- power spectrum of the calibrated data (``<ch>_cal``), after removing a
  straight line (part of any baseline model) and applying a Hann window;
  fitted with ``sigma**2 * (1 + (fknee / f)**alpha)``;
- coherence between LCP and RCP as a function of frequency. At frequencies
  well below the knee, it is the correlation of the 1/f part.

Subscans hit by residual RFI (rms of the detrended data above 3 times the
median, in either channel) are excluded from the spectra, otherwise a single
burst can dominate the average (SDSSJ133654 Obs0). Their fraction is
reported (``rfi_frac_<ch>``).

Input files are only read. Usage::

    python measure_real_noise.py /Volumes/Elements/SRT/*_dump.hdf5
"""

import os
import sys

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import detrend

from astropy.io import fits
from astropy.table import Table

CHANS = ["Feed0_LCP", "Feed0_RCP"]


def red_psd(f, logsigma2, fknee, alpha):
    return np.log(np.exp(logsigma2) * (1 + (np.abs(fknee) / f) ** alpha))


def subscan_spectra(t):
    """Windowed FFTs of the good subscans of both channels, grouped by length.

    Also returns the fraction of subscans rejected as RFI-contaminated, per
    channel.
    """
    subs = []
    for sid in np.unique(t["Scan_id"]):
        s = t[t["Scan_id"] == sid]
        if len(s) < 64:
            continue
        if not all(np.all(s[f"{ch}-filt"]) for ch in CHANS):
            continue
        subs.append((s, [detrend(np.asarray(s[f"{ch}_cal"])) for ch in CHANS]))
    rms = np.array([[np.std(y) for y in ys] for _, ys in subs])
    rfi = rms > 3 * np.median(rms, axis=0)
    groups = {}
    dts = []
    for (s, ys), bad in zip(subs, rfi):
        if np.any(bad):
            continue
        dts.append(np.median(np.diff(s["time"])) * 86400)
        n = len(s)
        win = np.hanning(n)
        groups.setdefault(n, []).append([np.fft.rfft(y * win) for y in ys])
    return groups, np.median(dts), rfi.mean(axis=0)


def analyze(fname):
    t = Table.read(fname)
    root = fname.replace("_dump.hdf5", "")
    hd = fits.getheader(root + "_cal.fits", 1)
    beam_fwhm_pix = hd["BMAJ"] / abs(hd["CDELT1"])
    nsub = len(np.unique(t["Scan_id"]))
    sizes = np.array([np.sum(t["Scan_id"] == s) for s in np.unique(t["Scan_id"])])
    span = []
    for sid in np.unique(t["Scan_id"])[:20]:
        s = t[t["Scan_id"] == sid]
        coord = s["x"][:, 0] if s["direction"][0] else s["y"][:, 0]
        span.append(np.ptp(coord))
    groups, dt, rfi_frac = subscan_spectra(t)

    # Use the most common subscan length, so that frequencies line up
    n = max(groups, key=lambda k: len(groups[k]))
    ffts = np.array(groups[n])  # (nsub, nchan, nfreq)
    freqs = np.fft.rfftfreq(n, dt)
    norm = dt / np.sum(np.hanning(n) ** 2)  # one-sided-ish PSD normalisation
    p = np.mean(np.abs(ffts) ** 2, axis=0) * norm
    cross = np.mean(ffts[:, 0] * np.conj(ffts[:, 1]), axis=0) * norm
    coh = np.abs(cross) / np.sqrt(p[0] * p[1])

    res = dict(
        name=os.path.basename(root),
        dt=dt,
        nsub=nsub,
        nsub_used=len(ffts),  # subscans with the most common length
        nsamp_median=int(np.median(sizes)),
        span_pix=np.median(span),
        beam_fwhm_pix=beam_fwhm_pix,
        **{f"rfi_frac_{ch}": f for ch, f in zip(CHANS, rfi_frac)},
    )
    fit_mask = freqs >= 2 / (n * dt)
    for i, ch in enumerate(CHANS):
        p0 = [np.log(np.median(p[i][fit_mask][-20:])), 1.0, 1.0]
        try:
            par, _ = curve_fit(
                red_psd, freqs[fit_mask], np.log(p[i][fit_mask]), p0=p0, maxfev=20000
            )
            # white sigma per sample: PSD level = sigma^2 * 2 * dt (one-sided)
            res[f"{ch}_sigma_mJy"] = np.sqrt(np.exp(par[0]) / (2 * dt)) * 1e3
            res[f"{ch}_fknee"] = abs(par[1])
            res[f"{ch}_alpha"] = par[2]
        except RuntimeError:
            res[f"{ch}_fknee"] = np.nan
        res[f"{ch}_rms_mJy"] = np.std(np.concatenate([np.fft.irfft(f[i], n) for f in ffts])) * 1e3
    low = fit_mask & (freqs < 0.5 * np.nanmin([res[f"{c}_fknee"] for c in CHANS]))
    res["coherence_low"] = np.mean(coh[low]) if np.any(low) else np.nan
    res["coherence_high"] = np.mean(coh[freqs > freqs[-1] / 2])
    return res, freqs, p, coh


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    files = sys.argv[1:]
    outdir = os.path.dirname(os.path.abspath(__file__))
    fig, axes = plt.subplots(1, 2, figsize=(7, 3), layout="constrained")
    rows = []
    for fname in sorted(files):
        res, freqs, p, coh = analyze(fname)
        rows.append(res)
        print(res)
        axes[0].loglog(freqs[1:], p[0][1:], lw=0.7, label=res["name"].split("_")[2])
        axes[1].semilogx(freqs[1:], coh[1:], lw=0.7)
    axes[0].set_xlabel("Frequency (Hz)", fontsize=7)
    axes[0].set_ylabel("PSD LCP (Jy$^2$/Hz)", fontsize=7)
    axes[1].set_xlabel("Frequency (Hz)", fontsize=7)
    axes[1].set_ylabel("LCP-RCP coherence", fontsize=7)
    axes[0].legend(fontsize=6)
    for ax in axes:
        ax.tick_params(labelsize=7)
    fig.savefig(os.path.join(outdir, "real_noise_spectra.png"), dpi=200)
    Table(rows=rows).write(os.path.join(outdir, "real_noise_properties.ecsv"), overwrite=True)
