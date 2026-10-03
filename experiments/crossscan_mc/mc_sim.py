"""Simulate one cross-scan observation as arrays, with the layout of the real data.

The numbers below come from ``measure_real_noise.py`` run on the SRT data in
``/Volumes/Elements/SRT`` (see ``NOTES.md``). Everything is in "pixel" units
along each arm (one pixel = ``PIX_DEG`` degrees) and in units of the white
noise rms per sample (``SIGMA_WHITE = 1``).

The output of ``simulate_observation`` is what srttools holds after loading the
scans: each subscan is already baseline-subtracted with
``srttools.fit.baseline_als``, as ``Scan.baseline_subtract`` does by default
(no excluded regions, as in the real configuration files).
"""

import numpy as np

from srttools.fit import baseline_als
from srttools.simulate import _create_baseline, _red_noise

# --- Layout of the real observations ---
DT = 0.0202  # s, sample time
NSAMP = 475  # samples per subscan
SPAN_PIX = 40.7  # scan length, pixels
NSUBSCANS = 80  # 40 per arm
PIX_DEG = 0.011778  # pixel size, deg
BEAM_FWHM_PIX = 3.933
BEAM_SIGMA_PIX = BEAM_FWHM_PIX / (2 * np.sqrt(2 * np.log(2)))
SRC_PIX = 22.0  # pixel coordinate of the nominal source position
CRPIX = 22.5  # CRPIX1 and CRPIX2 in the _cal.fits files
NPIX = 43
RA0, DEC0 = 248.26, 23.83  # deg, like HS1630+2355
CHANS = ["Feed0_LCP", "Feed0_RCP"]
# Arm and direction of motion, repeated: Dec ->, Dec <-, RA <-, RA ->
ARM_PATTERN = [("Dec", 1), ("Dec", -1), ("RA", -1), ("RA", 1)]

# --- Noise models ---
SIGMA_WHITE = 1.0
NOISE_MODELS = {
    "white": dict(fknee=None),
    # Typical values measured on the real data (NOTES.md)
    "1/f": dict(fknee=2.0, noise_index=1.5, pol_correlation=0.4),
}
# Baseline drift, common to both channels: random slope up to +-2 sigma over
# the scan, random offset, plus a random walk of rms ~2 sigma.
DRIFT_MODELS = {
    "flat": "flat",
    "drift": ((-4 / SPAN_PIX, 4 / SPAN_PIX), (-5, 5), 2.0),
}

# --- Residual RFI, modelled on the real data (NOTES.md) ---
# 1. An episode in LCP only: a contiguous block of 20% of the subscans, each
#    with Poisson(3) bursts of 1-3 s (FWHM), peak 10-50 sigma, positive.
# 2. One subscan with a huge burst in both channels: FWHM ~5 s, peak ~1000
#    sigma (+-20% between channels).
RFI_EPISODE_FRAC = 0.2
RFI_BURSTS_PER_SUBSCAN = 3
RFI_BURST_FWHM_S = (1.0, 3.0)
RFI_BURST_PEAK = (10.0, 50.0)
RFI_HUGE_PEAK = 1000.0
RFI_HUGE_FWHM_S = 5.0


def _burst(time, center, fwhm, peak):
    return peak * np.exp(-((time - center) ** 2) / (2 * (fwhm / 2.3548) ** 2))


def rfi_contamination(rng):
    """RFI to add to each subscan and channel: dict {(isub, ich): array}."""
    rfi = {}
    time = np.arange(NSAMP) * DT
    nblock = int(round(RFI_EPISODE_FRAC * NSUBSCANS))
    start = rng.integers(0, NSUBSCANS - nblock + 1)
    for isub in range(start, start + nblock):
        signal = np.zeros(NSAMP)
        for _ in range(rng.poisson(RFI_BURSTS_PER_SUBSCAN)):
            signal += _burst(
                time,
                rng.uniform(0, time[-1]),
                rng.uniform(*RFI_BURST_FWHM_S),
                rng.uniform(*RFI_BURST_PEAK),
            )
        rfi[(isub, 0)] = signal
    isub = rng.integers(0, NSUBSCANS)
    huge = _burst(time, rng.uniform(0.3, 0.7) * time[-1], RFI_HUGE_FWHM_S, RFI_HUGE_PEAK)
    for ich in range(len(CHANS)):
        rfi[(isub, ich)] = rfi.get((isub, ich), 0) + huge * rng.uniform(0.8, 1.2)
    return rfi


def subscan_layout():
    """Arm, direction and sample positions (pixels) of each subscan."""
    x = SRC_PIX + np.linspace(-SPAN_PIX / 2, SPAN_PIX / 2, NSAMP)
    layout = []
    for i in range(NSUBSCANS):
        arm, sign = ARM_PATTERN[i % len(ARM_PATTERN)]
        layout.append((arm, x if sign > 0 else x[::-1]))
    return layout


def beam_profile(x, center):
    return np.exp(-((x - center) ** 2) / (2 * BEAM_SIGMA_PIX**2))


def ideal_error():
    """Error of the forced-fit amplitude with white noise, no baseline, all data.

    Each sample contributes ``g_i**2 / sigma**2`` to the inverse variance,
    where ``g_i`` is the beam profile at the sample position; summed over all
    subscans and both channels.
    """
    info = sum(np.sum(beam_profile(x, SRC_PIX) ** 2) for _, x in subscan_layout())
    return SIGMA_WHITE / np.sqrt(len(CHANS) * info)


def simulate_observation(
    amplitude, noise="white", offset_beams=0.0, drift="flat", rfi=False, rng=None
):
    """Simulate one observation.

    Parameters
    ----------
    amplitude : float
        Peak of the source, in units of the white-noise rms per sample
    noise : str
        Key of ``NOISE_MODELS``
    offset_beams : float
        Pointing offset of the true source in RA, in beam FWHMs: the RA arm
        profile is shifted, the Dec arm profile is attenuated
    drift : str
        Key of ``DRIFT_MODELS``
    rfi : bool
        Add residual RFI (``rfi_contamination``)
    rng : `numpy.random.Generator`

    Returns
    -------
    subscans : list of dict
        Keys: ``arm`` ("RA" or "Dec"), ``x`` (pixel coordinate along the arm),
        ``time`` (s), and one array per channel, baseline-subtracted;
        ``<channel>_has_rfi`` tells whether RFI was added
    """
    rng = np.random.default_rng(rng)
    nmodel = NOISE_MODELS[noise]
    dx = offset_beams * BEAM_FWHM_PIX
    rfi_signal = rfi_contamination(rng) if rfi else {}
    subscans = []
    t0 = 0.0
    for isub, (arm, x) in enumerate(subscan_layout()):
        if arm == "RA":
            signal = amplitude * beam_profile(x, SRC_PIX + dx)
        else:
            signal = (
                amplitude * beam_profile(x, SRC_PIX) * np.exp(-(dx**2) / (2 * BEAM_SIGMA_PIX**2))
            )
        base = _create_baseline(x - SRC_PIX, DRIFT_MODELS[drift], rng=rng)
        if nmodel["fknee"] is not None:
            red = _red_noise(
                NSAMP,
                DT,
                nmodel["fknee"],
                SIGMA_WHITE,
                noise_index=nmodel["noise_index"],
                nchan=len(CHANS),
                correlation=nmodel["pol_correlation"],
                rng=rng,
            )
        else:
            red = np.zeros((len(CHANS), NSAMP))
        time = t0 + np.arange(NSAMP) * DT
        t0 = time[-1] + 1.0
        sub = dict(arm=arm, x=x, time=time)
        for ich, ch in enumerate(CHANS):
            raw = signal + base + red[ich] + rng.normal(0, SIGMA_WHITE, NSAMP)
            raw = raw + rfi_signal.get((isub, ich), 0)
            sub[ch] = baseline_als(time, raw)
            sub[ch + "_raw"] = raw
            sub[ch + "_has_rfi"] = (isub, ich) in rfi_signal
        subscans.append(sub)
    return subscans
