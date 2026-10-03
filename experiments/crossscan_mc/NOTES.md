# Monte Carlo characterisation of the current cross-scan methods

PR2, Stage A, step 2. Purpose: measure the bias, error coverage, false-detection
rate and upper-limit coverage of the three ways srttools users currently get a
flux (or upper limit) out of cross scans, before changing anything.

## Files

| File | What it does |
|---|---|
| `measure_real_noise.py` | Measures noise properties of the real data in `/Volumes/Elements/SRT` (read only). Writes `real_noise_properties.ecsv`, `real_noise_spectra.png`. |
| `mc_sim.py` | Simulates one observation as arrays, with the layout of the real data. |
| `methods.py` | The three current methods, (a), (b), (c), applied to simulated arrays. |
| `validate_e2e.py` | Checks that the arrays match what the real srttools pipeline produces from FITS files. |
| `run_mc.py` | Runs the grid (parallel). Writes `mc_results_current.ecsv`. |
| `analyze_mc.py` | Metrics (`summary_current.ecsv`) and figures (`fig_*_current.png`). |

Environment: Mac, `py313` micromamba env (numpy 2.4, scipy 1.17, astropy 7.2),
srttools from this worktree (`PYTHONPATH=../..:.`). scipy >= 1.15 is needed
because method (b) uses `make_splrep`.

```bash
cd experiments/crossscan_mc
PYTHONPATH=../..:. ~/mamba/envs/py313/bin/python measure_real_noise.py /Volumes/Elements/SRT/*_dump.hdf5
PYTHONPATH=../..:. ~/mamba/envs/py313/bin/python validate_e2e.py <scratch dir>
PYTHONPATH=../..:. ~/mamba/envs/py313/bin/python run_mc.py --ntrials 200 -o mc_results_current.ecsv
PYTHONPATH=../..:. ~/mamba/envs/py313/bin/python run_mc.py --rfi --ntrials 200 -o mc_results_rfi.ecsv
~/mamba/envs/py313/bin/python analyze_mc.py mc_results_current.ecsv mc_results_rfi.ecsv --tag current
```

Note: building a `ScanSet` without a config file makes srttools write a
`sample_config_file.ini` in the working directory. It is a by-product, deleted.

## 1. Real data (what the simulations reproduce)

From the six `*_dump.hdf5` files (`real_noise_properties.ecsv`):

- Sample time 0.0202 s; subscans of ~475-480 samples (9.6 s), spanning
  40.7 pixels = 10.3 beam FWHM (beam FWHM = 3.93 pixels, pixel = 0.0118 deg).
- 80 subscans per observation (40 per arm; HS1630 has 13 observations merged).
  Order: Dec ->, Dec <-, RA <-, RA ->, repeated.
- No samples are flagged (`<ch>-filt` is all True) in any of these dumps. RFI
  was handled in frequency: the bad spectral channels found with SDTrfistat
  were excluded with `--bad-intervals` before summing the band (see
  `/Volumes/Elements/README`). So the "flagged samples enter the averages as
  zeros" bug of `calculate_avg_cross` does not affect these data, but **the
  RFI that survives the channel cuts is still in the time series**.
- Residual RFI (subscans with detrended rms > 3x the median):
  - Mrk915 Obs0, LCP only: 11% of subscans (17% above 2x), in one
    contiguous stretch of ~16 subscans (~2.5 min); rms 3-8x normal,
    bursts lasting 1-3 s. RCP is clean. This explains the LCP/RCP
    disagreement found earlier for Mrk915.
  - SDSSJ133654 Obs0: one subscan with a ~24 Jy burst in both
    polarizations (120x the normal rms) over most of the subscan.
  - Others: below 1.3% of subscans.
- Pixel conventions: the nominal source sits at pixel coordinate 22.0 (on the
  edge between two histogram bins); `CRPIX1 = CRPIX2 = 22.5` (HS1630:
  `CRPIX1 = 23`, `CRPIX2 = 22.5`).
- Configuration: ALS baseline subtraction with default parameters and no
  excluded regions (the `.ini` files have none).
- Noise spectrum of the calibrated data (RFI subscans excluded, linear trend
  removed), fitted with sigma^2 (1 + (fknee/f)^alpha):
  - white level sigma = 11.5-14 mJy per sample (Mrk915 LCP: 20 mJy);
  - two regimes: HS1630, SDSSJ094715 (and Mrk915): fknee = 1.5-3 Hz,
    alpha = 1.2-1.7, LCP-RCP coherence below the knee 0.12-0.43;
    SDSSJ133654 (both days): fknee ~2 Hz, alpha = 2.2-2.9, coherence
    0.81-0.87;
  - LCP has a higher knee than RCP; above the knee the coherence is 0.06-0.2.
  - The source crosses the beam FWHM in ~0.9 s (~1 Hz), i.e. right where the
    1/f noise dominates.
- A first version of the spectra (`real_noise_spectra.png` before RFI
  rejection) made SDSSJ133654 Obs0 look anomalous: that was its single RFI
  subscan dominating the average.

## 2. Simulation set-up (`mc_sim.py`)

- Same layout as above: 80 subscans x 475 samples, alternating arms and
  directions, source at pixel 22.0, 43 pixels, `CRPIX` 22.5.
- Units: white noise rms per sample = 1, in each channel; two channels with
  equal gain.
- Noise models: `white`; `1/f` with fknee = 2 Hz, alpha = 1.5, LCP-RCP
  correlation of the 1/f part 0.4 (the first, more common regime above). 1/f noise is
  generated independently for each subscan.
- Baselines: `flat`; `drift`: per subscan, random slope (up to +-4 over the
  scan), random offset (+-5) and a random walk of amplitude 2, all common to
  both channels.
- Pointing offset: 0 or 0.3 beam FWHM in RA (RA-arm profile shifted by
  1.18 px; Dec-arm profile attenuated by a factor 0.78).
- The data are then baseline-subtracted with `srttools.fit.baseline_als`, as
  `Scan` does when loading.
- **S/N**: the source peak in units of the *ideal error*, the error of a
  forced-amplitude fit at the right position with white noise, no baseline,
  all subscans and both channels (`ideal_error()`: 0.0135 in units of the
  per-sample noise). With the real per-sample noise (~13 mJy), the ideal
  error is ~0.17 mJy. Same definition for all methods.
- Residual RFI (`rfi=True`), modelled on the two kinds seen in the real data:
  - an LCP-only episode: a contiguous block of 20% of the subscans, each
    with Poisson(3) Gaussian bursts of FWHM 1-3 s and peak 10-50 (white
    sigma);
  - one subscan with a huge burst in both channels (FWHM 5 s, peak ~1000,
    +-20% between channels).
  Check on one realization: 17/80 LCP subscans above 3x the median rms
  (Mrk915: 14/128), one subscan ~200x in both channels (SDSSJ133654 Obs0:
  120x).
- Grid: S/N in {0, 1, 2, 3, 5, 10, 30, 100} x noise {white, 1/f} x offset
  {0, 0.3} x baseline {flat, drift}, no RFI: 64 cells; plus the RFI cells
  (1/f noise, drift, no offset, 8 S/N values; added after Matteo pointed out
  the residual RFI). 200 trials each, one seed per trial (`[cell, trial]`,
  RFI cells numbered from 1000).

## 3. Methods and their outputs (`methods.py`)

None of the methods has a detection criterion; the one used here is what a
user would do with each output.

- **(a) SDTcal/SDTlcurve.** `fit_baseline_plus_bell` on each subscan, with the
  sky coordinate in degrees (RA ~ 248 deg), as `_treat_scan` does. Fits
  without a covariance matrix are dropped (`nfail`), as in `_treat_scan`;
  fits that raise an exception are counted in `ncrash` (in srttools this
  exception is *not* caught: SDTcal would stop). Flux = mean amplitude,
  error = RMS of the single-subscan errors (`calculate_src_flux`); the two
  channels are then averaged. No upper limit. Detected: flux > 3 x error.
- **(b) Matteo's script** (`plot_all_cross_scan_results.py`). The computing
  code is copied verbatim; inputs are the `_cal` columns after
  `calculate_avg_cross` (which de-spikes the data in place, as in SDTimage).
  Flux and error: TOTAL detrended profile in the bin labelled 0. Upper limit:
  the printed value, 3 x median of the TOTAL error outside the avoided region.
  Detected: flux > upper limit (the profile crosses the 3-sigma line).
- **(c) `ScanSet.calculate_avg_cross`**, the real method, on an in-memory
  `ScanSet`. Flux: mean of the four (arm x channel) profiles in the bin
  containing the source (pixel 22); error: the four `Sdev` combined in
  quadrature / 4. Upper limit: 3 x the median combined `Sdev` more than 6
  pixels from the source (the black "3 sigma u.l." line in Matteo's plots).
  Detected: flux > upper limit.

## 4. Validation against the real pipeline (`validate_e2e.py`)

One observation at S/N = 30 (white/flat, and 1/f/drift) written to FITS with
`save_scan` and loaded with srttools:

- `Scan` loading + ALS: channel data equal to the arrays used here to within
  3e-5 ideal errors (FITS float precision). Since (b) is verbatim and (c) is
  the real function, this validates their inputs.
- `CalibratorTable.from_scans` (SDTcal): amplitudes agree within 0.02 ideal
  errors; errors agree within 1% for all fits. In one of 150 fits a borderline
  fit fails in one pipeline and not in the other (float32 rounding).

## 5. Results

200 trials per cell; all numbers in units of the ideal error (S/N = true
flux). Figures: `fig_bias_current.png`, `fig_errratio_current.png`
(scatter of the measured flux / median reported error: >1 means the error
bars are too small), `fig_cov1_current.png`, `fig_detected_current.png`,
`fig_ulcov_current.png`, `fig_rfi_current.png`. Solid lines: no pointing
offset; dashed: 0.3 beam offset. Full table: `summary_current.ecsv`.

### (a) SDTcal/SDTlcurve Gaussian fit: unusable below S/N ~100

- Large positive bias at every S/N up to 30: median +40 (white noise, flat
  baseline) to +100 (1/f noise), mean up to +160. With pure noise
  (S/N = 0) it returns a "flux" of 40-100 ideal errors.
- Reported errors are absurd: 10^3-10^7 times the real scatter.
- 5-11 of the 160 single-subscan fits per observation are dropped (no
  covariance); in ~0.5% of the observations a fit raises an exception that
  SDTcal does not catch (the whole run would stop).
- Only at S/N = 100 with white noise and a flat baseline does it work
  (bias +3, detected), and even there errors are 9x too large.
- Likely causes, to be checked in Stage A3: each 9.6 s subscan is fitted on
  its own with free amplitude, position and width, so at low S/N the fit
  locks on to noise bumps (always positive); the width/amplitude degeneracy
  inflates the errors.

### (b) Cross-scan script: right flux, wrong errors and upper limits

- Flux: nearly unbiased at low S/N; at high S/N it is low by ~5% (-1.5 at
  S/N 30, -5.5 at S/N 100): the source sits on a bin edge and the bin value
  averages over a pixel on each side of the peak.
- Errors: 15-25% too small with white noise and a flat baseline (68% band
  contains the truth 58% of the time). With 1/f noise or a drifting
  baseline they are **2.5-3.3x too small**: the error is the scatter of the
  points in one bin, which does not see noise correlated on 1-10 s scales.
- False detections (flux above its 3-sigma line at S/N = 0): 1% with
  white/flat, **19-23%** with 1/f noise or a drifting baseline.
- The printed "upper limit" (3x median error away from the source) is a
  fixed number for a given data set, ~5-6, independent of the measured
  flux: below S/N 5 it is always above the true flux, from S/N 10 it is
  always below it. It is a 3-sigma *sensitivity* (itself underestimated
  ~3x with 1/f noise), not an upper limit.

### (c) `ScanSet.calculate_avg_cross`: decent flux, error bars 7-22x too large

- Real scatter of the flux is similar to (b) (1.6 with white/flat, 5-8
  with 1/f or drift), but the reported error is 37-55, so the 3-sigma
  "upper limit" is 110-160 ideal errors and **nothing is ever detected**,
  even at S/N = 100.
- Bias: +3 to +4 at S/N 0-3 with 1/f noise (+0.5 with white/flat),
  -10% at S/N 100 (bin dilution, as for (b), plus averaging four profiles).

### What the noise itself costs

With 1/f noise or a drifting baseline, the real scatter of both (b) and (c)
is 5-9 ideal errors, against 1.6-2.1 with white noise and a flat baseline.
This is the baseline (ALS) and 1/f noise leaking into the source bin, not a
bug of either method; a better method should reduce it, but the target is
not the ideal error.

### Pointing offset of 0.3 beams

All methods measure the flux at the nominal position, so an offset lowers
the flux: at S/N 30 by an extra ~4 ideal errors (~13%) for (b) and (c).
Neither method reports or corrects it.

### Residual RFI

`fig_rfi_current.png`: 1/f noise, drifting baseline, no offset, with
(dashed) and without (solid) the RFI model of section 2. The RFI dominates
everything:

- (a): bias ~+1300 at all S/N; errors still huge (no false detections only
  because of that).
- (b): bias ~+85 (median +56 to +74), even though the data were de-spiked
  by `calculate_avg_cross` before (b) uses them: the bursts last 1-3 s
  (~50-150 samples), too long for a spike filter. Errors grow to ~18 but are
  still 5x too small, and **71% of pure-noise observations are "detected"**.
- (c): bias ~+440; reported error ~1700, so it remains "consistent" and
  never detects.
- None of the methods notices that a minority of subscans is anomalous;
  rejecting outlier subscans (as `measure_real_noise.py` does, rms > 3x the
  median) or a robust average across subscans is needed.

### Summary for Stage A3

| | (a) Gaussian fit | (b) script | (c) calculate_avg_cross |
|---|---|---|---|
| Flux bias | +40 to +160 | ~0, -5% at high S/N | +3-4 at low S/N (1/f), -10% at high S/N |
| Errors (1/f or drift) | 10^3-10^7 x too large | 3x too small | 7-9x too large (22x white/flat) |
| False detections at S/N 0 | 0 | 19-23% | 0 |
| Detects S/N 10 | no | 65-98% | no |
| Upper limit | none | a sensitivity, not an UL | 110-160, always valid, useless |
| With residual RFI | bias +1300 | bias +85, 71% false detections | bias +440 |
