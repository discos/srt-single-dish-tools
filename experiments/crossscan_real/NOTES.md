# SDTcross (d) vs the cross-scan script (b) on the real SRT data

PR2, Stage B, step 4 (validation on real data), and the final demonstration
requested by Matteo: show, graphically and statistically, that the new
method works better on the real data.

## Files

| File | What it does |
|---|---|
| `validate_real.py` | Runs both methods on the six observations in `/Volumes/Elements/SRT` (read only): original data, channels, arms, time blocks, sign-flip null tests and fake-source injections. Writes `real_results.ecsv`. |
| `check_shape.py` | The joint-vs-window shape check (`CONS_SHP`) on the original and on 20 sign-flip null data sets per observation. Writes `shape_check.ecsv`. |
| `plot_real.py` | Statistics (`real_stats.ecsv`) and figures (`fig_real_profiles.png`, `fig_real_stats.png`). |

```bash
cd experiments/crossscan_real
PYTHONPATH=../..:../crossscan_mc ~/mamba/envs/py313/bin/python validate_real.py --npatterns 50
PYTHONPATH=../..:../crossscan_mc ~/mamba/envs/py313/bin/python plot_real.py real_results.ecsv
```

Environment: Mac, `py313` (scipy >= 1.15 is needed by method (b), which uses
`make_splrep`).

## Methods

- (b): the computing part of `plot_all_cross_scan_results.py`, verbatim
  (`experiments/crossscan_mc/methods.py`), on the `_cal` columns of the dump
  files as they are (flagged samples were stored as 0 by the SDTimage version
  that wrote them). Flux: TOTAL profile in the bin at CRPIX1 (labelled 0 in
  its plots); error: the error bar of that bin; detection: flux above its
  "3 sigma u.l." line.
- (d): `srttools.crossscan.measure_crosses` with the SDTcross defaults, on
  the same `_cal` columns, masking flagged samples (`-filt` columns). Beam
  FWHM from the band center (2.83 arcmin; the calibrators give 2.74-2.82).
  Detection: flux > 3 err.

## Tests

1. **Sign-flip null tests.** The source measured by (d) is subtracted (flux
   times the beam at the nominal position, good samples only), then each
   subscan (both channels) is multiplied by +1 or -1, with as many +1 as -1
   in each arm and direction of motion. The residual source cancels, the
   real noise (1/f, drifts, RFI, everything) stays. Over 50 random patterns
   per observation, the scatter of the flux is the real error, to compare
   with the reported one; and detections are false detections.
   - Without the subtraction, the reported errors of both methods are
     inflated by the source itself (each subscan has +-F): for SDSSJ133654
     (70 mJy) the error of (d) went from 1.6 to 7.9 mJy.
2. **Injection.** A fake point source of 5 times the error of (d) is added
   to the null data: detection rate and recovered flux.
3. **Time blocks.** The subscans are split in 8 consecutive blocks; the
   chi-square of the block fluxes around their weighted mean, per degree of
   freedom, should be ~1 (the source is assumed constant within an
   observation).
4. **LCP vs RCP, RA vs Dec**, with errors added in quadrature (for (d) the
   paired difference per subscan, `CONS_POL` in SDTcross, is more accurate,
   because the two channels share part of the noise).

Off-source positions along the arms are not used as null tests: they are at
most ~5 beams from the target, and a source outside the baseline window bends
the spline baseline (10-80% of its flux leaks into the measurement for a
source 2-4 beams away; see `srttools/tests/test_crossscan.py` and the SDTcross
tutorial).

## Results

`validate_real.py --npatterns 50` with the SDTcross defaults of commit
7ba805e (joint fit, knots every 3 FWHM, mean), `real_stats.ecsv`
(fluxes in mJy):

| obs | (b) flux | (d) flux | null scatter / error (b, d) | false det. (b, d) | blocks chi2/dof (b, d) | LCP-RCP [sigma] (b, d) | RA-Dec [sigma] (b, d) |
|---|---|---|---|---|---|---|---|
| HS1630+2355 0 | 6.58 +- 0.17 | 7.44 +- 0.33 | 2.68, 0.91 | 16%, 0% | 10.1, 0.32 | 3.9, 1.6 | 1.1, 0.3 |
| Mrk915 0 | 25.8 +- 1.5 | 26.1 +- 2.3 | 3.22, 0.97 | 30%, 0% | 5.8, 0.76 | 4.4, 2.0 | -2.6, -0.7 |
| SDSSJ094715 0 | 17.1 +- 0.6 | 15.9 +- 1.4 | 2.69, 1.12 | 14%, 0% | 12.0, 1.93 | 1.6, 0.6 | -1.1, -0.5 |
| SDSSJ094715 1 | 14.6 +- 0.6 | 15.4 +- 1.0 | 2.37, 1.04 | 10%, 0% | 9.0, 0.54 | 1.7, -0.8 | -1.4, -1.0 |
| SDSSJ133654 0 | 11.2 +- 9.5 | 68.9 +- 1.2 | 5.09, 1.02 | 46%, 0% | 9.0, 2.07 | 0.0, 1.2 | 5.1, 0.1 |
| SDSSJ133654 1 | 65.1 +- 0.6 | 71.8 +- 1.2 | 2.61, 0.98 | 14%, 0% | 12.5, 0.87 | -4.6, -2.8 | -3.2, -0.7 |

`fig_real_stats.png`: the four test statistics per observation.
`fig_real_profiles.png`: TOTAL profiles, (b) as plotted by the script (with
its "3 sigma u.l." line), (d) from SDTcross (baseline-subtracted, with the
fitted beam).

- **Errors.** On the sign-flip null data, the real scatter of (d) is
  0.91-1.12 times its reported error; for (b) 2.4-5.1 times. (b) finds a
  "source" in 10-46% of the null tests (its 3-sigma line is far too low),
  (d) in none of the 300.
- **Time blocks.** The fluxes of 8 consecutive blocks agree within the
  errors of (d) (chi2/dof 0.3-2.1, 7 degrees of freedom; the two values ~2
  have p ~ 0.05); with the errors of (b) they do not (chi2/dof 6-12).
- **Injection** of a source of 5 times the error of (d): detected in 94-100%
  of the cases by (d), bias < 0.3 error.
- **Channels and arms** agree within ~2 sigma for (d), except LCP-RCP for
  SDSSJ133654 Obs1 (-2.8 sigma, i.e. 4 mJy or 6% of 70 mJy: a relative
  calibration error between the channels would do it; the quadrature error
  used here also ignores the noise the channels share). For (b), 6 of 12
  comparisons exceed 2.5 sigma.
- **SDSSJ133654 Obs0**: (b) fails (11 +- 9.5 mJy, profile with a deep
  negative bowl and 130 mJy spikes). Removing the 7 subscans of about 1000
  flagged by the RFI rejection of SDTcross gives (b) 60.2 +- 0.6 mJy.
- (b) is lower than (d) for the bright SDSSJ133654 (60-65 vs 69-72 mJy),
  consistent with the -7% bias of (b) at high S/N in the Monte Carlo
  (`../crossscan_mc/NOTES.md`, section 7).

Caveat on the profiles: the (d) profiles are the data minus a baseline that
was fitted *together with* the beam template. Structure broader than the
beam is partly absorbed by the baseline, so these profiles look more like
the beam than the data justify. HS1630+2355 is visibly wider than the beam
in the (b) profile (~4' vs 2.8'), and the window fit (`joint=False`) gives
8.13 +- 0.46 mJy, against 7.44 +- 0.33 for the joint fit. For an extended
or mispointed source, both are lower limits to the total flux, and the
joint fit is the lower one (`../crossscan_mc/NOTES.md`, section 7).

## Shape check (joint vs window fit, `CONS_SHP`)

Added after the results above (decision with Matteo: keep the joint fit,
warn when the window fit disagrees). `check_shape.py --npatterns 20`:

| obs | original [sigma] | null: mean, std, max \|z\| |
|---|---|---|
| HS1630+2355 0 | -2.94 | -0.05, 0.72, 1.59 |
| Mrk915 0 | +0.22 | -0.11, 0.90, 1.83 |
| SDSSJ094715 0 | -0.88 | -0.30, 1.12, 2.60 |
| SDSSJ094715 1 | +0.50 | -0.07, 1.07, 2.46 |
| SDSSJ133654 0 | +1.95 | +0.17, 0.98, 2.38 |
| SDSSJ133654 1 | +3.87 | -0.16, 0.92, 2.04 |

- On null data the check is a standard normal (std 0.7-1.1), with no value
  beyond 3 in 120: no false alarms.
- HS1630+2355: -2.9 sigma (joint 7.44 +- 0.33, window 8.13 +- 0.46 mJy),
  just below the warning threshold: the source is probably wider than the
  beam or slightly mispointed, as its (b) profile suggests.
- SDSSJ133654 Obs1: +3.9 sigma (joint 71.8, window 68.7 mJy). Positive, so
  not a loss of the joint fit. The beam is narrower than the template (2.75'
  vs 2.83'), which costs the window fit 1.6% and the joint fit 0.4%
  (noiseless): with a 2.75' template the check goes to +2.6. The rest is
  probably the faint structure around the source, to which the window fit is
  more sensitive. The warning message depends on the sign.
