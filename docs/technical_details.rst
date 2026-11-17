Technical details
=================

This page collects implementation and maintenance notes that are useful to
developers, and the reasoning behind some non-obvious choices.

Supported Python and dependency versions
----------------------------------------

The package supports Python 3.12 and newer. Python 3.8 was the minimum until
2026, for old production machines at the telescope; it was raised to 3.12
because Python 3.8 to 3.10 have reached their end of life, and current releases of
the scientific stack (numpy, scipy, astropy, matplotlib) all support 3.12.

The lower bounds of the dependencies in ``pyproject.toml`` are not guesses:
they are the oldest versions that install from binary wheels on Python 3.12
and pass the test suite. They are exercised by the ``py312-test-oldestdeps``
tox environment, which uses ``uv`` with ``--resolution lowest-direct`` to
install exactly the lower bound of every direct dependency.

Most bounds (``numpy>=1.26``, ``scipy>=1.11.2``, ``astropy>=5.3.4``,
``matplotlib>=3.7.3``, ``h5py>=3.10``, ``pyyaml>=6.0.1``) are simply the first
releases with wheels for Python 3.12, rather than requirements of the code.
Matplotlib 3.8.0 is excluded because of a serious bug in that release.

Running the tests
-----------------

The recommended way is through tox, with the ``tox-uv`` plugin:

.. code-block:: console

    $ pip install tox tox-uv
    $ tox -e py313-test-alldeps     # latest versions, all optional dependencies
    $ tox -e py312-test-oldestdeps  # oldest supported versions
    $ tox -e py313-test-devdeps     # nightly builds of numpy, scipy, matplotlib, astropy
    $ tox -e codestyle              # pre-commit checks (ruff, ...)

Tox installs the package in a separate environment and runs the tests from
``.tmp/<envname>``, so that the source tree is not imported by mistake. On
macOS, if tox picks up a "universal2" Python from python.org, the oldest
versions of numpy will not install; ``UV_PYTHON_PREFERENCE=only-managed``
forces tox-uv to use a uv-managed Python instead.

The code style checks deliberately do not include ``codespell``: with
``--write-changes`` it silently "corrected" the MBFITS header keyword
``HESE`` into ``THESE``, and this code is full of FITS keywords and
abbreviations that a spell checker cannot tell from typos.

Time limits on hung tests
~~~~~~~~~~~~~~~~~~~~~~~~~

The ``SDTmonitor`` tests start background threads and processes, and a
hang there used to stall CI until GitHub cancelled the job after 6 hours,
with no hint of which test was stuck. Now:

* ``pytest-timeout`` (``timeout = 600`` in ``pyproject.toml``) stops any
  single test running longer than 10 minutes. With
  ``timeout_method = "thread"`` it prints the stack of *every* thread and
  ends the pytest process: the rest of the run is lost, but the report shows
  where the background threads are blocked, and a teardown that would hang
  again cannot do so. Override locally with e.g. ``pytest -o timeout=0``.
* ``timeout-minutes: 60`` on the CI test jobs cancels a job after one hour,
  as a backstop.

Warning policy in tests
-----------------------

``DeprecationWarning`` is turned into an error only when it is attributed to
``srttools`` itself (filter ``error::DeprecationWarning:srttools.*`` in
``pyproject.toml``). In practice, this catches calls from our code to
deprecated APIs of numpy, scipy, astropy etc., while ignoring deprecation
warnings raised *between* third-party packages (e.g. Pillow or pyparsing
warnings triggered inside old matplotlib versions), which we cannot fix.

The command line page of the docs
---------------------------------

The page listing the options of every script (``docs/scripts/cli.rst``) is
written again at each docs build by ``docs/conf.py``, which runs every script
with ``--help``. So it always matches the code, also on Read the Docs. Two
details:

* The file is not tracked by git (it is in ``.gitignore``). When it was
  tracked, every local docs build modified it, and the published page showed
  an old copy.
* Each script is run through its function (the one listed in
  ``[project.scripts]`` of ``pyproject.toml``) with the Python interpreter
  that builds the docs. The ``SDT*`` executables found on the ``PATH`` are not
  used: they might belong to another, older installation.

Dates and time zones
--------------------

Dates given on the command line (e.g. ``SDTinspect --only-after
20250101-000000``) are always interpreted as UTC, independent of the time
zone of the machine. Tests force a non-UTC time zone to check this, because
CI machines run in UTC and would not catch the problem otherwise.

Simulated noise and pointing errors
-----------------------------------

The simulator (``srttools.simulate``) is used to test the measurement of
fluxes from cross scans. Its noise model is:

- **White noise** of rms ``noise_amplitude``, independent in each channel and
  spectral bin.
- **1/f noise** (optional, ``fknee`` not None). The total power spectrum is
  :math:`\sigma^2 [1 + (f_{\rm knee}/f)^\alpha]`, where :math:`\sigma` is
  ``noise_amplitude`` and :math:`\alpha` is ``noise_index`` (default 1). It is
  generated in frequency space with random Gaussian amplitudes and phases,
  then transformed back to a time series. The zero-frequency term is set to
  zero (an offset is a baseline, not noise). The same 1/f series is added to
  all spectral bins of a channel, as gain drifts affect the whole band.
- **Correlation between the two channels** (``pol_correlation`` in
  ``sim_crossscans``, default 0.9). Part of the drift (atmosphere, common
  system noise) is shared by the two polarizations, and part (polarized
  emission, receiver noise) is not. Each channel gets
  :math:`\sqrt{\rho}\,n_{\rm common} + \sqrt{1-\rho}\,n_{\rm own}`, which has
  correlation coefficient :math:`\rho` and a noise level independent of it.
  Only the 1/f part is correlated.
- **Pointing offset** (``pointing_offset=(dra, ddec)`` in degrees on the sky,
  in ``sim_crossscans``): the true source is displaced from the center of the
  cross. Along each arm the profile is shifted; across the arm it is
  attenuated by :math:`\exp[-d^2/(2\sigma_{\rm beam}^2)]` (``beam_sigma``,
  default the 2.5' sigma of the simulated beam).

All random draws accept ``rng``: a seed or a ``numpy.random.Generator``.
With ``rng=None`` the global ``numpy.random`` state is used, as before.

Fluxes and upper limits from cross scans (``SDTcross``)
-------------------------------------------------------

``srttools.crossscan`` (command line: ``SDTcross``) measures the flux of a
point source of known position from cross scans. For each subscan and
channel, ``fit_subscan`` fits, in one linear least-squares fit, a cubic
spline baseline and a Gaussian of fixed width (the beam) centred at the
**nominal** source position. Only the amplitude of the Gaussian is free (a
*forced fit*), so it can be negative and is unbiased also when there is no
source. The spline has knots about every 3 beam FWHM, none within 1.5 FWHM
of the source, so that it cannot mimic the source.

The amplitudes of the two channels (LCP and RCP) are averaged subscan by
subscan, so that the noise they share is not counted twice. The flux is the
mean of these amplitudes over the subscans, and its error is their standard
deviation divided by the square root of the number of subscans. Since each
subscan is an independent measurement, this error includes 1/f noise and
baseline drifts, which make the samples within a subscan correlated, and
which an error computed from the per-sample noise would ignore. Subscans
with anomalous noise (residual RFI) are rejected first (see below).

A source is detected if its flux exceeds ``nsigma`` (default 3) times its
error. The upper limit is ``max(F, 0) + nsigma * sqrt(err**2 + (c *
max(F, 0))**2)``, where ``c`` is the relative calibration error (from the
calibrator table, or ``--cal-rel-err``). The arms (RA vs Dec) are compared
with errors added in quadrature; the polarizations (LCP vs RCP) with the
per-subscan difference, which accounts for the noise they share.

The settings were chosen with Monte Carlo simulations reproducing the noise
(white and 1/f), baseline drifts, RFI and pointing errors of real SRT cross
scans (``experiments/crossscan_mc``, section 7 of its notes). All the
variants tested (baseline fitted outside a window around the source and
amplitude measured afterwards, or fitted together with the Gaussian; knots
every 2, 3, 5 FWHM; mean or median) give correct errors: the scatter of the
flux over the simulations matches the reported error within ~12%, and there
are at most 1% false detections at 3 sigma. The joint fit gives errors ~30%
smaller than the window fit; the median costs ~25% more than the mean. Compared
with fitting a Gaussian to the profile averaged over all subscans, with an
error from the per-sample scatter, the real scatter of the flux is 10-40%
smaller, and the reported error is right instead of 2.5-3 times too small.

Limitations:

- Pointing errors and extended sources: the forced fit measures the peak at
  the nominal position. A pointing offset of 0.3 beams on one axis loses
  ~22% of the flux (17% when the baseline is fitted outside a window of 1.5
  FWHM, ``joint=False``); for a source 2 times wider than the beam, the
  joint fit gives 76% of the peak, because its baseline absorbs part of the
  wings. To detect this, ``SDTcross`` also measures each subscan with the
  baseline fitted outside a window of 1.5 FWHM, which is less sensitive to
  the shape of the source, and compares the two fits subscan by subscan
  (``CONS_SHP``, in standard deviations; a warning is printed beyond
  ``--nsigma``). For faint sources the difference is lost in the noise, as
  is the loss itself; for bright sources and calibrators, where a loss of a
  few percent dominates the error, the check becomes sensitive. The profiles
  plotted by ``SDTcross`` show whether the source is centred and as wide as
  the beam.
- Nearby sources: a source on an arm, outside the knot-free region, bends
  the baseline: from -50% to +18% of its flux, depending on its distance (1.5
  to 4.5 FWHM), leaks into the measurement.

Rejection of RFI-contaminated subscans in cross scans
-----------------------------------------------------

Radio-frequency interference (RFI) that survives the removal of bad spectral
channels can remain in the time series as bursts lasting a few seconds, too
long for the spike filter. Before averaging cross scans,
``ScanSet.calculate_avg_cross`` calls
``ScanSet.flag_rfi_subscans``, which uses
``srttools.imager.find_outlier_subscans`` on each channel:

1. for each arm, the median profile of all its subscans, in bins of one
   pixel, is subtracted; this removes the source, however bright;
2. the rms of the residuals is measured in each subscan;
3. subscans whose rms exceeds the median rms by more than ``rfi_nsigma``
   robust standard deviations (1.4826 times the median absolute deviation)
   are flagged: their ``<channel>-filt`` values are set to False.

Every flagged subscan is reported in a warning, with a further warning if
more than 20% of the subscans of a channel are flagged. The default
``rfi_nsigma = 5`` was chosen with Monte Carlo simulations reproducing the
noise and RFI of real SRT cross scans (``experiments/crossscan_mc``): it
flags 0.2-0.3% of the clean subscans with 1/f noise (none with white noise),
including when a bright source is observed with a pointing offset, and
misses only 0.1% of the subscans containing RFI bursts. Set ``rfi_nsigma=None``
(``--rfi-nsigma 0`` in ``SDTcross``) to disable the rejection.

This rejection is meant for cross scans: in a map, subscans at different
offsets see different parts of the sky and the median profile is not a good
model of each of them.

Rejection of outlier fits in calibration and light curves
---------------------------------------------------------

``SDTcal`` and ``SDTlcurve`` fit a Gaussian plus a baseline to each subscan
separately. RFI bursts can bias these fits, or the Gaussian can lock on to
the burst instead of the source. ``CalibratorTable.update`` therefore calls
``CalibratorTable.flag_outliers``, which compares each fit with the other
fits of the same observation (``Dir``), source, channel and arm (RA or Dec),
where the source is not expected to change. A fit is flagged in the new
``Outlier`` column if its residual rms (new ``Resid Std`` column; only
excesses count), amplitude or width is more than ``outlier_nsigma`` robust
standard deviations (1.4826 times the median absolute deviation) from the
median of its group. Groups with fewer than 6 fits are not checked. Flagged
fits are excluded from the calibration, from the source fluxes and from the
light curves written by ``SDTlcurve``, and are reported in warnings.

The default ``outlier_nsigma = 5`` was chosen with the Monte Carlo
simulations of ``experiments/crossscan_mc``. For a bright source (peak ~40
times the per-sample noise) it flags 0.4% of the clean fits; with RFI it
brings the flux back to its RFI-free value (without rejection, the RFI model
of the simulations biases the flux by ~20%). For faint sources (below S/N
~100 in a whole observation) many single-subscan fits lock on to noise:
20-45% of them are flagged, which reduces the bias of the mean flux but does
not remove it; the warning that more than 20% of the fits are flagged then
signals that the fits are unreliable. Set ``outlier_nsigma=None``
(``--outlier-nsigma 0``) to disable the rejection.

Tables written by previous versions are read with ``Resid Std`` set to NaN
(no fit is flagged for its residuals) and ``Outlier`` set to False.

Constant calibration and its error
----------------------------------

When the calibrators cannot constrain a line in elevation (fewer than 3
fits), or the target lies outside the elevation range of the calibrators
(by more than half that range), the calibration factor (Jy per count) is a
constant, from ``CalibratorTable.Jy_over_counts_rough``. It rejects fits more
than 5 formal errors away from the mean and averages the others, weighting
them with their formal errors. If fewer than 5 fits would survive (few fits,
or fits that scatter well beyond their formal errors, e.g. because of
pointing errors or a bad calibrator scan) it warns "Calibration fit is
shaky" and uses the median of the fits. The error of the median is
1.2533 sigma / sqrt(N), where sigma is the robust standard deviation of the
fits (from the median absolute deviation), or their median formal error if
larger. Before this, the error was the median formal error of one fit.

On a Medicina cross scan of NGC7027 (2022-12-13, 8.18 GHz, 14 fits per
channel scattering by 26-41% with formal errors of ~1%), the error of the
calibration factor goes from ~1% to 12% (RCP) and 18% (LCP), and the flux
densities of Cas A from the two polarizations, which differ by 20%, become
consistent.

Calibrator maps and off-axis feeds
----------------------------------

``SDTcal`` fits each subscan as if it crossed the calibrator. This is true
for the central feed of a cross scan, but not for most subscans of a
calibrator map, nor for the off-axis feeds of a multi-feed receiver. The
fit of a subscan that misses the source gives a small or zero amplitude, and
a Flux/Counts many times too large. Therefore, a channel is fitted only if
its feed passes within 1/4 of the beam FWHM (1.22 lambda / D at the centre
of the band) of the catalogue position of the calibrator; there the peak is
attenuated by at most 16% (by ~3% for the closest subscan of a map with the
usual spacing of 1/3 of a beam).

The Gaussian fit starts at the point of the subscan closest to the
calibrator, with the width of the beam. Starting at the middle of the
subscan, as ``fit_baseline_plus_bell`` does by default, the fit converges to
zero amplitude when a narrow source is far from the middle, as in most
subscans of a map.

With few calibrator fits per channel (e.g. one per map, or the two arms of a
single cross scan) a line in elevation cannot be fitted: with fewer than 3
valid fits the conversion factor is the constant from
``Jy_over_counts_rough``.

On a K-band, 7-feed session (2018-07-24) with five 3C147 maps as
calibrators, these changes bring the flux density of 3C84 from ~235 Jy to
23.7 Jy, consistent with an independent 2D fit of the same maps (~25 Jy,
with the same 10% map-to-map scatter).

Excluded regions in the baseline subtraction
--------------------------------------------

``Scan.baseline_subtract`` (and ``SDTimage --exclude``) can exclude regions
of the sky, e.g. an extended source, from the baseline fit. Subscans with at
least 300 samples are fitted on a median-filtered (31 samples) version of
the data: the samples in the excluded regions are first replaced by a linear
interpolation of the others, so that the filtered data contain no source.
The same replacement is done inside the asymmetric least squares fit. The
final zero-level correction (the median of the residuals) uses only the
samples outside the excluded regions, provided there are at least 20 of them.

On the W44 maps of 2016-02-25 (K band, ~870 samples per subscan), before
this fix the exclusion had no effect: the images were identical with and
without it, and the sky along the scans crossing W44 was 0.3-0.6 times the
pixel noise below zero, on average. With the exclusion, the average is zero.
Wide exclusions leave few samples for the baseline: the subscan-to-subscan
offsets (stripes) along those scans become larger.

Pixel size and empty pixels
---------------------------

Each image pixel is the mean of the samples that fall in it. A pixel with no
samples has value 0 (and NaN standard deviation). If the pixel is smaller
than the spacing between subscans, whole columns (or rows) of pixels inside
the map get no samples: they show up as dark stripes through a bright source,
and any sum over the image, e.g. the flux density of an extended source in
Jy/pixel, is too low. The pixel should be at least as large as the subscan
spacing.

``ScanSet.calculate_images`` warns when more than 1% of the filled pixels are
empty pixels inside the map: empty pixels with filled pixels on both sides
along one axis that also lie within the map along the other axis. Gaps in the
ragged outer rows of a map, where subscans end at slightly different
positions, do not count.

On a Medicina map of Cas A at 8.18 GHz (2022-12-13; subscans 1.25' apart),
1' pixels leave 24% of the pixels empty and the flux density in a 12' aperture
comes out 19% lower than with 1.5' pixels, which leave none.

.. _image-errors:

Image maps: Sdev, ERR and calibration
-------------------------------------

What each map contains
~~~~~~~~~~~~~~~~~~~~~~

For each channel (e.g. ``Feed0_LCP``), ``SDTimage`` saves several maps on the
same pixel grid. If a pixel contains the :math:`N` samples
:math:`x_1, \dots, x_N`:

+ ``<channel>``: the mean, :math:`\bar{x} = \frac{1}{N}\sum_i x_i`;
+ ``<channel>-EXPO``: the number of samples, :math:`N`;
+ ``<channel>-Sdev``: the standard deviation of the samples,
  :math:`s = \sqrt{\frac{1}{N-1}\sum_i (x_i - \bar{x})^2}`. It is NaN when
  :math:`N < 2`, where a scatter cannot be measured;
+ ``<channel>-ERR``: the error of the mean, :math:`s/\sqrt{N}`, plus, in
  calibrated maps, the calibration error (see below);
+ ``<channel>-Outliers``: a score that is larger when the samples of the
  pixel are more likely to contain outliers.

Sdev describes the samples, ERR describes the pixel value: with more samples
per pixel, Sdev stays about the same and ERR decreases.

Two ways of calibrating
~~~~~~~~~~~~~~~~~~~~~~~

The calibration (``SDTimage --calibrate``) converts counts to Jy/beam or
Jy/pixel with a factor :math:`f` (Jy per count) that depends on elevation,
with a relative error :math:`\epsilon = \sigma_f / f`.

+ **Per-sample calibration** (default; ``calibrate_scans=True`` in
  ``ScanSet.calculate_images``): each sample is multiplied by the factor at
  its own elevation, then the samples are binned. Sdev is the standard
  deviation of the calibrated samples.
+ **Quick calibration** (``SDTimage --quick``;
  ``ScanSet.calibrate_images``): the images are made in counts and then
  multiplied by a single factor, the one at the mean elevation of the map.
  Sdev is multiplied by :math:`|f|` (a scatter is never negative). The
  uncalibrated maps are kept with ``-RAW`` in their name (e.g.
  ``Feed0_LCP-RAW-Sdev``). It is faster, but wrong if the gain changes over
  the elevation range of the map.

In both cases the calibration error of a pixel with value :math:`F` (in Jy)
is :math:`\epsilon F`, added in quadrature to the statistical error:

.. math::

    \mathrm{ERR} = \sqrt{\left(\frac{s}{\sqrt{N}}\right)^2 + (\epsilon F)^2}

(with per-sample calibration, :math:`\epsilon` is the average over the
samples).

The calibration error is the *same* error for the whole map: if the factor is
5% too high, every pixel is 5% too high. It does not average down when pixels
are added. For example, take the flux density of an extended source as the
sum of 100 pixels of 0.1 Jy/pixel, each with a statistical error of 0.02 Jy,
and :math:`\epsilon = 5\%`. The total is 10 Jy; its statistical error is
:math:`\sqrt{100} \times 0.02 = 0.2` Jy, and its calibration error is
:math:`5\% \times 10 = 0.5` Jy, giving :math:`\sqrt{0.2^2 + 0.5^2} = 0.54`
Jy. Adding the ERR of the pixels in quadrature would instead give only
:math:`\sqrt{100} \times 0.005 = 0.05` Jy of calibration error. So, for sums
over pixels, take the statistical part from Sdev and EXPO
(:math:`s/\sqrt{N}`), sum it in quadrature, and add the calibration error of
the total once.

Combined maps
~~~~~~~~~~~~~

+ **Destriped maps** (``SDTimage --destripe``) combine the map of the RA
  subscans and the map of the Dec subscans. Their ERR is the error of the mean
  of the two, weighted with the number of samples:
  :math:`\sqrt{(N_h e_h)^2 + (N_v e_v)^2} / (N_h + N_v)`, where :math:`e_h`
  and :math:`e_v` are the ERR of the two maps. Where one direction has a
  single sample (no Sdev), the ERR of the map of all subscans is used.
+ **Maps averaged over channels** (``SDTimage --scrunch-channels``):
  ``TOTAL-Sdev`` and ``TOTAL-ERR`` are :math:`\sqrt{\sum_c \sigma_c^2} / n`,
  for :math:`n` channels, as for the mean of independent measurements.
  Channels of the same feed share the atmosphere and the baseline drifts, so
  they are not fully independent: the true error is somewhat larger.

What ERR does not include
~~~~~~~~~~~~~~~~~~~~~~~~~

ERR is computed only from the scatter of the samples *inside* each pixel.
Tests on real SRT maps (3C84 and empty sky at K band, 2018; W44, 2016) show
what it misses:

+ **Noise correlated along the scans.** Baseline drifts and residual
  offsets between subscans (stripes) move all the samples of a subscan
  together, so they do not appear in the scatter inside a pixel. On empty
  sky (K band), the robust standard deviation of value / ERR is 1.38-1.39, not 1
  (robust standard deviation: 1.4826 times the median absolute deviation from
  the median, which ignores the tails). The difference between the W44 maps
  made with RA and with Dec subscans, divided by its ERR, has a robust
  standard deviation of 1.60. The scatter of the row (or column) averages is
  2-6 times what uncorrelated noise of the size of ERR would give.
  On W44, destriping with the source excluded from the baseline fit
  (``--exclude``) removes most of the stripes, reduces the sky noise by 14%,
  and the empty-sky ratio from 1.50 to 1.38.
+ **Pixels with few samples.** With 3-6 samples, Sdev itself is very
  uncertain: some pixels, whose samples agree by chance, get a tiny ERR.
  Because of these, the plain standard deviation of value / ERR on empty sky
  is 2.1, instead of 1.38.
+ **Bright sources.** Near the peak of a source, the samples in a pixel
  differ because the beam response changes across the pixel, not because of
  noise: ERR becomes large (up to 80 times the empty-sky value on the peak of
  3C84). It is not the error of the peak value.
+ **Systematic effects.** Consecutive maps of 3C84 agree to 6-10%, about 30
  times their statistical error: pointing errors, atmosphere and gain changes
  dominate, and no single map can measure them.

In practice:

+ for the noise of a map, use the robust standard deviation of an empty
  region, not ERR;
+ fit sources with uniform weights, not with weights :math:`1/\mathrm{ERR}^2`,
  which would make the fit ignore the core of the source;
+ add the repeatability of repeated observations (or of the calibrators) to
  the error of a flux density.

Where the products are written
------------------------------

There are two output directories, both in the ``[local]`` section of the
config file:

+ ``productdir``: the final products, i.e. the images of ``SDTimage``, the
  calibration tables of ``SDTcal`` and ``SDTlcurve``, the results of
  ``SDTcross``. Their names are built by
  ``srttools.read_config.output_file_name``: the name of the config file
  without ``.ini``, plus a tail (``.fits``, ``_cal.hdf5``, ``_cross.fits``,
  ...). An explicit output name on the command line (``-o``, ``--output``,
  ``--outroot``) is used as it is.

+ ``intermediatedir``: the intermediate products, i.e. the processed scans
  (one HDF5 file per subscan, in the same directory structure as the raw
  data below ``datadir``; raw files outside ``datadir`` go to a subdirectory
  with the name of their own directory), the diagnostic plots (e.g. the
  calibrator fits of ``SDTcal --show``, in ``<subscan>_scanfit`` next to the
  processed subscan), the ``_dump.hdf5`` files of ``SDTimage`` and
  the ``_cal_unfilt.hdf5`` cache of ``SDTcal``. Default: ``productdir``
  followed by ``intermediate``
  (``srttools.read_config.DEFAULT_INTERMEDIATE_SUBDIR``).

``srttools.read_config.read_config`` never returns an empty ``productdir``:
if the config file has none (or ``None``), it sets it to
``srttools_products`` in the current directory
(``srttools.read_config.DEFAULT_PRODUCTDIR``) and gives a warning. Since
every tool gets its settings from ``read_config``, this holds for all of
them, and also for a ``Scan`` loaded without a config file. The reason is
that the raw data directory should never be written to: archives are often
read-only or shared, and files appearing next to the raw data were easy to
mistake for data. Relative paths in the config file start from the directory
of the config file; relative paths on the command line start from the
current directory.

The reason for two directories is that the processed scans are many small
files, useful only to speed up later runs and to keep the edits made by hand,
while the final products are the few files the user looks at.

Older versions behaved differently: the default was to write next to the raw
data (HDF5 files left there are not read any more); the processed scans were
written directly in ``productdir``; images, calibration tables and dumps were
written next to the config file. With an old config file the processed scans
are therefore not found and are made again. ``read_config`` warns when the
config file has no ``intermediatedir``, the default one does not exist yet,
and ``productdir`` contains HDF5 files at most one directory level down
(files named after the config file, i.e. final products, are ignored). The
warning points to ``SDTupgradeconfig`` (``srttools.read_config.upgrade_config``),
which adds ``intermediatedir``, equal to ``productdir``, right after the
``productdir`` line, editing the file as text so that comments survive, and
keeps the original as ``<file>.bak``.

``read_config(fname, productdir=..., intermediatedir=...)`` replaces the
directories of the config file. ``read_config`` keeps the last config it read
in memory, and the replacements are stored there, so the scans loaded later
in the same run, which call ``read_config`` again, get them too. If only
``productdir`` is replaced, and the config file has no ``intermediatedir``,
the default ``intermediatedir`` moves with it. The command-line options
``--productdir`` and ``--intermediatedir`` of ``SDTpreprocess``,
``SDTimage``, ``SDTcal``, ``SDTlcurve``, ``SDTcross`` and ``SDTrfistat`` are
added by
``srttools.read_config.add_directory_arguments`` and applied, before any scan
is loaded, by ``apply_directory_arguments``. When the input is a dump file
without a config file, there is no config to change: the directories are
written in the ``meta`` of the loaded ``ScanSet`` instead, which is where
``save_ds9_images``, ``SDTcross`` and the interactive filter look for them.
Dump files of older versions have no ``productdir`` in their ``meta``: their
products are written next to the config file (``SDTimage``) or the dump
(``SDTcross``), as before.

The name of the processed file of a raw file is given by
``srttools.scan.processed_file_name``, used by the interactive filter (to
save) and by ``SDTrfistat -c`` (to find the processed scans of the
directories of the config file; it writes its RFI tables and plots to
``productdir``).

Interactive filter web page (``--interactive``)
-----------------------------------------------

The interactive filter has two layers:

+ ``srttools.filter_session.FilterSession`` holds the edits (zapped ranges,
  flags, baseline corrections, alignments) as a list, and has no graphics.
  Each edit changes the channel columns of the ``ScanSet`` (and their
  ``-filt`` masks, ``True`` for good samples) at once, so that new images can
  be computed straight away. A copy of the original columns is kept: *undo*
  restores it and repeats all edits but the last one. Nothing is written to
  disk until ``save``, which copies the edited columns into the ``Scan`` of
  each edited subscan and writes it to a directory: the one given by the
  caller ("save as"), or the intermediate product directory. There is always
  somewhere to save: without a config file, the default is
  ``srttools_products/intermediate`` in the current directory.
  ``check_intermediatedir`` refuses a directory that cannot be
  created or written to, and one where some scan would end up in the
  directory of its raw data file, so the raw data directory is never written
  to. ``config_hint`` gives the line to put in the config file when the
  directory is not its ``intermediatedir``: otherwise later runs would not find
  the edited scans.

+ ``srttools.web_filter.FilterApp`` builds the page with Bokeh and calls the
  session. ``serve`` runs a Bokeh server (a small web server that keeps the
  page and the Python program in contact, so that clicks run Python code).
  Every browser page gets its own widgets, but all of them share one session,
  so reloading the page loses nothing. The server listens only on
  ``localhost``: from another computer it can only be reached through an ssh
  tunnel, which requires a login on the machine. There is no password: other
  users logged in to the same machine can open the page too.

The page has two views, chosen with the ``view`` argument of ``FilterApp``
and ``serve``:

+ ``"map"`` (``SDTimage --interactive``): the map, and the subscans through
  the clicked pixel.

+ ``"scans"`` (``SDTpreprocess --interactive``): no map is calculated. A
  table lists all the subscans (``FilterSession.scan_summary``: name,
  direction, fraction of masked samples in the current channel, flagged or
  not) and the subscan chosen in the table is plotted. It uses the same two
  plots as the map view, showing only the one of the direction of the
  subscan, so that all the editing code is shared. There is no "Apply" (no
  map to recompute) and no "Align all" (it needs the other subscans through
  a pixel). Since here unsaved edits have no other use, "Close and continue"
  with unsaved edits only warns the first time, and closes the second time.

Both views have the "Save to" box, filled with the default directory when the
page opens and checked with ``check_intermediatedir`` every time it changes, so a
directory that cannot be used is reported before any work is lost.

The map looks much as in the old matplotlib display: the ds9-like logarithmic
scale (``srttools.utils.ds9_like_log_scale``), colours from the 20th
percentile to the maximum, and the ``gnuplot2`` colour map, whose blue low end
shows faint stripes better than darker colour maps. The "Outliers" view shows
the ``-Outliers`` image when the pixels have more than 6 samples on average,
and the standard deviation image otherwise, on the same scale, after setting
empty pixels to the mean of the others.

One difference from the old display: pixels with no samples do not enter the
scaling or the percentile. They are drawn at the minimum of the filled pixels.
In the old display they count as zeros, so in a map made without baseline
subtraction (raw counts well above zero) they stretch the scale, and the sky
background turns black. With this rule, adding a constant to all the samples
does not change the picture.

Speed. Measured on an Apple M1 Pro laptop, with bokeh 3.9, on the W44
maps of 2016-02-25 (K band, two channels, 267,627 samples in 339 subscans,
map of 180 x 151 pixels):

=========================================  =========
Step                                       Time
=========================================  =========
Load the scan set (HDF5 dump)              0.4 s
Build the page                             54 ms
Click a pixel (5 subscans, 3881 samples)   8 ms
Zap a range and redraw the subscans        9 ms
Apply (recompute the images and the map)   0.35 s
Undo, with 20 edits to repeat              28 ms
=========================================  =========

The first image calculation of a run takes ~3 s longer, because the binning
code is compiled (numba). The subscans through a pixel are sent to the browser
in full (0.1 MB here), with no thinning. On the Medicina map of Cas A (25,000
samples) all steps take less than 30 ms.
