v1.0.0alpha (2026-11-19)
------------------------

Breaking Changes
^^^^^^^^^^^^^^^^

- Python >= 3.12 is now required, and minimum versions are set for all dependencies (numpy 1.26, astropy 5.3.4, scipy 1.11.2, matplotlib 3.7.3, h5py 3.10, pyyaml 6.0.1; pytest 7.3.2 for the tests). The unused ``stingray`` dependency was dropped, and ``tqdm`` (progress bars) is now optional. (`#279 <https://github.com/discos/srt-single-dish-tools/pull/279>`__)
- ``SDTimage --crosses-only`` has been removed: use the new ``SDTcross`` instead. ``SDTimage`` no longer saves cross-scan profiles.

  Image ``-Sdev`` maps are now the standard deviation of the samples in each pixel normalized by N - 1 (it was N), NaN where a pixel has fewer than two samples (it was zero or, by rounding, NaN). Calibrated Sdev maps no longer include the calibration error: with quick calibration (``calibrate_images``) they are the raw Sdev times the absolute value of the calibration factor (the old error formula added relative errors linearly and was negative for negative pixels); with ``calibrate_scans`` they are the standard deviation of the calibrated samples (they were multiplied by 1 + the relative calibration error). The calibration error is now in the new ERR maps.

  The old matplotlib interactive filter has been removed: the module ``srttools.interactive_filter``, ``ScanSet.interactive_display``, ``ScanSet.update_scan``, ``ScanSet.find_scans_through_pixel``, ``ScanSet.reprocess_scans_through_pixel``, ``Scan.interactive_filter``, and the ``interactive`` argument of ``Scan``. ``SDTimage --interactive`` and ``SDTpreprocess --interactive`` now open the Bokeh web page only, and stop with an error if bokeh is not installed.

  Products are never written next to the raw data by default. Without a ``productdir`` in the config file (or with ``productdir : None``), all tools write processed scans and plots to ``srttools_products`` in the current directory, with a warning. HDF5 files left next to the raw data by older versions are no longer read: copy them to a directory with the same structure and set ``productdir`` to it, or let the scans be processed again.

  Final products (images of ``SDTimage``, calibration tables of ``SDTcal`` and ``SDTlcurve``, results of ``SDTcross``) are written to ``productdir`` instead of next to the config file (or the dump file). Processed scans and ``SDTimage`` dump files go to the new ``intermediatedir``, by default ``productdir/intermediate``: with old config files, the scans processed by older versions directly in ``productdir`` are not found (a warning says so). Run ``SDTupgradeconfig`` on the config file to keep using them.

  Processed scans and diagnostic plots mirror the directory structure of the raw data below ``datadir``; before, below the directory that the raw files have in common with ``workdir``, which gave long paths (the whole absolute path of the data, when they were outside ``workdir``) and failed on Windows with data and ``workdir`` on different drives. Raw files outside ``datadir`` go to a subdirectory named as their own directory. When ``workdir`` and ``datadir`` differ, the scans processed by older versions are in other subdirectories, and are processed again. (`#280 <https://github.com/discos/srt-single-dish-tools/pull/280>`__)


New Features
^^^^^^^^^^^^

- Enhanced quicklook webpage so that it allows the user to select the main feed to be displayed on screen (`#256 <https://github.com/discos/srt-single-dish-tools/pull/256>`__)
- Multiple monitor fix and additions.
  SDTmonitor is now fully compatible with Windows OSes.
  It also can display a dynamic number of feeds and channels and is not limited to 14 channels as it was previously.
  The PollingObserver used internally has been replaced with PollingObserverVFS (*Virtual File System* like NFS),
  which is our most common use case for polling.
  The tests coverage has also been drastically increased. (`#266 <https://github.com/discos/srt-single-dish-tools/pull/266>`__)
- Power spectrum added to preprocessing diagnostics (`#271 <https://github.com/discos/srt-single-dish-tools/pull/271>`__)
- Added the ``basename_only`` configuration keyword.

  The said keyword allow to write output files directly inside the desired
  productdir avoiding replicating the full original path. (`#274 <https://github.com/discos/srt-single-dish-tools/pull/274>`__)
- The simulator accepts a random seed or generator (``rng``), 1/f noise with a configurable knee frequency and partial correlation between channels, and a pointing offset for cross scans.

  ``ScanSet.calculate_avg_cross`` now flags whole subscans with anomalous noise, e.g. residual RFI, with a warning; the threshold is set with ``rfi_nsigma`` (default 5, ``None`` to disable).

  ``SDTcal`` and ``SDTlcurve`` flag fits whose residual rms, amplitude or width are outliers with respect to the other subscans of the same observation, source, channel and arm (e.g. RFI), with a warning, and exclude them from calibration and fluxes; the threshold is set with ``--outlier-nsigma`` (default 5, 0 to disable). The calibrator table has two new columns, ``Resid Std`` and ``Outlier``.

  New tool ``SDTcross`` (module ``srttools.crossscan``): flux, error and upper limit of a point source from cross scans. In each subscan a smooth baseline is removed and the amplitude of a Gaussian beam at the source position is measured; the flux is the average over subscans and its error comes from the scatter between subscans, which is correct also with 1/f noise and baseline drifts. Results per channel, arm and in total, consistency checks between arms and between polarizations, profiles and plots.

  New ``<channel>-ERR`` image maps: the error of the mean in each pixel, Sdev / sqrt(N), with the relative error of the calibration factor times the flux added in quadrature for calibrated maps (the same with and without ``calibrate_scans``). The calibration error is the same for the whole map: it does not average down when pixels are summed. They are saved together with the Sdev maps (``save_sdev``); after destriping, they are the error of the exposure-weighted mean of the horizontal and vertical maps.

  ``ScanSet.calculate_images`` (and ``SDTimage``) warns when more than 1% of the pixels inside the map have no samples, which happens when the pixel is smaller than the spacing between subscans: those pixels are set to 0 and bias sums over the image, e.g. the flux density of an extended source in Jy/pixel.

  ``SDTpreprocess --interactive`` opens a web page (Bokeh) listing all subscans, to inspect and clean them one by one without a map: zap ranges, flag, unflag, fit a baseline, undo. New option ``--port``. It needs the config file; given single files too, it lists only those subscans.

  The interactive filter page has a "Save to" box: edits can be saved to any directory ("save as") except the raw data one. It always starts with a usable default, and a directory that cannot be used is reported at once. ``FilterSession.save`` takes an ``intermediatedir`` argument and returns the directory used.

  New ``intermediatedir`` in the ``[local]`` section of config files (default: ``productdir/intermediate``), for the intermediate products: processed scans, diagnostic plots, ``SDTimage`` dump files, the ``SDTcal`` cache. ``SDTpreprocess``, ``SDTimage``, ``SDTcal``, ``SDTlcurve``, ``SDTcross`` and ``SDTrfistat`` take ``--productdir`` and ``--intermediatedir`` from the command line. New ``SDTupgradeconfig`` adds ``intermediatedir`` to old config files, so that the scans processed by older versions are used again.

  ``SDTcal`` corrects the amplitude of each calibrator subscan for the pointing offset across it, measured by the perpendicular subscans of the same observation: a subscan passing beside the peak gives a lower amplitude, by exp(-d^2/(2 sigma^2)) for an offset d and beam width sigma. The correction factor is in the new ``Offset Corr`` column of the calibration table, and the position of the track across the subscan in ``Track Offset``. Subscans without perpendicular ones (e.g. maps scanned in one direction) are not corrected. Cached calibrator fits of older versions are recalculated.

  ``SDTcross --fit-position``, for bright sources such as calibrators: fits the position of the source along each arm, measures the flux there and corrects it for the pointing offset across each arm, as ``SDTcal`` does for calibrators. The offsets and correction factors are written in the header (``OFF_RA``, ``OFF_DEC``, ``COR_RA``, ``COR_DEC``). New function ``srttools.crossscan.fit_arm_positions``. (`#280 <https://github.com/discos/srt-single-dish-tools/pull/280>`__)
- Interactive filtering of maps (``SDTimage --interactive``) is now a web page built with Bokeh, served locally (also through an ssh tunnel) or shown in a Jupyter notebook. Edits are kept in memory, can be undone, and are saved only to the product directory. Install with ``pip install srttools[interactive]``. (`#281 <https://github.com/discos/srt-single-dish-tools/pull/281>`__)


Bug Fixes
^^^^^^^^^

- Replaced ``commonprefix`` with ``commonpath`` in *scan.py*. (`#265 <https://github.com/discos/srt-single-dish-tools/pull/265>`__)
- Transformed relative paths in absolute in *scan.py product_path_from_file_name*. (`#268 <https://github.com/discos/srt-single-dish-tools/pull/268>`__)
- Fixed broken linkcheck for stackoverflow.com. Now a required cookie is
  retrieved before checking links. The cookie will only be used for requests
  directed to stackoverflow.com. (`#276 <https://github.com/discos/srt-single-dish-tools/pull/276>`__)
- Dates given to ``SDTinspect --only-after/--only-before`` are now always interpreted as UTC, independent of the local time zone.
  Fix import of ``TableMergeError``, which was removed from ``astropy.table.np_utils`` in astropy 8.
  Fix compatibility with numpy 2.5 (setting array shapes is deprecated) and astropy 8 (SDFITS ``TDIM`` given as a tuple).
  Use scipy sparse arrays instead of the deprecated sparse matrices in the ALS baseline fit. (`#279 <https://github.com/discos/srt-single-dish-tools/pull/279>`__)
- SDTcal no longer stops when the Gaussian fit of a single subscan raises an error (e.g. on a strong spike): the subscan is skipped with a warning.

  In the ``<channel>_cal`` columns that ``ScanSet.calculate_avg_cross`` writes (saved in the dump file), flagged samples are now NaN instead of zero.

  srttools no longer forces the matplotlib backend (it used TkAgg whenever DISPLAY was set, Agg otherwise): the backend chosen by the user (``MPLBACKEND``, ``matplotlibrc``) or the platform default is used, e.g. the native macOS one.

  Plots that are only saved to files (debug plots, calibration summaries, ``SDTcross`` profiles, RFI statistics, ...) are now drawn on figures not attached to pyplot: they never open a window with interactive backends. Some were also left open (RFI statistics, destriping) or saved before adding the legend (outlier removal).

  The alt-az diagnostic plots (``delta_altaz.png``, ``altaz_with_src.png``) are written to the current directory only in debug mode.

  Quick calibration of images (``calibrate_images``) no longer gives a flux of exactly one calibration factor to pixels with NaN Sdev; the Sdev is calculated around the pixel mean, so it is never NaN by rounding and keeps its precision with large offsets.

  ``scrunch_images``: ``TOTAL-Sdev`` took the square root twice; it is now the square root of the sum of the squared Sdev of the channels divided by their number, and there is a ``TOTAL-ERR`` computed in the same way.

  Average cross profiles in ``ScanSet`` (``calculate_avg_cross``, used by ``save_ds9_images`` with ``save_crosses``): quick calibration (without ``calibrate_scans``) no longer fails, calibrated profiles are stored under the right names, and Sdev, ERR and the calibration error are computed as for images (new ``err`` column in the saved tables).

  SDTcal on multi-feed receivers: each channel is fitted once, with the coordinates of its own feed. Before, every channel was fitted once per feed (e.g. 7 times with the K-band 7-feed receiver), six times with the coordinates of another feed, which multiplied the calibration and light-curve points and gave wrong fitted positions.

  The calibration (Jy/counts conversion) ignores fits with zero amplitude, which give an infinite Flux/Counts. Before, a single such fit made the elevation fit raise an error for the whole channel (e.g. with calibrator maps, where many subscans miss the source).

  SDTcal with calibrator maps and multi-feed receivers: only the channels whose feed passes within 1/4 of the beam FWHM of the calibrator are fitted; the others (most subscans of a map, the off-axis feeds of a cross scan) gave meaningless Flux/Counts. The Gaussian fit starts at the calibrator position with the beam width: from the middle of the subscan it often ended at zero amplitude. Channels with fewer than 3 calibrator fits get a constant calibration instead of failing the elevation fit (``ZeroDivisionError`` from statsmodels).

  Baseline subtraction (``baseline_als``) now honours the excluded regions (``avoid_regions``, ``SDTimage --exclude``) also in subscans with 300 samples or more, where they were silently ignored; the final zero-level correction uses only the samples outside the excluded regions. Before, the sky around an extended source excluded from the fit ended up below zero ("bowls" along the scans crossing it).

  ``ScanSet.calculate_images`` with ``direction`` (images of one scan direction) and quick calibration returned the uncalibrated images; only the copy stored in ``images_hor``/``images_ver`` was calibrated (the destriping, which uses those, was not affected).

  The destriping (``SDTimage --destripe``, ``ScanSet.destripe_images``) writes its diagnostic images (``<channel>_hor.png``, ``<channel>_destr.png``, ...) to the current directory only in debug mode; they were always written.

  Destriped images: the ERR of pixels where one scan direction has a single sample (whose Sdev is undefined) is the ERR of the image of all subscans, instead of NaN.

  Constant calibration (``Jy_over_counts_rough``, used with fewer than 3 calibrator fits or outside their elevation range): when the fits scatter well beyond their formal errors, or there are fewer than 5, the error of the calibration factor is now the error of the median from the robust scatter of the fits. It used to be the median formal error of a single fit, e.g. 1% instead of 17-27% on a Medicina calibrator scan whose fits scatter by 26-41%. When enough fits survive the rejection of outliers, the error is never smaller than the robust scatter of all the fits divided by the square root of their number: the few fits that survive can agree among themselves much better than the calibration is known (4% instead of 17% on the same scan).

  ``SDTrfistat -c`` looked for the processed scans next to the raw data, where they are no longer written: it now finds them in ``intermediatedir``, and writes its RFI tables and plots to ``productdir``. ``--outroot`` was ignored when the tables were written, and the tables could then not be plotted: fixed.

  ``CalibratorTable.from_scans`` loads the subscans with its ``config_file``: called from Python without reading the config file first, its temporary files went to ``./srttools_products`` regardless of ``productdir``. The calibrator fit plots of ``SDTcal --show`` (``<subscan>_scanfit``) go to ``intermediatedir`` instead of next to the raw data.

  ``ScanSet.read_images_from_meta`` read back only the last saved image of each kind (images, horizontal and vertical scan images): it now reads them all.

  Simulated cross scans (``sim_crossscans``): in the odd crosses, scanned backwards, the positions were reversed but the data were not, so the source was one sample (0.08 arcmin) away from its nominal position, and on the wrong side when a ``pointing_offset`` was given. Positions and data are now reversed together.

  When spectral data are cleaned, only part of the band is kept (by default 10-90% of it, or the ``--splat`` range): the start frequency of the channels now moves with it, not only the bandwidth. Before, the band looked lower than it was (by 150 MHz with the default range on a 1500 MHz band), and calibrator fluxes from the catalog and beam widths were computed at the wrong frequency (about 1% and 2% off at C band).

  ``SDTcal`` reused its cache of calibrator fits (``<config>_cal_unfilt.hdf5``) even when the calibrator directories, ``--splat`` or ``--nofilt`` had changed, silently giving the old calibration. The cache now records them, and is recalculated with a warning when they differ. ``read_config`` also reads a config file again if it was modified after the last reading in the same Python session (e.g. in a notebook).

  The pointing error of RA subscans in the ``SDTcal`` table (``RA err``, plotted by ``SDTcal --show``) is now an offset on the sky, like ``Dec err``: it was the difference of RA coordinates, too large by 1/cos(Dec) (by 1.6 times at Dec 50 degrees). (`#280 <https://github.com/discos/srt-single-dish-tools/pull/280>`__)


Documentation
^^^^^^^^^^^^^

- New "Technical details" page (supported versions, how to run the tests, warning policy); fixed the Python and dependency versions in the installation instructions. (`#279 <https://github.com/discos/srt-single-dish-tools/pull/279>`__)
- The command line page of the docs is now rebuilt at every docs build (also on Read the Docs) from the code being documented, and is no longer tracked by git; fixed the author list of the docs. (`#280 <https://github.com/discos/srt-single-dish-tools/pull/280>`__)


Internal Changes
^^^^^^^^^^^^^^^^

- Use padded number for feed in output images (`#259 <https://github.com/discos/srt-single-dish-tools/pull/259>`__)
- Fixed ``coverage xml`` tox command. (`#264 <https://github.com/discos/srt-single-dish-tools/pull/264>`__)
- Fixed the macOS runner to version 13 in order to have an Intel processor based
  architecture for tests. (`#270 <https://github.com/discos/srt-single-dish-tools/pull/270>`__)
- Introduced the ``fail-fast: false`` strategy on CI tests. This will allow
  runners to go on if tests on another runner failed. This will allow to spot
  OS/Python version related issues from common issues without having to commit
  and wait another test failure on another platform. (`#273 <https://github.com/discos/srt-single-dish-tools/pull/273>`__)
- Modernized packaging, tox and CI (oldest- and development-version test jobs, ruff-based code style checks, time limits on hung tests); applied code style fixes to the whole repository. (`#279 <https://github.com/discos/srt-single-dish-tools/pull/279>`__)
- Require ``jinja2>=3.1.6`` for the documentation build, instead of pinning the
  vulnerable ``jinja2==3.1.3``.

  The old experiment scripts (``experiments/``) are no longer in the repository.

  Publish the package to PyPI automatically when a GitHub release is published.

  Tests that take several seconds and duplicate a faster test are marked ``slow`` and run only with ``--run-slow``: in CI, in a separate Linux coverage job, in parallel with the others, and in the weekly jobs. The coverage reports of all jobs are sent to Codecov together by a new workflow, ``codecov_upload.yml``, which also works for pull requests from forks. All CI jobs print the 20 slowest tests. ``SDTrfistat`` collects the flagged channels 4 times faster (the tables were stacked row by row). The imager test setups no longer draw debug plots, halving their time (the plots are tested elsewhere). (`#280 <https://github.com/discos/srt-single-dish-tools/pull/280>`__)


v0.8 (unreleased)
-----------------

New Features
^^^^^^^^^^^^

- Add date to info string in quicklook (`#252 <https://github.com/discos/srt-single-dish-tools/pull/252>`__)


v0.7.0 (2025-01-24)
-------------------

New Features
^^^^^^^^^^^^

- Save average cross scans in FITS extensions, with uncertainties (`#247 <https://github.com/discos/srt-single-dish-tools/pull/247>`__)
- New script to analyze RFI content of observations, and way to apply filters (`#248 <https://github.com/discos/srt-single-dish-tools/pull/248>`__)


Bug Fixes
^^^^^^^^^

- Fix incompatibility with Astropy 6+ (`#238 <https://github.com/discos/srt-single-dish-tools/pull/238>`__)
- Support new naming for summary file (`#242 <https://github.com/discos/srt-single-dish-tools/pull/242>`__)
- Use rough calibration if the calibrator's elevation was too far from the source's. This avoids a bug where a linear calibration vs elevation curve was extrapolated too far (and even gave negative values in some cases) (`#246 <https://github.com/discos/srt-single-dish-tools/pull/246>`__)


Internal Changes
^^^^^^^^^^^^^^^^

- Update dependencies (`#241 <https://github.com/discos/srt-single-dish-tools/pull/241>`__)


v0.6.4 (2023-10-03)
-------------------

New Features
^^^^^^^^^^^^

- Add AZ/El information to the header of images; create functions for circular statistics (`#235 <https://github.com/discos/srt-single-dish-tools/pull/235>`__)


Internal Changes
^^^^^^^^^^^^^^^^

- Use towncrier for changelog; update test infrastructure; eliminate dependency on pyregion and use Astropy Regions instead (`#236 <https://github.com/discos/srt-single-dish-tools/pull/236>`__)


v0.6.3 (2023-09-08)
-------------------

New features
~~~~~~~~~~~~

-  Make SDTinspect ignore suffixes and prefixes from source names
   `#235 <https://github.com/discos/srt-single-dish-tools/pull/235>`__

v0.6.2 (2023-08-16)
-------------------

Internal changes
~~~~~~~~~~~~~~~~

-  Fix docs and link check tests
   `#232 <https://github.com/discos/srt-single-dish-tools/pull/232>`__

New features
~~~~~~~~~~~~

-  Be more permissive with calibrator names
   `#233 <https://github.com/discos/srt-single-dish-tools/pull/235>`__

v0.6.1 (2023-08-02)
-------------------

Internal changes
~~~~~~~~~~~~~~~~

-  More informative warnings when invalid data are present, and only if
   data are in valid range
   `#229 <https://github.com/discos/srt-single-dish-tools/pull/229>`__
-  Use ``black`` for automatic formatting
-  Update testing infrastructure
-  Optimize ALS baseline algorithm
   `#225 <https://github.com/discos/srt-single-dish-tools/pull/225>`__
-  Suppress calculation of sun coordinates and offsets if feeds not in
   file - Take II
   `#217 <https://github.com/discos/srt-single-dish-tools/pull/217>`__
-  Add creator info to header
   `#213 <https://github.com/discos/srt-single-dish-tools/pull/213>`__
-  Eliminate all pickle files produced by cleaning.
   `#208 <https://github.com/discos/srt-single-dish-tools/pull/208>`__


New features
~~~~~~~~~~~~

-  Refresh scanset structure and inner workings, add observing info
   `#212 <https://github.com/discos/srt-single-dish-tools/pull/212>`__

Bug fixes
~~~~~~~~~

-  Fix interactive display not working with some matplotlib backends
   `#231 <https://github.com/discos/srt-single-dish-tools/pull/231>`__

v0.6 (2021-04-23)
-----------------


New features
~~~~~~~~~~~~

-  Sun coordinates
   `#205 <https://github.com/discos/srt-single-dish-tools/pull/205>`__
-  A robust and working quicklook infrastructure
-  New bandwidth format for SARDARA
   `#105 <https://github.com/discos/srt-single-dish-tools/pull/105>`__
-  Parse logs – update CAL information more easily
   `#103 <https://github.com/discos/srt-single-dish-tools/pull/103>`__
-  SDFITS converter
   `#101 <https://github.com/discos/srt-single-dish-tools/pull/101>`__
-  Allow multiple onoff cycles in a single scan
   `#99 <https://github.com/discos/srt-single-dish-tools/pull/99>`__
-  Better filtering of calibrators
   `#202 <https://github.com/discos/srt-single-dish-tools/pull/202>`__


Internal changes
~~~~~~~~~~~~~~~~

-  Fixes to test infrastructure - eliminate Appveyor and Travis, Send
   info to codecov instead of coveralls
-  Some tweaks for the SDTmonitor


Bug fixes
~~~~~~~~~

-  Code style refresh
   `#204 <https://github.com/discos/srt-single-dish-tools/pull/204>`__

0.5 (2018-07-10)
----------------


Internal changes
~~~~~~~~~~~~~~~~

-  Create Code of Conduct and Contribution Guidelines


New features
~~~~~~~~~~~~

-  Simulator script
   `#54 <https://github.com/discos/srt-single-dish-tools/pull/54>`__
-  Destriping algorithm - basket weaving
   `#56 <https://github.com/discos/srt-single-dish-tools/pull/56>`__
   `#57 <https://github.com/discos/srt-single-dish-tools/pull/57>`__
-  Allow to mask out the zeros in the image
   `#58 <https://github.com/discos/srt-single-dish-tools/pull/58>`__
-  Add mask option for baseline fitting, and ds9 region filtering
   `#59 <https://github.com/discos/srt-single-dish-tools/pull/59>`__
-  Characterize the beam with Zernike coeffs
   `#60 <https://github.com/discos/srt-single-dish-tools/pull/60>`__
-  New subdivision of channels in feed + polarization
   `#61 <https://github.com/discos/srt-single-dish-tools/pull/61>`__
-  Opacity calculation from Skydip scans
   `#62 <https://github.com/discos/srt-single-dish-tools/pull/62>`__
-  Read temperatures and use them for calibration curves
   `#64 <https://github.com/discos/srt-single-dish-tools/pull/64>`__
-  Add Monitoring scripts
   `#88 <https://github.com/discos/srt-single-dish-tools/pull/88>`__
-  Convert the coordinates in multifeed observations
   `#90 <https://github.com/discos/srt-single-dish-tools/pull/90>`__
-  Class (fits) converter
   `#97 <https://github.com/discos/srt-single-dish-tools/pull/97>`__

Documentation
~~~~~~~~~~~~~

-  Draft tutorial notebook
   `#67 <https://github.com/discos/srt-single-dish-tools/pull/67>`__


Internal changes
~~~~~~~~~~~~~~~~

-  Uniform calculation of baseline between cross scans and maps in sims
   `#68 <https://github.com/discos/srt-single-dish-tools/pull/68>`__
-  Make opacity calculation aware of summary files
   `#71 <https://github.com/discos/srt-single-dish-tools/pull/71>`__
-  Better format checks; close open fits file; avoid summary.fits
   `#74 <https://github.com/discos/srt-single-dish-tools/pull/74>`__
-  Fix info not saved
   `#77 <https://github.com/discos/srt-single-dish-tools/pull/77>`__
-  Use the MAD for outlier detection
   `#78 <https://github.com/discos/srt-single-dish-tools/pull/78>`__
-  Better cleaning
   `#80 <https://github.com/discos/srt-single-dish-tools/pull/80>`__
-  New tutorial, including destriping and regions
   `#81 <https://github.com/discos/srt-single-dish-tools/pull/81>`__
-  Better smoothing of rms spectra; better handling of invaliddata
   `#86 <https://github.com/discos/srt-single-dish-tools/pull/86>`__
-  Allow different file formats for debugging
   `#87 <https://github.com/discos/srt-single-dish-tools/pull/87>`__
-  Understand XARCOS files
   `#89 <https://github.com/discos/srt-single-dish-tools/pull/89>`__
-  First working draft of MBFITS converter
   `#91 <https://github.com/discos/srt-single-dish-tools/pull/91>`__
-  Fix detrend tests
   `#92 <https://github.com/discos/srt-single-dish-tools/pull/92>`__
-  Raise if watchdog not installed
   `#93 <https://github.com/discos/srt-single-dish-tools/pull/93>`__
-  Fix center of mass
   `#94 <https://github.com/discos/srt-single-dish-tools/pull/94>`__
   `#95 <https://github.com/discos/srt-single-dish-tools/pull/95>`__


0.4 (2017-09-06)
----------------

-  Rework interactive interface on_click and help
   `#38 <https://github.com/discos/srt-single-dish-tools/pull/38>`__
-  Method to convert times to TDB
   `#41 <https://github.com/discos/srt-single-dish-tools/pull/41>`__
-  Faster imager tests
   `#42 <https://github.com/discos/srt-single-dish-tools/pull/42>`__
-  Apply user filters to scanset
   `#43 <https://github.com/discos/srt-single-dish-tools/pull/43>`__
-  Calibrate scans, not images (slower but more correct)
   `#44 <https://github.com/discos/srt-single-dish-tools/pull/44>`__
-  Cleanup code and increase test coverage (#46)
-  Refactor calibrate (#48)
-  Fix coordinates
   `#50 <https://github.com/discos/srt-single-dish-tools/pull/50>`__
