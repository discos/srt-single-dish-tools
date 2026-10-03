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

Rejection of RFI-contaminated subscans in cross scans
-----------------------------------------------------

Radio-frequency interference (RFI) that survives the removal of bad spectral
channels can remain in the time series as bursts lasting a few seconds, too
long for the spike filter. Before averaging cross scans,
``ScanSet.calculate_avg_cross`` (used by ``SDTimage --crosses-only``) calls
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
(``--rfi-nsigma 0`` in ``SDTimage``) to disable the rejection.

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
