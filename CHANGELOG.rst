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
