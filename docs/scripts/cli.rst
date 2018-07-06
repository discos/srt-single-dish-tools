Command line interface
======================

SDTbulkchange
-------------

.. code-block:: none

    usage: SDTbulkchange [-h] -k KEY [-v VALUE] [--debug] [files [files ...]]

    Change all values of a given column or header keyword in fits files

    positional arguments:
      files                 Single files to preprocess

    optional arguments:
      -h, --help            show this help message and exit
      -k KEY, --key KEY     Path to key or data column. E.g. "EXT,header,KEY" to
                            change key KEY in the headerin extension EXT;
                            EXT,data,COL to change columnCOL in the data of
                            extension EXT
      -v VALUE, --value VALUE
                            Value to be written
      --debug               Plot stuff and be verbose


SDTcal
------

.. code-block:: none

    usage: SDTcal [-h] [--sample-config] [--nofilt] [-c CONFIG] [--splat SPLAT]
                  [-o OUTPUT] [--show] [--check]
                  [file]

    Load a series of cross scans from a config file and use them as calibrators.

    positional arguments:
      file                  Input calibration file

    optional arguments:
      -h, --help            show this help message and exit
      --sample-config       Produce sample config file
      --nofilt              Do not filter noisy channels
      -c CONFIG, --config CONFIG
                            Config file
      --splat SPLAT         Spectral scans will be scrunched into a single channel
                            containing data in the given frequency range, starting
                            from the frequency of the first bin. E.g. '0:1000'
                            indicates 'from the first bin of the spectrum up to
                            1000 MHz above'. ':' or 'all' for all the channels.
      -o OUTPUT, --output OUTPUT
                            Output file containing the calibration
      --show                Show calibration summary
      --check               Check consistency of calibration


SDTconvert
----------

.. code-block:: none

    usage: SDTconvert [-h] [-f FORMAT] [--test] [--detrend] [files [files ...]]

    Load a series of scans and convert them to variousformats

    positional arguments:
      files                 Single files to process or directories

    optional arguments:
      -h, --help            show this help message and exit
      -f FORMAT, --format FORMAT
                            Format of output files (options: mbfits, indicating
                            MBFITS v. 1.65; mbfitsw, indicating MBFITS v. 1.65
                            wrapped in asingle file for each FEBE; fitsmod
                            (default), indicating a fitszilla with converted
                            coordinates for feed number *n* in a separate COORDn
                            extensions); classfits, indicating a FITS file
                            readable into CLASS, calibrated when possible
      --test                Only to be used in tests!
      --detrend             Detrend data before converting to MBFITS


SDTfake
-------

.. code-block:: none

    usage: SDTfake [-h] [-s SOURCE_FLUX] [-n NOISE_AMPLITUDE] [-b BASELINE]
                   [-g GEOMETRY GEOMETRY GEOMETRY GEOMETRY]
                   [--beam-width BEAM_WIDTH] [--spacing SPACING] [-o OUTDIR_ROOT]
                   [--scan-speed SCAN_SPEED] [--integration-time INTEGRATION_TIME]
                   [--no-cal] [--debug]

    Simulate a single scan or a map with a point source.

    optional arguments:
      -h, --help            show this help message and exit
      -s SOURCE_FLUX, --source-flux SOURCE_FLUX
                            Source flux in Jy
      -n NOISE_AMPLITUDE, --noise-amplitude NOISE_AMPLITUDE
                            White noise amplitude
      -b BASELINE, --baseline BASELINE
                            Baseline kind: "flat", "slope" (linearly
                            increasing/decreasing), "messy" (random walk) or a
                            number (which gives an amplitude to the random-walk
                            baseline, that would be 20 for "messy")
      -g GEOMETRY GEOMETRY GEOMETRY GEOMETRY, --geometry GEOMETRY GEOMETRY GEOMETRY GEOMETRY
                            Geometry specification: length_ra, length_dec,
                            width_ra, width_dec, in arcmins. A square map of 2
                            degrees would be specified as 120 120 120 120. A
                            cross-like map, 2x2 degrees wide but only along
                            1-degree stripes, is specified as 120 120 60 60
      --beam-width BEAM_WIDTH
                            Gaussian beam width in arcminutes
      --spacing SPACING     Spacing between scans in arcminutes (default 0.5)
      -o OUTDIR_ROOT, --outdir-root OUTDIR_ROOT
                            Output directory root. Here, source and calibrator
                            scans/maps will be saved in outdir/gauss_ra,
                            outdir/gauss_dec, outdir/calibrator1,
                            outdir/calibrator2, where outdir is the outdir root
      --scan-speed SCAN_SPEED
                            Scan speed in arcminutes/second
      --integration-time INTEGRATION_TIME
                            Integration time in seconds
      --no-cal              Don't simulate calibrators
      --debug               Plot stuff and be verbose


SDTimage
--------

.. code-block:: none

    usage: SDTimage [-h] [--sample-config] [-c CONFIG] [--refilt] [--altaz]
                    [--sub] [--interactive] [--calibrate CALIBRATE] [--nofilt]
                    [-g] [-e EXCLUDE [EXCLUDE ...]] [--chans CHANS] [-o OUTFILE]
                    [-u UNIT] [--destripe] [--npix-tol NPIX_TOL] [--debug]
                    [--quick] [--scrunch-channels] [--bad-chans BAD_CHANS]
                    [--splat SPLAT]
                    [file]

    Load a series of scans from a config file and produce a map.

    positional arguments:
      file                  Load intermediate scanset from this file

    optional arguments:
      -h, --help            show this help message and exit
      --sample-config       Produce sample config file
      -c CONFIG, --config CONFIG
                            Config file
      --refilt              Re-run the scan filtering
      --altaz               Do images in Az-El coordinates
      --sub                 Subtract the baseline from single scans
      --interactive         Open the interactive display
      --calibrate CALIBRATE
                            Calibration file
      --nofilt              Do not filter noisy channels
      -g, --global-fit      Perform global fitting of baseline
      -e EXCLUDE [EXCLUDE ...], --exclude EXCLUDE [EXCLUDE ...]
                            Exclude region from global fitting of baseline and
                            baseline subtraction. It can be specified as X1, Y1,
                            radius1, X2, Y2, radius2 in image coordinates or as a
                            ds9-compatible region file in image or fk5 coordinates
                            containing circular regions to be excluded. Currently,
                            baseline subtraction only takes into account fk5
                            coordinates and global fitting image coordinates. This
                            will change in the future.
      --chans CHANS         Comma-separated channels to include in global fitting
                            (Feed0_RCP, Feed0_LCP, ...)
      -o OUTFILE, --outfile OUTFILE
                            Save intermediate scanset to this file.
      -u UNIT, --unit UNIT  Unit of the calibrated image. Jy/beam or Jy/pixel
      --destripe            Destripe the image
      --npix-tol NPIX_TOL   Number of pixels with zero exposure to tolerate when
                            destriping the image, or the full row or column is
                            discarded. Default None, meaning that the image will
                            be destriped as a whole
      --debug               Plot stuff and be verbose
      --quick               Calibrate after image creation, for speed (bad when
                            calibration depends on elevation)
      --scrunch-channels    Sum all the images from the single channels into one.
      --bad-chans BAD_CHANS
                            Channels to be discarded when scrunching, separated by
                            a comma (e.g. --bad-chans Feed2_RCP,Feed3_RCP )
      --splat SPLAT         Spectral scans will be scrunched into a single channel
                            containing data in the given frequency range, starting
                            from the frequency of the first bin. E.g. '0:1000'
                            indicates 'from the first bin of the spectrum up to
                            1000 MHz above'. ':' or 'all' for all the channels.


SDTinspect
----------

.. code-block:: none

    usage: SDTinspect [-h] [-g GROUP_BY [GROUP_BY ...]] [--options OPTIONS] [-d]
                      [--only-after ONLY_AFTER] [--only-before ONLY_BEFORE]
                      directories [directories ...]

    From a given list of directories, read the relevant information and link
    observations to calibrators. A single file is read for each directory.

    positional arguments:
      directories           Directories to inspect

    optional arguments:
      -h, --help            show this help message and exit
      -g GROUP_BY [GROUP_BY ...], --group-by GROUP_BY [GROUP_BY ...]
      --options OPTIONS     Options to be written in config files; they have to be
                            specified as a string defining a dictionary. For
                            example,'{"pixel_size": 0.6, "noise_threshold": 5}'
      -d, --dump-config-files
      --only-after ONLY_AFTER
                            Only after a certain date and time, e.g. ``--only-
                            after 20150510-111020`` to indicate scans done after
                            11:10:20 UTC on May 10th, 2015
      --only-before ONLY_BEFORE
                            Only before a certain date and time, e.g. ``--only-
                            before 20150510-111020`` to indicate scans done before
                            11:10:20 UTC, May 10th, 2015


SDTlcurve
---------

.. code-block:: none

    usage: SDTlcurve [-h] [-s SOURCE [SOURCE ...]] [--sample-config] [--nofilt]
                     [-c CONFIG] [--splat SPLAT] [-o OUTPUT]
                     [file]

    Load a series of cross scans from a config file and obtain a calibrated curve.

    positional arguments:
      file                  Input calibration file

    optional arguments:
      -h, --help            show this help message and exit
      -s SOURCE [SOURCE ...], --source SOURCE [SOURCE ...]
                            Source or list of sources
      --sample-config       Produce sample config file
      --nofilt              Do not filter noisy channels
      -c CONFIG, --config CONFIG
                            Config file
      --splat SPLAT         Spectral scans will be scrunched into a single channel
                            containing data in the given frequency range, starting
                            from the frequency of the first bin. E.g. '0:1000'
                            indicates 'from the first bin of the spectrum up to
                            1000 MHz above'. ':' or 'all' for all the channels.
      -o OUTPUT, --output OUTPUT
                            Output file containing the calibration


SDTmonitor
----------

.. code-block:: none

    usage: SDTmonitor [-h] [-c CONFIG] [--test] directory

    Run the SRT quicklook in a given directory.

    positional arguments:
      directory             Directory to monitor

    optional arguments:
      -h, --help            show this help message and exit
      -c CONFIG, --config CONFIG
                            Config file
      --test                Only to be used in tests!


SDTopacity
----------

.. code-block:: none

    usage: SDTopacity [-h] [--tatm TATM] [--tau0 TAU0] [--t0 T0] files [files ...]

    Calculate opacity from a skydip scan and plot the fit results

    positional arguments:
      files        File to inspect

    optional arguments:
      -h, --help   show this help message and exit
      --tatm TATM  Atmospheric temperature
      --tau0 TAU0  Initial value for tau (to be fit)
      --t0 T0      Initial value for Tsys (to be fitted)


SDTpreprocess
-------------

.. code-block:: none

    usage: SDTpreprocess [-h] [-c CONFIG] [--sub] [--interactive] [--nofilt]
                         [--debug] [--splat SPLAT] [-e EXCLUDE [EXCLUDE ...]]
                         [files [files ...]]

    Load a series of scans from a config file and preprocess them, or preprocess a
    single scan.

    positional arguments:
      files                 Single files to preprocess

    optional arguments:
      -h, --help            show this help message and exit
      -c CONFIG, --config CONFIG
                            Config file
      --sub                 Subtract the baseline from single scans
      --interactive         Open the interactive display for each scan
      --nofilt              Do not filter noisy channels
      --debug               Plot stuff and be verbose
      --splat SPLAT         Spectral scans will be scrunched into a single channel
                            containing data in the given frequency range, starting
                            from the frequency of the first bin. E.g. '0:1000'
                            indicates 'from the first bin of the spectrum up to
                            1000 MHz above'. ':' or 'all' for all the channels.
      -e EXCLUDE [EXCLUDE ...], --exclude EXCLUDE [EXCLUDE ...]
                            Exclude region from global fitting of baseline and
                            baseline subtraction. It can be specified as X1, Y1,
                            radius1, X2, Y2, radius2 in image coordinates or as a
                            ds9-compatible region file in image or fk5 coordinates
                            containing circular regions to be excluded. Currently,
                            baseline subtraction only takes into account fk5
                            coordinates and global fitting image coordinates. This
                            will change in the future.


