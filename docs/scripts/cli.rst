Command line interface
======================

SDTbulkchange
-------------

.. code-block:: none

    usage: SDTbulkchange [-h] [-k KEY] [-v VALUE] [--apply-cal-mark] [--recursive]
                         [--debug]
                         [files ...]

    Change all values of a given column or header keyword in fits files

    positional arguments:
      files                 Single files to preprocess

    options:
      -h, --help            show this help message and exit
      -k KEY, --key KEY     Path to key or data column. E.g. "EXT,header,KEY" to
                            change key KEY in the headerin extension EXT;
                            EXT,data,COL to change columnCOL in the data of
                            extension EXT
      -v VALUE, --value VALUE
                            Value to be written
      --apply-cal-mark      Short for -k "DATA TABLE,data,flag_cal" -v 1
      --recursive           Look for file in up to two subdirectories
      --debug               Plot stuff and be verbose


SDTcal
------

.. code-block:: none

    usage: SDTcal [-h] [--sample-config] [--nofilt] [-c CONFIG] [--splat SPLAT]
                  [-o OUTPUT] [--snr-min SNR_MIN] [--show] [--check]
                  [file]

    Load a series of cross scans from a config file and use them as calibrators.

    positional arguments:
      file                  Input calibration file

    options:
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
      --snr-min SNR_MIN     Minimum SNR for calibrator measurements to be
                            considered valid
      --show                Show calibration summary
      --check               Check consistency of calibration


SDTconvert
----------

.. code-block:: none

    usage: SDTconvert [-h] [-f FORMAT] [--test] [--detrend] [--save-locally]
                      [files ...]

    Load a series of scans and convert them to variousformats

    positional arguments:
      files                 Single files to process or directories

    options:
      -h, --help            show this help message and exit
      -f FORMAT, --format FORMAT
                            Format of output files (options: mbfits, indicating
                            MBFITS v. 1.65; mbfitsw, indicating MBFITS v. 1.65
                            wrapped in asingle file for each FEBE; fitsmod
                            (default), indicating a fitszilla with converted
                            coordinates for feed number *n* in a separate COORDn
                            extensions); classfits, indicating a FITS file
                            readable into CLASS, calibrated when possible;sdfits,
                            for the SDFITS convention
      --test                Only to be used in tests!
      --detrend             Detrend data before converting to MBFITS
      --save-locally        Save all data in the current directory, notalongside
                            the original data.


SDTfake
-------

.. code-block:: none

    usage: SDTfake [-h] [-s SOURCE_FLUX] [-n NOISE_AMPLITUDE] [-b BASELINE]
                   [-g GEOMETRY GEOMETRY GEOMETRY GEOMETRY]
                   [--beam-width BEAM_WIDTH] [--spacing SPACING] [-o OUTDIR_ROOT]
                   [--scan-speed SCAN_SPEED] [--integration-time INTEGRATION_TIME]
                   [--spectral-bins SPECTRAL_BINS] [--no-cal] [--sun] [--debug]

    Simulate a single scan or a map with a point source.

    options:
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
      --spectral-bins SPECTRAL_BINS
                            Simulate a spectrum with this number of bins
      --no-cal              Don't simulate calibrators
      --sun                 Simulate a map of the Sun
      --debug               Plot stuff and be verbose


SDTimage
--------

.. code-block:: none

    usage: SDTimage [-h] [--sample-config] [-c CONFIG] [--refilt] [--altaz]
                    [--sub] [--interactive] [--crosses-only]
                    [--calibrate CALIBRATE] [--nofilt] [-g]
                    [-e EXCLUDE [EXCLUDE ...]] [--chans CHANS] [-o OUTFILE]
                    [-u UNIT] [--frame {icrs,altaz,sun}] [--destripe]
                    [--npix-tol NPIX_TOL] [--debug] [--quick] [--scrunch-channels]
                    [--nosave] [--noplot] [--bad-chans BAD_CHANS] [--splat SPLAT]
                    [--bad-intervals BAD_INTERVALS]
                    [file]

    Load a series of scans from a config file and produce a map.

    positional arguments:
      file                  Load intermediate scanset from this file

    options:
      -h, --help            show this help message and exit
      --sample-config       Produce sample config file
      -c CONFIG, --config CONFIG
                            Config file
      --refilt              Re-run the scan filtering
      --altaz               Do images in Az-El coordinates (deprecated in favor of
                            --frame altaz)
      --sub                 Subtract the baseline from single scans
      --interactive         Open the interactive display
      --crosses-only        Only save cross scan results (no images)
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
      --frame {icrs,altaz,sun}
                            Reference frame for the image. One of icrs, altaz, sun
      --destripe            Destripe the image
      --npix-tol NPIX_TOL   Number of pixels with zero exposure to tolerate when
                            destriping the image, or the full row or column is
                            discarded. Default None, meaning that the image will
                            be destriped as a whole
      --debug               Plot stuff and be verbose
      --quick               Calibrate after image creation, for speed (bad when
                            calibration depends on elevation)
      --scrunch-channels    Sum all the images from the single channels into one.
      --nosave              Do not save the hdf5 intermediate files whenloading
                            subscans.
      --noplot              Do not produce diagnostic plots for data processing
      --bad-chans BAD_CHANS
                            Channels to be discarded when scrunching, separated by
                            a comma (e.g. --bad-chans Feed2_RCP,Feed3_RCP )
      --splat SPLAT         Spectral scans will be scrunched into a single channel
                            containing data in the given frequency range, starting
                            from the frequency of the first bin. E.g. '0:1000'
                            indicates 'from the first bin of the spectrum up to
                            1000 MHz above'. ':' or 'all' for all the channels.
      --bad-intervals BAD_INTERVALS
                            Comma-separated list of frequencies to avoid in the
                            analysis, e.g. '5000:5100,5500:5550' will avoid the
                            frequency intervals 5000-5100 and 5500-5550 MHz. Note:
                            if data were already filtered, you need to specify
                            --refilt as well


SDTinspect
----------

.. code-block:: none

    usage: SDTinspect [-h] [-g GROUP_BY [GROUP_BY ...]] [--options OPTIONS] [-d]
                      [--only-after ONLY_AFTER] [--only-before ONLY_BEFORE]
                      [--ignore-suffix IGNORE_SUFFIX]
                      [--ignore-prefix IGNORE_PREFIX] [--save-calibrator-config]
                      directories [directories ...]

    From a given list of directories, read the relevant information and link
    observations to calibrators. A single file is read for each directory.

    positional arguments:
      directories           Directories to inspect

    options:
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
      --ignore-suffix IGNORE_SUFFIX
                            Suffix, or comma-separated list of suffixes, to be
                            removed from source name. E.g. --ignore-suffix
                            _ra,_dec,_k
      --ignore-prefix IGNORE_PREFIX
                            Prefix, or comma-separated list of prefixes, to be
                            removed from source name. E.g. --ignore-prefix
                            ra_,dec_,k_
      --save-calibrator-config
                            Save calibrator config files as if they were targets


SDTlcurve
---------

.. code-block:: none

    usage: SDTlcurve [-h] [-s SOURCE [SOURCE ...]] [--sample-config] [--nofilt]
                     [-c CONFIG] [--splat SPLAT] [-o OUTPUT]
                     [file]

    Load a series of cross scans from a config file and obtain a calibrated curve.

    positional arguments:
      file                  Input calibration file

    options:
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

    usage: SDTmonitor [-h] [-c CONFIG] [--polling] [-p PORT] [-v] [-w WORKERS]
                      directories [directories ...]

    Run the SRT quicklook in a given directory.

    positional arguments:
      directories           Directories to monitor

    options:
      -h, --help            show this help message and exit
      -c CONFIG, --config CONFIG
                            Configuration file
      --polling             Use a platform-independent, polling watchdog
      -p PORT, --port PORT  The port on which the server will be listening
      -v, --verbosity       Set the verbosity level
      -w WORKERS, --workers WORKERS
                            The maximum number of worker processes to spawn


SDTopacity
----------

.. code-block:: none

    usage: SDTopacity [-h] [--tatm TATM] [--tau0 TAU0] [--t0 T0] files [files ...]

    Calculate opacity from a skydip scan and plot the fit results

    positional arguments:
      files        File to inspect

    options:
      -h, --help   show this help message and exit
      --tatm TATM  Atmospheric temperature
      --tau0 TAU0  Initial value for tau (to be fit)
      --t0 T0      Initial value for Tsys (to be fitted)


SDTparselog
-----------

.. code-block:: none

    usage: SDTparselog [-h] [--to-csv] [--list-calon] [files ...]

    Read ACS logs and return useful information

    positional arguments:
      files         Single files to preprocess

    options:
      -h, --help    show this help message and exit
      --to-csv      Save a CSV file with the results
      --list-calon  List files with calibration mark on


SDTpreprocess
-------------

.. code-block:: none

    usage: SDTpreprocess [-h] [-c CONFIG] [--sub] [--interactive] [--nofilt]
                         [--debug] [--plot] [--nosave] [--splat SPLAT]
                         [--bad-intervals BAD_INTERVALS]
                         [-e EXCLUDE [EXCLUDE ...]]
                         [files ...]

    Load a series of scans from a config file and preprocess them, or preprocess a
    single scan.

    positional arguments:
      files                 Single files to preprocess

    options:
      -h, --help            show this help message and exit
      -c CONFIG, --config CONFIG
                            Config file
      --sub                 Subtract the baseline from single scans
      --interactive         Open the interactive display for each scan
      --nofilt              Do not filter noisy channels
      --debug               Be verbose
      --plot                Plot stuff
      --nosave              Do not save the hdf5 intermediate files whenloading
                            subscans.
      --splat SPLAT         Spectral scans will be scrunched into a single channel
                            containing data in the given frequency range, starting
                            from the frequency of the first bin. E.g. '0:1000'
                            indicates 'from the first bin of the spectrum up to
                            1000 MHz above'. ':' or 'all' for all the channels.
      --bad-intervals BAD_INTERVALS
                            Comma-separated list of frequencies to avoid in the
                            analysis, e.g. '5000:5100,5500:5550' will avoid the
                            frequency intervals 5000-5100 and 5500-5550 MHz.Note:
                            if data were already filtered, you need to specify
                            --refilt as well
      -e EXCLUDE [EXCLUDE ...], --exclude EXCLUDE [EXCLUDE ...]
                            Exclude region from global fitting of baseline and
                            baseline subtraction. It can be specified as X1, Y1,
                            radius1, X2, Y2, radius2 in image coordinates or as a
                            ds9-compatible region file in image or fk5 coordinates
                            containing circular regions to be excluded. Currently,
                            baseline subtraction only takes into account fk5
                            coordinates and global fitting image coordinates. This
                            will change in the future.


SDTrfistat
----------

.. code-block:: none

    usage: SDTrfistat [-h] [--threshold THRESHOLD] [-c CONFIG] [--outroot OUTROOT]
                      [files ...]

    Calculate statistics on the RFI filtered out by SDTpreprocess.

    positional arguments:
      files                 List of files produced by SDTimage or SDTpreprocess
                            (HDF5 format).

    options:
      -h, --help            show this help message and exit
      --threshold THRESHOLD
                            Threshold (% from maximum) for RFI flagging
      -c CONFIG, --config CONFIG
                            Config file
      --outroot OUTROOT     Root for output files


