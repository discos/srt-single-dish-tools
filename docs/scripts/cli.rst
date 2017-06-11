Command line interface
======================

SDTcal
------

::

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



SDTimage
--------

::

    usage: SDTimage [-h] [--sample-config] [-c CONFIG] [--refilt] [--altaz]
                    [--sub] [--interactive] [--calibrate CALIBRATE] [--nofilt]
                    [-g] [-e EXCLUDE [EXCLUDE ...]] [--chans CHANS] [-o OUTFILE]
                    [-u UNIT] [--debug] [--splat SPLAT]
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
                            Exclude region from global fitting of baseline
      --chans CHANS         Comma-separated channels to include in global fitting
                            (Ch0, Ch1, ...)
      -o OUTFILE, --outfile OUTFILE
                            Save intermediate scanset to this file.
      -u UNIT, --unit UNIT  Unit of the calibrated image. Jy/beam or Jy/pixel
      --debug               Plot stuff and be verbose
      --splat SPLAT         Spectral scans will be scrunched into a single channel
                            containing data in the given frequency range, starting
                            from the frequency of the first bin. E.g. '0:1000'
                            indicates 'from the first bin of the spectrum up to
                            1000 MHz above'. ':' or 'all' for all the channels.



SDTinspect
----------

::

    usage: SDTinspect [-h] [-g GROUP_BY [GROUP_BY ...]] [-d]
                      directories [directories ...]

    From a given list of directories, read the relevant information and link
    observations to calibrators. A single file is read for each directory.

    positional arguments:
      directories           Directories to inspect

    optional arguments:
      -h, --help            show this help message and exit
      -g GROUP_BY [GROUP_BY ...], --group-by GROUP_BY [GROUP_BY ...]
      -d, --dump-config-files



SDTlcurve
---------

::

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



SDTpreprocess
-------------

::

    usage: SDTpreprocess [-h] [-c CONFIG] [--sub] [--interactive] [--nofilt]
                         [--debug] [--splat SPLAT]
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



