Command line interface
======================

SDTcal
------

::

    usage: SDTcal [-h] [--sample-config] [--nofilt] [-c CONFIG] [--splat SPLAT]
                  [-o OUTPUT] [--show]
                  [file]

    Load a series of scans from a config file and produce a map.

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



SDTimage
--------

::

    usage: SDTimage [-h] [--sample-config] [-c CONFIG] [--refilt] [--sub]
                    [--interactive] [--calibrate CALIBRATE] [--nofilt] [-g]
                    [-e EXCLUDE [EXCLUDE ...]] [--chans CHANS] [-o OUTFILE]
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
      --sub                 Subtract the baseline from single scans
      --interactive         Open the interactive display
      --calibrate CALIBRATE
                            Calibration file
      --nofilt              Do not filter noisy channels
      -g, --global-fit      Perform global fitting of baseline
      -e EXCLUDE [EXCLUDE ...], --exclude EXCLUDE [EXCLUDE ...]
                            Exclude region from global fitting of baseline
      --chans CHANS         Comma-separated hannels to include in global fitting
                            (Ch0, Ch1, ...)
      -o OUTFILE, --outfile OUTFILE
                            Save intermediate scanset to this file.
      --splat SPLAT         Spectral scans will be scrunched into a single channel
                            containing data in the given frequency range, starting
                            from the frequency of the first bin. E.g. '0:1000'
                            indicates 'from the first bin of the spectrum up to
                            1000 MHz above'. ':' or 'all' for all the channels.



SDTlcurve
---------

::

    usage: SDTlcurve [-h] [--sample-config] [-c CONFIG]
                     [--pickle-file PICKLE_FILE] [--splat SPLAT] [--refilt]

    Load a series of scans from a config file and produce a map.

    optional arguments:
      -h, --help            show this help message and exit
      --sample-config       Produce sample config file
      -c CONFIG, --config CONFIG
                            Config file
      --pickle-file PICKLE_FILE
                            Name for the intermediate pickle file
      --splat SPLAT         Spectral scans will be scrunched into a single channel
                            containing data in the given frequency range, starting
                            from the frequency of the first bin. E.g. '0:1000'
                            indicates 'from the first bin of the spectrum up to
                            1000 MHz above'. ':' or 'all' for all the channels.
      --refilt              Re-run the scan filtering



