Command line interface
======================

SDTimage
--------

::

    usage: SDTimage [-h] [--sample-config] [-c CONFIG] [--refilt] [--interactive]
                    [--splat SPLAT]

    Load a series of scans from a config file and produce a map.

    optional arguments:
      -h, --help            show this help message and exit
      --sample-config       Produce sample config file
      -c CONFIG, --config CONFIG
                            Config file
      --refilt              Re-run the scan filtering
      --interactive         Open the interactive display
      --splat SPLAT         Spectral scans will be scrunched into a single channel
                            containing data in the given frequency range, starting
                            from the frequency of the first bin. E.g. '0:1000'
                            indicates 'from the first bin of the spectrum up to
                            1000 MHz above'. ':' or 'all' for all the channels.


.. automodule:: srttools.core.scan

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


.. automodule:: srttools.core.calibration

