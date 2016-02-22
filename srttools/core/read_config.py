"""Read the configuration file."""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import os
import glob
# For Python 2 and 3 compatibility
try:
    import configparser
except ImportError:
    import ConfigParser as configparser


SRT_tools_config_file = None
SRT_tools_config = None


def sample_config_file(fname='sample_config_file.ini'):
    """Create a sample config file, to be modified by hand."""
    string = """
[local]
    workdir : .
    datadir : ../../TEST_DATASET

[analysis]
    projection : ARC
    interpolation : spline
    prefix : test_
    list_of_directories :
;;Two options: either a list of directories:
;        dir1
;        dir2
;; or a star symbol for all directories
;         *

;; Coordinates have to be specified in decimal degrees. ONLY use if different
;; from target coordinates!
;    reference_ra : 10.5
;    reference_dec : 5.3

;; Number of pixels, specified as pair x y

    npix : 30 30

;; Channels to save from RFI filtering. It might indicate known strong spectral
;; lines
    goodchans :

;; Percentage of channels to filter out for rough RFI filtering (Spectral data
;; only)
    filtering_factor : 0.3
    """
    with open(fname, 'w') as fobj:
        print(string, file=fobj)


def get_config_file():
    """Get the current config file."""
    return SRT_tools_config_file


def read_config(fname=None):
    """Read a config file and return a dictionary of all entries."""
    global SRT_tools_config_file, SRT_tools_config

    # --- If already read, use existing config ---

    if fname == SRT_tools_config_file and SRT_tools_config is not None:
        return SRT_tools_config

    if fname is None and SRT_tools_config is not None:
        return SRT_tools_config

    # ---------------------------------------------

    config_output = {}

    Config = configparser.ConfigParser()

    if fname is None:
        fname = 'config.ini'

    SRT_tools_config_file = fname
    Config.read(fname)

    # ---------- Set default values --------------------------------------

    config_output['projection'] = 'ARC'
    config_output['interpolation'] = 'linear'
    config_output['workdir'] = './'
    config_output['datadir'] = './'
    config_output['list_of_directories'] = '*'
    config_output['npix'] = '32 32'
    config_output['goodchans'] = None
    config_output['filtering_factor'] = '0'

    # --------------------------------------------------------------------

    # Read local information

    local_params = dict(Config.items('local'))

    config_output.update(local_params)

    # Read analysis information
    analysis_params = dict(Config.items('analysis'))

    config_output.update(analysis_params)

    config_output['list_of_directories'] = \
        [s for s in analysis_params['list_of_directories'].splitlines()
         if s.strip()]  # This last instruction eliminates blank lines

    # If the list of directories is not specified, or if a '*' symbol is used,
    # use glob in the datadir to determine the list

    if config_output['list_of_directories'] in ([], ['*'], '*'):
        config_output['list_of_directories'] = \
            [os.path.split(f)[1]  # return name without path
             for f in glob.glob(os.path.join(config_output['datadir'], '*'))
             if os.path.isdir(f)]  # only if it's a directory

    config_output['npix'] = [int(n) for n in config_output['npix'].split()]
    if config_output['goodchans'] is not None:
        config_output['goodchans'] = \
            [int(n) for n in config_output['goodchans']]

    config_output['filtering_factor'] = \
        float(config_output['filtering_factor'])

    SRT_tools_config = config_output
    return config_output


def test_read_config():
    """Test that config file are read."""
    import os
    curdir = os.path.abspath(os.path.dirname(__file__))
    datadir = os.path.join(curdir, '..', '..', 'TEST_DATASET')

    fname = os.path.join(datadir, 'test_config.ini')

    config = read_config(fname)

    print('\n --- Test config file ---')
    for k in config.keys():
        print("{}: {}".format(k, config[k]))


def test_read_incomplete_config():
    """Test that config file are read."""
    import os
    curdir = os.path.abspath(os.path.dirname(__file__))
    datadir = os.path.join(curdir, '..', '..', 'TEST_DATASET')

    fname = os.path.join(datadir, 'test_config_incomplete.ini')

    config = read_config(fname)

    print('\n --- Test incomplete config file ---')
    for k in config.keys():
        print("{}: {}".format(k, config[k]))
