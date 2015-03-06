# For Python 2 and 3 compatibility
try:
    import configparser
except ImportError:
    import ConfigParser as configparser

def read_config(fname):
    '''Read a config file and return a dictionary of all entries'''

    config_output = {}

    Config = configparser.ConfigParser()

    dum = Config.read(fname)

    # Read local information

    local_params = dict(Config.items('local'))
    print(local_params)
    config_output['workdir'] = local_params['work_directory']
    config_output['datadir'] = local_params['data_directory']

    # Read analysis information
    analysis_params = dict(Config.items('analysis'))
    config_output['list_of_directories'] = \
        [s for s in analysis_params['list_of_directories'].splitlines()
         if s.strip()]  # This last instruction eliminates blank lines

    config_output['projection'] = analysis_params['projection']
    config_output['interpolation'] = analysis_params['interpolation']
    config_output['prefix'] = analysis_params['prefix']

    return config_output


def test_read_config():
    '''Test that config file are read.'''
    import os
    curdir = os.path.abspath(os.path.dirname(__file__))
    datadir = os.path.join(curdir, '..', '..', 'TEST_DATASET')

    fname = os.path.join(datadir, 'test_config.ini')

    config = read_config(fname)

    for k in config.keys():
        print("{}: {}".format(k, config[k]))
