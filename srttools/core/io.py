import astropy.io.fits as fits
from astropy.table import Table
import numpy as np


def detect_data_kind(fname):
    '''Placeholder for function that recognizes data format.'''
    return 'fitszilla'


def print_obs_info_fitszilla(fname):
    '''Placeholder for function that prints out oberving information.'''
    lchdulist = fits.open(fname)
    section_table_data = lchdulist['SECTION TABLE'].data
    sample_rates = section_table_data['sampleRate']

    print('Sample rates:', sample_rates)

    rf_input_data = lchdulist['RF INPUTS'].data
    print('Feeds          :', rf_input_data['feed'])
    print('IFs            :', rf_input_data['ifChain'])
    print('Polarizations  :', rf_input_data['polarization'])
    print('Frequencies    :', rf_input_data['frequency'])
    print('Bandwidths     :', rf_input_data['bandWidth'])

    lchdulist.close()


def read_data_fitszilla(fname):
    '''Open a fitszilla FITS file and read all relevant information.'''
    lchdulist = fits.open(fname)

    section_table_data = lchdulist['SECTION TABLE'].data
    chan_ids = section_table_data['id']

    rf_input_data = lchdulist['RF INPUTS'].data
    feeds = rf_input_data['feed']
    IFs = rf_input_data['ifChain']
    polarizations = rf_input_data['polarization']

    feed_input_data = lchdulist['FEED TABLE'].data
    xoffsets = feed_input_data['xOffset']
    yoffsets = feed_input_data['yOffset']
    relpowers = feed_input_data['relativePower']

    data_table_data = lchdulist['DATA TABLE'].data

    info_to_retrieve = ['time', 'raj2000', 'decj2000', 'az', 'el',
                        'derot_angle']

    new_table = Table()
    for info in info_to_retrieve:
        new_table[info] = data_table_data[info]

    for i in chan_ids:
        new_table['Ch{}'.format(i)] = data_table_data['Ch{}'.format(i)]
        new_table['Ch{}'.format(i)].meta = {'polarization': polarizations[i],
                                            'feed': feeds[i],
                                            'IF': IFs[i],
                                            'xoffset': xoffsets[feeds[i]],
                                            'yoffset': yoffsets[feeds[i]],
                                            'relpower': relpowers[feeds[i]],
                                            }

    lchdulist.close()
    return new_table


def read_data(fname):
    '''Read the data, whatever the format, and return them'''
    kind = detect_data_kind(fname)
    if kind == 'fitszilla':
        return read_data_fitszilla(fname)


def test_open_data_fitszilla():
    '''Test that data are read.'''
    import os
    import matplotlib.pyplot as plt
    curdir = os.path.abspath(os.path.dirname(__file__))
    datadir = os.path.join(curdir, '..', '..', 'TEST_DATASET')

    fname = os.path.join(datadir, '20140603-103246-scicom-3C157',
                         '20140603-103246-scicom-3C157_003_003.fits')
    print_obs_info_fitszilla(fname)
    table = read_data(fname)
    print(table)
    print(table['Ch0'].meta)
    for i in range(2):
        plt.plot(table.field('time'), table.field('Ch{}'.format(i))[:])
    plt.show()

