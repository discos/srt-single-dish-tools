from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import astropy.io.fits as fits
from astropy.table import Table
import numpy as np
import astropy.units as u


def detect_data_kind(fname):
    '''Placeholder for function that recognizes data format.'''
    if fname.endswith('.hdf5'):
        return 'hdf5'
    else:
        return 'fitszilla'


def correct_coordinates(ra, dec, derot_angle, xoffset, yoffset):
    '''Correct coordinates for feed position.

    Uses the metadata in the channel columns xoffset and yoffset'''

    # Clockwise rotation of angle derot_angle
    new_ra = ra + \
        xoffset * np.cos(derot_angle) - yoffset * np.sin(derot_angle)
    new_dec = dec + \
        xoffset * np.sin(derot_angle) + yoffset * np.cos(derot_angle)

    return new_ra, new_dec


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
    frequencies = rf_input_data['frequency']
    bandwidths = rf_input_data['bandWidth']

    feed_input_data = lchdulist['FEED TABLE'].data
    xoffsets = feed_input_data['xOffset']
    yoffsets = feed_input_data['yOffset']
    relpowers = feed_input_data['relativePower']

    data_table_data = lchdulist['DATA TABLE'].data

    info_to_retrieve = ['time', 'az', 'el', 'derot_angle']

    new_table = Table()
    for info in info_to_retrieve:
        new_table[info] = data_table_data[info]

    # Duplicate raj and decj columns (in order to be corrected later)
    new_table['raj2000'] = \
        np.tile(data_table_data['raj2000'],
                (np.max(feeds) + 1, 1)).transpose()
    new_table['decj2000'] = \
        np.tile(data_table_data['decj2000'],
                (np.max(feeds) + 1, 1)).transpose()

    for f in list(set(feeds)):
        ra = new_table['raj2000'][:, f]
        dec = new_table['decj2000'][:, f]
        new_ra, new_dec = \
            correct_coordinates(ra, dec, new_table['derot_angle'],
                                xoffsets[f], yoffsets[f])

        new_table['raj2000'][:, f] = new_ra
        new_table['decj2000'][:, f] = new_dec

    for info in ['raj2000', 'decj2000', 'az', 'el', 'derot_angle']:
        new_table[info].unit = u.radian

    for ic, ch in enumerate(chan_ids):
        new_table['Ch{}'.format(ch)] = data_table_data['Ch{}'.format(ch)]

        new_table['Ch{}'.format(ch)].meta = {'polarization': polarizations[ic],
                                             'feed': feeds[ic],
                                             'IF': IFs[ic],
                                             'frequency': frequencies[ic],
                                             'bandwidth': bandwidths[ic],
                                             'xoffset': xoffsets[feeds[ic]],
                                             'yoffset': yoffsets[feeds[ic]],
                                             'relpower': relpowers[feeds[ic]],
                                             }
        new_table['Ch{}_feed'.format(ch)] = \
            np.zeros(len(data_table_data), dtype=np.uint8) + feeds[ic]

    lchdulist.close()
    return new_table


def read_data(fname):
    '''Read the data, whatever the format, and return them'''
    kind = detect_data_kind(fname)
    if kind == 'fitszilla':
        return read_data_fitszilla(fname)
    elif kind == 'hdf5':
        return Table.read(fname)


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
    for i in range(2):
        plt.plot(table.field('time'), table.field('Ch{}'.format(i))[:])
    plt.show()


def root_name(fname):
    return fname.replace('.fits', '').replace('.hdf5', '')
