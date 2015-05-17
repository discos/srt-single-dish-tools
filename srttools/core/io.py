from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import astropy.io.fits as fits
from astropy.table import Table
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import os
from astropy.time import Time
import matplotlib.pyplot as plt
import unittest

DEBUG_MODE = False

locations={'SRT': EarthLocation(4865182.7660, 791922.6890, 4035137.1740,
                                unit=u.m),
           'Greenwich': EarthLocation(lat=51.477*u.deg, lon=0*u.deg)}


def standard_offsets():
    # Case with fake multibeam data in files. Damn it
    radius = 0.0006671036
    # 0 for feed 0, radius for the other six
    radii = np.array([0] + [radius]*6)
    # Feeds 1--6 are at angles -60, -120, etc. Here I use angle 0 for
    # convenience for feed 0, but it has no effect since radii[0] is 0
    feed_angles = -np.arange(0, 7, 1) * np.pi * 2/6

    xoffsets = radii * np.cos(feed_angles)
    yoffsets = radii * np.sin(feed_angles)
    return xoffsets, yoffsets


def detect_data_kind(fname):
    '''Placeholder for function that recognizes data format.'''
    if fname.endswith('.hdf5'):
        return 'hdf5'
    else:
        return 'fitszilla'


def correct_offsets(derot_angle, xoffset, yoffset):
    '''Correct feed offsets for derotation angle'''

    # Clockwise rotation of angle derot_angle
    new_xoff = xoffset * np.cos(derot_angle) - yoffset * np.sin(derot_angle)
    new_yoff = xoffset * np.sin(derot_angle) + yoffset * np.cos(derot_angle)

    return new_xoff, new_yoff


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


#@profile
def read_data_fitszilla(fname):
    '''Open a fitszilla FITS file and read all relevant information.'''
    global DEBUG_MODE
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

    if DEBUG_MODE and len(xoffsets) > 1:
        xoffsets, yoffsets = standard_offsets()

    relpowers = feed_input_data['relativePower']

    data_table_data = lchdulist['DATA TABLE'].data

    info_to_retrieve = ['time', 'derot_angle']

    new_table = Table()
    for info in info_to_retrieve:
        new_table[info] = data_table_data[info]

    if DEBUG_MODE:
        # Case with fake multibeam data in files. Still damn it
        new_table['derot_angle'][:] = 0

    # Duplicate raj and decj columns (in order to be corrected later)
    new_table['raj2000'] = \
        np.tile(data_table_data['raj2000'],
                (np.max(feeds) + 1, 1)).transpose()
    new_table['decj2000'] = \
        np.tile(data_table_data['decj2000'],
                (np.max(feeds) + 1, 1)).transpose()
    new_table['el'] = \
        np.tile(data_table_data['el'],
                (np.max(feeds) + 1, 1)).transpose()
    new_table['az'] = \
        np.tile(data_table_data['az'],
                (np.max(feeds) + 1, 1)).transpose()


    for info in ['raj2000', 'decj2000', 'az', 'el', 'derot_angle']:
        new_table[info].unit = u.radian

    # Coordinate correction. Will it work?
    for i in range(0, len(xoffsets)):
        # offsets < 0.001 arcseconds: don't correct (usually feed 0)
        if xoffsets[i] < np.radians(0.001 / 60.) and \
           yoffsets[i] < np.radians(0.001 / 60.):
               continue
        xoffs, yoffs = correct_offsets(new_table['derot_angle'],
                                       xoffsets[i],
                                       yoffsets[i])

        new_table['el'][:, i] += yoffs
        # TODO: Not sure about this cosine factor
        new_table['az'][:, i] += xoffs / np.cos(new_table['el'][:, i])

        obstimes = Time(new_table['time'] * u.day, format='mjd', scale='utc')
        coords = AltAz(az = new_table['az'][:, i],
                       alt = new_table['el'][:, i], unit= u.radian,
                       location=locations['SRT'],
                       obstime=obstimes)

        # According to line_profiler, coords.icrs is *by far* the longest
        # operation in this function, taking between 80 and 90% of the
        # execution time. Need to study a way to avoid this.
        coords_deg = coords.icrs
        new_table['raj2000'][:, i] = np.radians(coords_deg.ra)
        new_table['decj2000'][:, i] = np.radians(coords_deg.dec)

    for ic, ch in enumerate(chan_ids):
        new_table['Ch{}'.format(ch)] = \
            data_table_data['Ch{}'.format(ch)] * relpowers[feeds[ic]]

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

#@profile
def read_data(fname):
    '''Read the data, whatever the format, and return them'''
    kind = detect_data_kind(fname)
    if kind == 'fitszilla':
        return read_data_fitszilla(fname)
    elif kind == 'hdf5':
        return Table.read(fname)


#def test_open_data_fitszilla():
#    '''Test that data are read.'''
#    import os
#    import matplotlib.pyplot as plt
#    curdir = os.path.abspath(os.path.dirname(__file__))
#    datadir = os.path.join(curdir, '..', '..', 'TEST_DATASET')
#
#    fname = os.path.join(datadir, '20140603-103246-scicom-3C157',
#                         '20140603-103246-scicom-3C157_003_003.fits')
#    print_obs_info_fitszilla(fname)
#    table = read_data(fname)
#    for i in range(2):
#        plt.plot(table.field('time'), table.field('Ch{}'.format(i))[:])
#    plt.show()


class TestCoords(unittest.TestCase):
    @classmethod
    def setup_class(klass):
        global DEBUG_MODE
        DEBUG_MODE = True
        curdir = os.path.abspath(os.path.dirname(__file__))
        datadir = os.path.join(curdir, '..', '..', 'TEST_DATASET')

        fname = os.path.join(datadir, '20150410-001307-scicom-W44',
                             '20150410-001307-scicom-W44_002_003.fits')
        klass.table = read_data(fname)

    def step_coordinate_conversion(self):
        new_table = self.table

        probe_location = SkyCoord(ra = new_table['raj2000'][:, 0],
                                  dec = new_table['decj2000'][:, 0],
                                  unit= u.radian)
        print(new_table['time'])
        print((new_table['time'] * u.day).unit)
        obstimes = Time(new_table['time'] * u.day, format='mjd', scale='utc')

        print(obstimes)
        altaz = probe_location.transform_to(AltAz(location=locations['SRT'],
                                                  obstime=obstimes))
        print(altaz.alt, altaz.alt.unit)
        print(new_table['el'][:, 0].to(u.deg))
        delta_alt = (altaz.alt - new_table['el'][:, 0].to(u.deg))
        delta_az = (altaz.az - new_table['az'][:, 0].to(u.deg))

        print(delta_alt.to(u.arcsecond), delta_az.to(u.arcsecond))
        plt.hist2d(delta_alt.to(u.arcsecond).value,
                   delta_az.to(u.arcsecond).value, bins=20)
        plt.colorbar()
        plt.show()

    def test_coordinates(self):

        self.step_coordinate_conversion()


def root_name(fname):
    return fname.replace('.fits', '').replace('.hdf5', '')


#@profile
def profile_coords():
        '''Same test above, with profiling'''
        global DEBUG_MODE
        DEBUG_MODE = True
        curdir = os.path.abspath(os.path.dirname(__file__))
        datadir = os.path.join(curdir, '..', '..', 'TEST_DATASET')

        fname = os.path.join(datadir, '20150410-001307-scicom-W44',
                             '20150410-001307-scicom-W44_002_003.fits')
        new_table = read_data(fname)

#profile_coords()