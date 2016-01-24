# -*- coding: utf-8 -*-

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz
import os
from astropy.time import Time
import matplotlib.pyplot as plt
import unittest
from .io import print_obs_info_fitszilla, read_data, locations


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

        probe_location = SkyCoord(ra = new_table['ra'][:, 0],
                                  dec = new_table['dec'][:, 0],
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
