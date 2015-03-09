from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
from .io import read_data
import glob
from .read_config import read_config
import os
import numpy as np
from astropy import wcs
from astropy.table import Table, vstack
import logging
from scipy.optimize import curve_fit


def linear_fun(x, m, q):
    return m * x + q


def rough_baseline_sub(time, lc):
    '''Rough function to subtract the baseline'''
    m0 = 0
    q0 = min(lc)

    # only consider start and end quarters of image
    nbin = len(time)
    bins = np.arange(nbin, dtype=int)
    good = np.logical_or(bins <= nbin / 4, bins >= nbin / 4 * 3)
#    print(good, bins, bins[good], time[good], lc[good])
    par, pcov = curve_fit(linear_fun, time[good], lc[good], [m0, q0])

    return lc - linear_fun(time, *par)


class Scan():
    '''Class containing a single scan'''
    def __init__(self, fname, config_file=None, config=None):
        self.filename = fname
        assert config_file is not None or config is not None, \
            'Please specify either a config file or a config dictionary'
        self.config_file = config_file
        self.config = config
        if config is None:
            self.config = read_config(config_file)
        self.table = read_data(fname)
        self.chan_columns = [i for i in self.table.columns
                             if i.startswith('Ch')]

        self.nrows = len(self.table)

        self.baseline_subtract()
        self._update()

    def _update(self):
        '''Internally copy selected columns to attributes of the class.'''
        for col in self.table.columns:
            setattr(self, col, self.table.field(col))

    def update(self, name, value):
        '''Update the table and the attributes.'''

        if len(value) == self.nrows:
            self.table[name] = value
            self._update()
        else:
            self.__dict__[name] = value

    def baseline_subtract(self, kind='rough'):
        '''Subtract the baseline.'''
        if kind == 'rough':
            for col in self.chan_columns:
                self.table[col] = rough_baseline_sub(self.table['time'],
                                                     self.table[col])

    def zap_birdies(self):
        '''Zap bad intervals.'''
        pass

    def __repr__(self):
        '''Give the print() function something to print.'''
        reprstring = '\n\n----Scan from file {} ----\n'.format(self.filename)
        reprstring += repr(self.table)
        return reprstring

    def correct_coordinates(self):
        '''Correct coordinates for feed position.

        Uses the metadata in the channel columns xoffset and yoffset'''
        pass


class ScanSet():
    '''Class containing a set of scans'''
    def __init__(self, config_file=None, tablefile=None):
        self.config = read_config(config_file)
        self.config_file = config_file

        if tablefile is None:
            self.scans = self.load_scans()
            self.table = Table()
            self.update_table()
        else:
            self.table = Table()
            self.update_table(tablefile=tablefile)

    def __repr__(self):
        '''Give the print() function something to print.'''
        reprstring = ''
        for s in self.scans:
            reprstring += repr(s)

        return reprstring

    def load_scans(self):
        scans = []

        datadir = self.config['datadir']
        dirlist = self.config['list_of_directories']
        for d in dirlist:
            for f in glob.glob(os.path.join(datadir, d, '*')):
                scans.append(Scan(f, self.config_file, self.config))
        return scans

    def update_table(self, overwrite=False, tablefile=None):
        '''Update the table with all the subscans.'''

        if overwrite:
            self.table = Table()

        if tablefile is None:
            for s in self.scans:
                self.table = vstack((self.table, s.table))
        else:
            self.table = Table.read(tablefile, path='table')
            print(self.table)

        self.chan_columns = [i for i in self.table.columns
                             if i.startswith('Ch')]
        allras = self.table['raj2000']
        alldecs = self.table['decj2000']
        self.coordinates = np.array(list(zip(self.table['raj2000'],
                                             self.table['decj2000'])))
        self.mean_ra = np.mean(allras)
        self.mean_dec = np.mean(alldecs)
        self.min_ra = np.min(allras)
        self.min_dec = np.min(alldecs)
        self.max_ra = np.max(allras)
        self.max_dec = np.max(alldecs)

    def convert_coordinates(self):
        '''Convert the coordinates from sky to pixel.'''

        npix = np.array([int(n) for n in self.config['npix']])
        w = wcs.WCS(naxis=2)
        w.wcs.crpix = npix / 2
        delta_ra = self.max_ra - self.min_ra
        delta_dec = self.max_dec - self.min_dec
        w.wcs.cdelt = np.array([delta_ra / npix[0],
                                delta_dec / npix[1]])
        if self.config['reference_ra'] is None:
            ref_ra = self.mean_ra
        if self.config['reference_dec'] is None:
            ref_dec = self.mean_dec
        w.wcs.crval = np.array([ref_ra, ref_dec])
        w.wcs.ctype = ["RA---{}".format(self.config['projection']),
                       "DEC--{}".format(self.config['projection'])]

        pixcrd = w.wcs_world2pix(self.coordinates, 1)

        self.x = pixcrd[:, 0]
        self.y = pixcrd[:, 1]

    def save(self, filename=None):
        if filename is None:
            filename = 'bu.hd5'
        self.table.write(filename, path='table', overwrite=True)

    def baseline_subtract(self, kind='rough'):
        '''Subtract the baseline from all scans.'''
        for s in self.scans:
            s.baseline_subtract(kind)
        self.update_table(overwrite=True)

    def zap_birdies(self):
        '''Zap bad intervals.'''
        for s in self.scans:
            s.zap_birdies()
        self.update_table(overwrite=True)

    def calculate_images(self):
        '''Obtain image from all scans'''
        expomap, xedges, yedges = np.histogram2d(self.x, self.y,
                                                 bins=self.config['npix'])
        self.images = {}
        for ch in self.chan_columns:
            img, xedges, yedges = np.histogram2d(self.x, self.y,
                                                 bins=self.config['npix'],
                                                 weights=self.table[ch])
            self.images[ch] = img / expomap


def test_01_scan():
    '''Test that data are read.'''
    import os
    import matplotlib.pyplot as plt
    curdir = os.path.abspath(os.path.dirname(__file__))
    datadir = os.path.join(curdir, '..', '..', 'TEST_DATASET')

    fname = os.path.join(datadir, '20140603-103246-scicom-3C157',
                         '20140603-103246-scicom-3C157_003_003.fits')

    config = os.path.join(curdir, '..', '..', 'TEST_DATASET',
                          'test_config.ini')

    scan = Scan(fname, config)

    plt.ion()
    for col in scan.chan_columns:
        plt.plot(scan.time, getattr(scan, col))
    plt.draw()
    print('Drawn')


def test_02_scanset():
    '''Test that sets of data are read.'''
    import os
    curdir = os.path.abspath(os.path.dirname(__file__))
    config = os.path.join(curdir, '..', '..', 'TEST_DATASET',
                          'test_config.ini')

    scanset = ScanSet(config)

    scanset.save('test.hdf5')
#

def test_03_convert_coords():
    '''Test coordinate conversion.'''

    curdir = os.path.abspath(os.path.dirname(__file__))
    config = os.path.join(curdir, '..', '..', 'TEST_DATASET',
                          'test_config.ini')

    scanset = ScanSet(config, tablefile='test.hdf5')

    scanset.convert_coordinates()
    print(np.array(list(zip(scanset.x, scanset.y))))
    scanset.save('test_coord.hdf5')


def test_04_rough_image():
    '''Test coordinate conversion.'''

    curdir = os.path.abspath(os.path.dirname(__file__))
    config_file = os.path.join(curdir, '..', '..', 'TEST_DATASET',
                               'test_config.ini')

    config = read_config(config_file)
    scanset = ScanSet(config_file, tablefile='test_coord.hdf5')

    scanset.convert_coordinates()

    import matplotlib.pyplot as plt

    scanset.calculate_images()

    img = scanset.images['Ch0']

    plt.imshow(img)
    plt.ioff()
    plt.show()
