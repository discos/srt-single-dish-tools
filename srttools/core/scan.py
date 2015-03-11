from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
from .io import read_data
import glob
from .read_config import read_config, get_config_file
import os
import numpy as np
from astropy import wcs
from astropy.table import Table, vstack
import logging
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def linear_fun(x, m, q):
    '''A linear function'''
    return m * x + q


def rough_baseline_sub(time, lc):
    '''Rough function to subtract the baseline'''
    m0 = 0
    q0 = min(lc)

    # only consider start and end quarters of image
    nbin = len(time)
#    bins = np.arange(nbin, dtype=int)

    for percentage in [0.8, 0.15]:
        sorted_els = np.argsort(lc)
        # Select the lowest half elements
        good = sorted_els[: int(nbin * percentage)]
    #    good = np.logical_or(bins <= nbin / 4, bins >= nbin / 4 * 3)

        time_filt = time[good]
        lc_filt = lc[good]
        back_in_order = np.argsort(time_filt)
        lc_filt = lc_filt[back_in_order]
        time_filt = time_filt[back_in_order]
        par, pcov = curve_fit(linear_fun, time_filt, lc_filt, [m0, q0],
                              maxfev=6000)

        lc -= linear_fun(time, *par)

    return lc


class Scan(Table):
    '''Class containing a single scan'''
    def __init__(self, data=None, config_file=None,
                 **kwargs):

        if config_file is None:
            config_file = get_config_file()

        if isinstance(data, Table):
            Table.__init__(self, data, **kwargs)
        elif data is None:
            Table.__init__(self, **kwargs)
            self.meta['config_file'] = config_file
            self.meta.update(read_config(self.meta['config_file']))
        else:  # if data is a filename
            print('Loading file {}'.format(data))
            table = read_data(data)
            Table.__init__(self, table, **kwargs)
            self.meta['filename'] = data
            self.meta['config_file'] = config_file

            self.meta.update(read_config(self.meta['config_file']))

        self.baseline_subtract()

    def chan_columns(self):
        '''List columns containing samples'''
        return [i for i in self.columns
                if i.startswith('Ch')]

    def baseline_subtract(self, kind='rough'):
        '''Subtract the baseline.'''
        if kind == 'rough':
            for col in self.chan_columns():
                self[col] = rough_baseline_sub(self['time'],
                                               self[col])

    def zap_birdies(self):
        '''Zap bad intervals.'''
        pass

    def __repr__(self):
        '''Give the print() function something to print.'''
        reprstring = '\n\n----Scan from file {} ----\n'.format(self.filename)
        reprstring += repr(self)
        return reprstring

    def correct_coordinates(self):
        '''Correct coordinates for feed position.

        Uses the metadata in the channel columns xoffset and yoffset'''
        pass

    def write(self, fname, **kwargs):
        t = Table(self)
        t.write(fname, path='scan', **kwargs)


class ScanSet(Table):
    '''Class containing a set of scans'''
    def __init__(self, data=None, **kwargs):

        if isinstance(data, Table):
            Table.__init__(self, data, **kwargs)
        else:  # data is a config file
            config_file = data
            config = read_config(config_file)
            self.meta['scan_list'] = \
                self.list_scans(config['datadir'],
                                config['list_of_directories'])

            for i_s, s in enumerate(self.load_scans()):
                if i_s == 0:
                    scan_table = Table(s)
                else:
                    scan_table = vstack([scan_table, s],
                                        metadata_conflicts='silent')

            Table.__init__(self, scan_table)
            self.meta.update(config)
            self.meta['config_file'] = get_config_file()

            self.meta['chan_columns'] = [i for i in self.columns
                                         if i.startswith('Ch')]
            allras = self['raj2000']
            alldecs = self['decj2000']

            self.meta['mean_ra'] = np.mean(allras)
            self.meta['mean_dec'] = np.mean(alldecs)
            self.meta['min_ra'] = np.min(allras)
            self.meta['min_dec'] = np.min(alldecs)
            self.meta['max_ra'] = np.max(allras)
            self.meta['max_dec'] = np.max(alldecs)

            self.convert_coordinates()

    def list_scans(self, datadir, dirlist):
        '''List all scans contained in the directory listed in config'''
        scan_list = []

        for d in dirlist:
            for f in glob.glob(os.path.join(datadir, d, '*')):
                scan_list.append(f)
        return scan_list

    def load_scans(self):
        '''Load the scans in the list one by ones'''
        for f in self.meta['scan_list']:
            yield Scan(f)

    def get_coordinates(self):
        '''Give the coordinates as pairs of RA, DEC'''
        return np.array(list(zip(self['raj2000'],
                                 self['decj2000'])))

    def convert_coordinates(self):
        '''Convert the coordinates from sky to pixel.'''

        npix = np.array([int(n) for n in self.meta['npix']])
        w = wcs.WCS(naxis=2)
        w.wcs.crpix = npix / 2
        delta_ra = self.meta['max_ra'] - self.meta['min_ra']
        delta_dec = self.meta['max_dec'] - self.meta['min_dec']
        w.wcs.cdelt = np.array([delta_ra / npix[0],
                                delta_dec / npix[1]])
        if not hasattr(self.meta, 'reference_ra'):
            self.meta['reference_ra'] = self.meta['mean_ra']
        if not hasattr(self.meta, 'reference_dec'):
            self.meta['reference_dec'] = self.meta['mean_dec']

        w.wcs.crval = np.array([self.meta['reference_ra'],
                                self.meta['reference_dec']])

        w.wcs.ctype = ["RA---{}".format(self.meta['projection']),
                       "DEC--{}".format(self.meta['projection'])]

        pixcrd = w.wcs_world2pix(self.get_coordinates(), 1)

        self['x'] = pixcrd[:, 0]
        self['y'] = pixcrd[:, 1]

    def calculate_images(self):
        '''Obtain image from all scans'''
        expomap, xedges, yedges = np.histogram2d(self['x'], self['y'],
                                                 bins=self.meta['npix'])
        images = {}
        for ch in self.meta['chan_columns']:
            img, xedges, yedges = np.histogram2d(self['x'], self['y'],
                                                 bins=self.meta['npix'],
                                                 weights=self[ch])
            images[ch] = img / expomap

        return images

    def write(self, fname, **kwargs):
        t = Table(self)
        t.write(fname, path='scanset', **kwargs)


def test_01_scan():
    '''Test that data are read.'''
    import os
    curdir = os.path.dirname(__file__)
    datadir = os.path.join(curdir, '..', '..', 'TEST_DATASET')

    fname = \
        os.path.abspath(
            os.path.join(datadir, '20140603-103246-scicom-3C157',
                         '20140603-103246-scicom-3C157_003_003.fits'))

    config_file = \
        os.path.abspath(os.path.join(curdir, '..', '..', 'TEST_DATASET',
                                     'test_config.ini'))

    read_config(config_file)

    scan = Scan(fname)

    scan.write('scan.hdf5', overwrite=True)
    print('Drawn')


def test_01b_read_scan():
    scan = Scan('scan.hdf5')
    plt.ion()
    for col in scan.chan_columns():
        plt.plot(scan['time'], scan[col])
    plt.draw()

    return scan


def test_02_scanset():
    '''Test that sets of data are read.'''
    import os
    curdir = os.path.abspath(os.path.dirname(__file__))
    config = os.path.join(curdir, '..', '..', 'TEST_DATASET',
                          'test_config.ini')

    scanset = ScanSet(config)

    scanset.write('test.hdf5', overwrite=True)


def test_03_rough_image():
    '''Test image production.'''

    scanset = ScanSet(Table.read('test.hdf5', path='scanset'))

    import matplotlib.pyplot as plt

    images = scanset.calculate_images()

    img = images['Ch0']

    plt.figure('img')
    plt.imshow(img)
    plt.colorbar()
    plt.ioff()
    plt.show()
