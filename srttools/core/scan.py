from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
from .io import read_data, root_name
import glob
from .read_config import read_config, get_config_file
import os
import numpy as np
from astropy import wcs
from astropy.table import Table, vstack
import astropy.io.fits as fits
import logging
import matplotlib.pyplot as plt
from .fit import rough_baseline_sub, linear_fun


def list_scans(datadir, dirlist):
    '''List all scans contained in the directory listed in config'''
    scan_list = []

    for d in dirlist:
        for f in glob.glob(os.path.join(datadir, d, '*.fits')):
            scan_list.append(f)
    return scan_list


class Scan(Table):
    '''Class containing a single scan'''
    def __init__(self, data=None, config_file=None, norefilt=True,
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
            if os.path.exists(root_name(data) + '.hdf5'):
                data = root_name(data) + '.hdf5'
            print('Loading file {}'.format(data))
            table = read_data(data)
            Table.__init__(self, table, masked=True, **kwargs)
            self.meta['filename'] = os.path.abspath(data)
            self.meta['config_file'] = config_file

            self.meta.update(read_config(self.meta['config_file']))

            self.check_order()

            if 'ifilt' not in self.meta.keys() \
                    or not self.meta['ifilt'] \
                    or not norefilt:
                self.interactive_filter()
            if 'backsub' not in self.meta.keys() \
                    or not self.meta['backsub'] \
                    or not norefilt:
                print('Subtracting the baseline')
                self.baseline_subtract()
            self.save()

    def chan_columns(self):
        '''List columns containing samples'''
        return np.array([i for i in self.columns
                         if i.startswith('Ch') and not i.endswith('filt')])

    def baseline_subtract(self, kind='rough'):
        '''Subtract the baseline.'''
        if kind == 'rough':
            for col in self.chan_columns():
                self[col] = rough_baseline_sub(self['time'],
                                               self[col])
        self.meta['backsub'] = True

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
        '''Set default path and call Table.write'''
        t = Table(self)
        t.write(fname, path='scan', **kwargs)

    def check_order(self):
        '''Check that times in a scan are monotonically increasing'''
        assert np.all(self['time'] == np.sort(self['time'])), \
            'The order of times in the table is wrong'

    def interactive_filter(self, save=True):
        '''Run the interactive filter'''
        from .interactive_filter import select_data
        for ch in self.chan_columns():
            info = select_data(self['time'], self[ch])
            # Treat zapped intervals
            xs = info['zap'].xs
            good = np.ones(len(self['time']), dtype=bool)
            if len(xs) >= 2:
                intervals = list(zip(xs[:-1:2], xs[1::2]))
                for i in intervals:
                    good[np.logical_and(self['time'] >= i[0],
                                        self['time'] <= i[1])] = False
            self['{}-filt'.format(ch)] = good

            if len(info['fitpars']) > 1:
                self[ch] -= linear_fun(self['time'], *info['fitpars'])
            # TODO: make it channel-independent
                self.meta['backsub'] = True

            # TODO: make it channel-independent
            if info['FLAG']:
                self.meta['FLAG'] = True
        if save:
            self.save()
        self.meta['ifilt'] = True

    def save(self, fname=None):
        '''Call self.write with a default filename, or specify it.'''
        if fname is None:
            fname = root_name(self.meta['filename']) + '.hdf5'
        self.write(fname, overwrite=True)


class ScanSet(Table):
    '''Class containing a set of scans'''
    def __init__(self, data=None, norefilt=True, config_file=None,**kwargs):

        self.norefilt = norefilt
        if isinstance(data, Table):
            Table.__init__(self, data, **kwargs)
            if config_file is not None:
                config = read_config(config_file)
                self.meta.update(config)

            self.create_wcs()
        else:  # data is a config file
            config_file = data
            config = read_config(config_file)
            self.meta['scan_list'] = \
                self.list_scans(config['datadir'],
                                config['list_of_directories'])

            for i_s, s in enumerate(self.load_scans()):
                if 'FLAG' in s.meta.keys() and s.meta['FLAG']:
                    continue
                if i_s == 0:
                    scan_table = Table(s)
                else:
                    scan_table = vstack([scan_table, s],
                                        metadata_conflicts='silent')

            Table.__init__(self, scan_table)
            self.meta.update(config)
            self.meta['config_file'] = get_config_file()

            allras = self['raj2000']
            alldecs = self['decj2000']

            self.meta['mean_ra'] = np.mean(allras)
            self.meta['mean_dec'] = np.mean(alldecs)
            self.meta['min_ra'] = np.min(allras)
            self.meta['min_dec'] = np.min(alldecs)
            self.meta['max_ra'] = np.max(allras)
            self.meta['max_dec'] = np.max(alldecs)

            self.convert_coordinates()

        self.chan_columns = [i for i in self.columns
                             if i.startswith('Ch')]

    def list_scans(self, datadir, dirlist):
        '''List all scans contained in the directory listed in config'''
        return list_scans(datadir, dirlist)

    def load_scans(self):
        '''Load the scans in the list one by ones'''
        for f in self.meta['scan_list']:
            yield Scan(f, norefilt=self.norefilt)

    def get_coordinates(self):
        '''Give the coordinates as pairs of RA, DEC'''
        return np.array(list(zip(self['raj2000'],
                                 self['decj2000'])))

    def create_wcs(self):
        '''Create a wcs object from the pointing information'''
        npix = np.array(self.meta['npix'])
        self.wcs = wcs.WCS(naxis=2)

        self.wcs.wcs.crpix = npix / 2
        delta_ra = self.meta['max_ra'] - self.meta['min_ra']
        delta_dec = self.meta['max_dec'] - self.meta['min_dec']

        if not hasattr(self.meta, 'reference_ra'):
            self.meta['reference_ra'] = self.meta['mean_ra']
        if not hasattr(self.meta, 'reference_dec'):
            self.meta['reference_dec'] = self.meta['mean_dec']
        self.wcs.wcs.crval = np.array([self.meta['reference_ra'],
                                       self.meta['reference_dec']])
        self.wcs.wcs.cdelt = np.array([-delta_ra / npix[0],
                                       delta_dec / npix[1]])

        self.wcs.wcs.ctype = \
            ["RA---{}".format(self.meta['projection']),
             "DEC--{}".format(self.meta['projection'])]

    def convert_coordinates(self):
        '''Convert the coordinates from sky to pixel.'''
        self.create_wcs()

        pixcrd = self.wcs.wcs_world2pix(self.get_coordinates(), 1)

        self['x'] = pixcrd[:, 0]
        self['y'] = pixcrd[:, 1]

    def calculate_images(self):
        '''Obtain image from all scans'''
        images = {}
        for ch in self.chan_columns:
            if '{}-filt'.format(ch) in self.keys():
                good = self['{}-filt'.format(ch)]
            else:
                good = np.ones(len(self[ch]), dtype=bool)

            expomap, _, _ = np.histogram2d(self['x'][good], self['y'][good],
                                           bins=self.meta['npix'])

            img, _, _ = np.histogram2d(self['x'][good], self['y'][good],
                                       bins=self.meta['npix'],
                                       weights=self[ch][good])
            img_sq, _, _ = np.histogram2d(self['x'][good], self['y'][good],
                                          bins=self.meta['npix'],
                                          weights=self[ch][good] ** 2)
            mean = img / expomap
            images[ch] = mean
            images['{}-Sdev'.format(ch)] = img_sq / expomap - mean ** 2

        return images

    def write(self, fname, **kwargs):
        '''Set default path and call Table.write'''
        t = Table(self)
        t.write(fname, path='scanset', **kwargs)

    def save_ds9_images(self):
        '''Save a ds9-compatible file with one image per extension.'''
        images = self.calculate_images()
        self.create_wcs()

        hdulist = fits.HDUList()

        header = self.wcs.to_header()

        hdu = fits.PrimaryHDU(header=header)
        hdulist.append(hdu)
        for ic, ch in enumerate(self.chan_columns):
            hdu = fits.ImageHDU(images[ch], header=header, name='IMG' + ch)
            hdulist.append(hdu)
            hdu = fits.ImageHDU(images['{}-Sdev'.format(ch)], header=header,
                                name='{}-Sdev'.format(ch))
            hdulist.append(hdu)
        hdulist.writeto('img.fits', clobber=True)


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


def test_01b_read_scan():
    scan = Scan('scan.hdf5')
    plt.ion()
    for col in scan.chan_columns():
        plt.plot(scan['time'], scan[col])
    plt.show()

    return scan


def test_02_scanset():
    '''Test that sets of data are read.'''
    import os
    plt.ioff()
    curdir = os.path.abspath(os.path.dirname(__file__))
    config = os.path.join(curdir, '..', '..', 'TEST_DATASET',
                          'test_config.ini')

#    scanset = ScanSet(config, norefilt=False)
    scanset = ScanSet(config)

    scanset.write('test.hdf5', overwrite=True)


def test_03_rough_image():
    '''Test image production.'''

    plt.ion()
    curdir = os.path.abspath(os.path.dirname(__file__))
    config = os.path.join(curdir, '..', '..', 'TEST_DATASET',
                          'test_config.ini')
    scanset = ScanSet(Table.read('test.hdf5', path='scanset'),
                      config_file=config)

    images = scanset.calculate_images()

    img = images['Ch0']

    plt.figure('img')
    plt.imshow(img)
    plt.colorbar()
    plt.show()


def test_03_image_stdev():
    '''Test image production.'''

    curdir = os.path.abspath(os.path.dirname(__file__))
    config = os.path.join(curdir, '..', '..', 'TEST_DATASET',
                          'test_config.ini')
    scanset = ScanSet(Table.read('test.hdf5', path='scanset'),
                      config_file=config)

    images = scanset.calculate_images()

    img = images['Ch0-Sdev']

    plt.figure('log(img-Sdev)')
    plt.imshow(np.log10(img))
    plt.colorbar()
    plt.ioff()
    plt.show()


def test_04_ds9_image():
    '''Test image production.'''

    scanset = ScanSet.read('test.hdf5')

    scanset.save_ds9_images()
