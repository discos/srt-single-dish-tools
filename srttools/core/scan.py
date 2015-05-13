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
import astropy.units as u
import re


chan_re = re.compile(r'^Ch[0-9]+$')


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

#            if 'ifilt' not in self.meta.keys() \
#                    or not self.meta['ifilt'] \
#                    or not norefilt:
#                self.interactive_filter()
            if 'backsub' not in self.meta.keys() \
                    or not self.meta['backsub'] \
                    or not norefilt:
                print('Subtracting the baseline')
                self.baseline_subtract()

            self.save()

    def chan_columns(self):
        '''List columns containing samples'''
        return np.array([i for i in self.columns
                         if chan_re.match(i)])

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

    def write(self, fname, **kwargs):
        '''Set default path and call Table.write'''
        print('Saving to {}'.format(fname))
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
            # Temporary, waiting for AstroPy's metadata handling improvements
            feed = self[ch + '_feed'][0]

            selection = self['raj2000'][:, feed]

            ravar = np.abs(selection[-1] -
                           selection[0])

            selection = self['decj2000'][:, feed]
            decvar = np.abs(selection[-1] -
                            selection[0])

            # Choose if plotting by R.A. or Dec.
            if ravar > decvar:
                dim = 'raj2000'
            else:
                dim = 'decj2000'

            # ------- CALL INTERACTIVE FITTER ---------
            info = select_data(self[dim][:, feed], self[ch],
                               xlabel=dim)
            # -----------------------------------------

            # Treat zapped intervals
            xs = info['zap'].xs
            good = np.ones(len(self[dim]), dtype=bool)
            if len(xs) >= 2:
                intervals = list(zip(xs[:-1:2], xs[1::2]))
                for i in intervals:
                    good[np.logical_and(self[dim][:, feed] >= i[0],
                                        self[dim][:, feed] <= i[1])] = False
            self['{}-filt'.format(ch)] = good

            if len(info['fitpars']) > 1:
                self[ch] -= linear_fun(self[dim][:, feed], *info['fitpars'])
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
    def __init__(self, data=None, norefilt=True, config_file=None, **kwargs):

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
            scan_list = \
                self.list_scans(config['datadir'],
                                config['list_of_directories'])

            scan_list.sort()
            nscans = len(scan_list)

            for i_s, s in enumerate(self.load_scans(scan_list)):
                print('{}/{}'.format(i_s, nscans))
                if 'FLAG' in s.meta.keys() and s.meta['FLAG']:
                    continue
                s['Scan_id'] = i_s + np.zeros(len(s['time']), dtype=np.long)
                if i_s == 0:
                    scan_table = Table(s)
                else:
                    scan_table = vstack([scan_table, s],
                                        metadata_conflicts='silent')

            Table.__init__(self, scan_table)
            self.meta['scan_list'] = scan_list
            self.meta.update(config)
            self.meta['config_file'] = get_config_file()

            self.meta['scan_list'] = np.array(self.meta['scan_list'],
                                              dtype='S')

            allras = self['raj2000']
            alldecs = self['decj2000']

            self.meta['mean_ra'] = np.mean(allras)
            self.meta['mean_dec'] = np.mean(alldecs)
            self.meta['min_ra'] = np.min(allras)
            self.meta['min_dec'] = np.min(alldecs)
            self.meta['max_ra'] = np.max(allras)
            self.meta['max_dec'] = np.max(alldecs)

            self.convert_coordinates()

        self.chan_columns = np.array([i for i in self.columns
                                      if chan_re.match(i)])

    def list_scans(self, datadir, dirlist):
        '''List all scans contained in the directory listed in config'''
        return list_scans(datadir, dirlist)

    def load_scans(self, scan_list):
        '''Load the scans in the list one by ones'''
        for f in scan_list:
            yield Scan(f, norefilt=self.norefilt)

    def get_coordinates(self):
        '''Give the coordinates as pairs of RA, DEC'''

        return np.array(np.dstack([self['raj2000'],
                                   self['decj2000']]))

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

        # TODO: check consistency of units
        # Here I'm assuming all angles are radians
        crval = np.array([self.meta['reference_ra'],
                          self.meta['reference_dec']])
        self.wcs.wcs.crval = np.degrees(crval)

        cdelt = np.array([-delta_ra / npix[0],
                          delta_dec / npix[1]])
        self.wcs.wcs.cdelt = np.degrees(cdelt)

        self.wcs.wcs.ctype = \
            ["RA---{}".format(self.meta['projection']),
             "DEC--{}".format(self.meta['projection'])]

#    def scrunch_channels(self, feeds=None, polarizations=None,
#                         chan_names=None):
#        '''Scrunch channels and reduce their number.
#
#        POLARIZATIONS NOT IMPLEMENTED YET!
#        2-D lists of channels NOT IMPLEMENTED YET!
#
#        feed and polarization filters can be given as:
#
#        None:          all channels are to be summed in one
#        list of chans: channels in this list are summed, the others are deleted
#                       only one channel remains
#        2-d array:     the channels arr[0, :] will go to chan 0, arr[1, :] to
#                       chan 1, and so on.
#
#        At the end of the process, all channels have been eliminated but the
#        ones selected.
#        The axis-1 length of feeds and polarizations MUST be the same, unless
#        one of them is None.
#        '''
#        # TODO: Implement polarizations
#        # TODO: Implement 2-d arrays
#
#        allfeeds = np.array([self[ch + '_feed'][0]
#                             for ch in self.chan_columns])
#        if feeds is None:
#            feeds = list(set(allfeeds))
#
#        feed_mask = np.in1d(allfeeds, feeds)

    def convert_coordinates(self):
        '''Convert the coordinates from sky to pixel.'''
        self.create_wcs()

        self['x'] = np.zeros_like(self['raj2000'])
        self['y'] = np.zeros_like(self['decj2000'])
        coords = np.degrees(self.get_coordinates())
        for f in range(len(self['raj2000'][0, :])):
            pixcrd = self.wcs.wcs_world2pix(coords[:, f], 1)

            self['x'][:, f] = pixcrd[:, 0]
            self['y'][:, f] = pixcrd[:, 1]

    def calculate_images(self, scrunch=False):
        '''Obtain image from all scans'''
        images = {}
        xbins = np.linspace(np.min(self['x']),
                            np.max(self['x']),
                            self.meta['npix'][0])
        ybins = np.linspace(np.min(self['y']),
                            np.max(self['y']),
                            self.meta['npix'][1])

        total_expo = 0
        total_img = 0
        total_sdev = 0
        for ch in self.chan_columns:
            feeds = self[ch+'_feed']
            allfeeds = list(set(feeds))
            assert len(allfeeds) == 1, 'Feeds are mixed up in channels'
            feed = feeds[0]
            if '{}-filt'.format(ch) in self.keys():
                good = self['{}-filt'.format(ch)]
            else:
                good = np.ones(len(self[ch]), dtype=bool)

            expomap, _, _ = np.histogram2d(self['x'][:, feed][good],
                                           self['y'][:, feed][good],
                                           bins=[xbins, ybins])

            img, _, _ = np.histogram2d(self['x'][:, feed][good],
                                       self['y'][:, feed][good],
                                       bins=[xbins, ybins],
                                       weights=self[ch][good])
            img_sq, _, _ = np.histogram2d(self['x'][:, feed][good],
                                          self['y'][:, feed][good],
                                          bins=[xbins, ybins],
                                          weights=self[ch][good] ** 2)
            good = expomap > 0
            mean = img.copy()
            total_img += mean
            mean[good] /= expomap[good]
            images[ch] = mean
            img_sdev = img_sq
            total_sdev += img_sdev
            img_sdev[good] = img_sdev[good] / expomap[good] - mean[good] ** 2

            images['{}-Sdev'.format(ch)] = img_sdev
            total_expo += expomap

        self.images = images
        if scrunch:
            # Filter the part of the image whose value of exposure is higher
            # than the 10% percentile (avoid underexposed parts)
            good = total_expo > np.percentile(total_expo, 10)
            bad = np.logical_not(good)
            total_img[bad] = 0
            total_sdev[bad] = 0
            total_img[good] /= total_expo[good]
            total_sdev[good] = total_sdev[good] / total_expo[good] - total_img[good] ** 2

            images = {self.chan_columns[0]: total_img,
                      '{}-Sdev'.format(self.chan_columns[0]): total_sdev,
                      '{}-EXPO'.format(self.chan_columns[0]): total_expo}

        return images

    def interactive_display(self):
        '''Modify original scans from the image display'''
        from .interactive_filter import ImageSelector

        if not hasattr(self, 'images'):
            self.calculate_images()

        for ch in self.chan_columns:
            fig = plt.figure('Imageactive Display')
            ax = fig.add_subplot(111)
            img = self.images[ch]
            imagesel = ImageSelector(img, ax, fun=self.rerun_scan_analysis)
            plt.show()

    def rerun_scan_analysis(self, x, y, key):
        print(x, y, key)
        if key == 'a':
            good_entries = np.logical_and(self['x'].astype(int) == int(x),
                                          self['y'].astype(int) == int(y))
            scans = list(set(self['Scan_id'][good_entries]))
            for s in scans:
                sname = self.meta['scan_list'][s].decode()
                rescan = Scan(sname, norefilt=False)
                rescan.save()

        elif key == 'h':
            pass
        elif key == 'v':
            pass

    def write(self, fname, **kwargs):
        '''Set default path and call Table.write'''
        t = Table(self)
        t.write(fname, path='scanset', **kwargs)

    def save_ds9_images(self, fname=None, save_sdev=False, scrunch=False):
        '''Save a ds9-compatible file with one image per extension.'''
        if fname is None:
            fname = 'img.fits'
        images = self.calculate_images(scrunch=scrunch)
        self.create_wcs()

        hdulist = fits.HDUList()

        header = self.wcs.to_header()

        hdu = fits.PrimaryHDU(header=header)
        hdulist.append(hdu)
        
        for ic, ch in enumerate(images.keys()):
            is_sdev = ch.endswith('Sdev')
            
            if is_sdev and not save_sdev:
                continue
            
            hdu = fits.ImageHDU(images[ch], header=header, name='IMG' + ch)
            hdulist.append(hdu)

        hdulist.writeto(fname, clobber=True)


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

    scanset = ScanSet(config, norefilt=True)
#    scanset = ScanSet(config)

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
#    plt.ioff()
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


def test_03b_image_scrunch():
    '''Test image production.'''

    plt.ion()
    curdir = os.path.abspath(os.path.dirname(__file__))
    config = os.path.join(curdir, '..', '..', 'TEST_DATASET',
                          'test_config.ini')
    scanset = ScanSet(Table.read('test.hdf5', path='scanset'),
                      config_file=config)

    images = scanset.calculate_images(scrunch=True)

    img = images['Ch0']

    plt.figure('img - scrunched')
    plt.imshow(img)
    plt.colorbar()
    img = images['Ch0-Sdev']

    plt.figure('img - scrunched - sdev')
    plt.imshow(img)
    plt.colorbar()
    plt.ioff()
    plt.show()


def test_04_interactive_image():
    '''Test image production.'''

    curdir = os.path.abspath(os.path.dirname(__file__))
    config = os.path.join(curdir, '..', '..', 'TEST_DATASET',
                          'test_config.ini')
    scanset = ScanSet(Table.read('test.hdf5', path='scanset'),
                      config_file=config)

    scanset.interactive_display()


def test_04_ds9_image():
    '''Test image production.'''

    scanset = ScanSet.read('test.hdf5')

    scanset.save_ds9_images()


def test_05_all_multifeed():
    '''Test that sets of data are read also with multifeed.'''
    import os
    plt.ioff()
    curdir = os.path.abspath(os.path.dirname(__file__))
    config = os.path.join(curdir, '..', '..', 'TEST_DATASET',
                          'test_config_w44.ini')

    scanset = ScanSet(config, norefilt=True)
#    scanset = ScanSet(config)

    scanset.write('test_multifeed.hdf5', overwrite=True)

    scanset = ScanSet(Table.read('test_multifeed.hdf5', path='scanset'),
                      config_file=config)

    images = scanset.calculate_images()

    scanset.save_ds9_images('multifeed.fits')

    img = images['Ch0']

    plt.figure('img 0')
    plt.imshow(img)
    plt.colorbar()

    img = images['Ch0-Sdev']

    plt.figure('log(img 0 -Sdev)')
    plt.imshow(np.log10(img))
    plt.colorbar()

    img = images['Ch4']

    plt.figure('img 4')
    plt.imshow(img)
    plt.colorbar()

    img = images['Ch4-Sdev']

    plt.figure('log(img 4 -Sdev)')
    plt.imshow(np.log10(img))
    plt.colorbar()

    images = scanset.calculate_images(scrunch=True)
    scanset.save_ds9_images('scrunch.fits', scrunch=True)

    img = images['Ch0']

    plt.figure('img - scrunched')
    plt.imshow(img)
    plt.colorbar()
    img = images['Ch0-Sdev']

    plt.figure('log(img - scrunched - sdev)')
    plt.imshow(np.log10(img))
    plt.colorbar()
    plt.ioff()
    plt.show()
#
