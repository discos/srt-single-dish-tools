from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
from .io import read_data, root_name, DEBUG_MODE
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
from .interactive_filter import select_data
import astropy.units as u
import re
import unittest

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
                 interactive=False,
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

            if interactive:
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
        for ch in self.chan_columns():
            # Temporary, waiting for AstroPy's metadata handling improvements
            feed = self[ch + '_feed'][0]

            selection = self['ra'][:, feed]

            ravar = np.abs(selection[-1] -
                           selection[0])

            selection = self['dec'][:, feed]
            decvar = np.abs(selection[-1] -
                            selection[0])

            # Choose if plotting by R.A. or Dec.
            if ravar > decvar:
                dim = 'ra'
            else:
                dim = 'dec'

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

            tables = []

            for i_s, s in enumerate(self.load_scans(scan_list)):
                print('{}/{}'.format(i_s, nscans))
                if 'FLAG' in s.meta.keys() and s.meta['FLAG']:
                    continue
                s['Scan_id'] = i_s + np.zeros(len(s['time']), dtype=np.long)

                tables.append(s)
            scan_table = Table(vstack(tables))

            Table.__init__(self, scan_table)
            self.meta['scan_list'] = scan_list
            self.meta.update(config)
            self.meta['config_file'] = get_config_file()

            self.meta['scan_list'] = np.array(self.meta['scan_list'],
                                              dtype='S')
            self.analyze_coordinates(altaz=False)
            self.analyze_coordinates(altaz=True)

            self.convert_coordinates()

        self.chan_columns = np.array([i for i in self.columns
                                      if chan_re.match(i)])
        self.current = None

    def analyze_coordinates(self, altaz=False):
        if altaz:
            hor, ver = 'az', 'el'
        else:
            hor, ver = 'ra', 'dec'

        allhor = self[hor]
        allver = self[ver]

        self.meta['mean_' + hor] = np.mean(allhor)
        self.meta['mean_' + ver] = np.mean(allver)
        self.meta['min_' + hor] = np.min(allhor)
        self.meta['min_' + ver] = np.min(allver)
        self.meta['max_' + hor] = np.max(allhor)
        self.meta['max_' + ver] = np.max(allver)

    def list_scans(self, datadir, dirlist):
        '''List all scans contained in the directory listed in config'''
        return list_scans(datadir, dirlist)

    def load_scans(self, scan_list):
        '''Load the scans in the list one by ones'''
        for f in scan_list:
            yield Scan(f, norefilt=self.norefilt)

    def get_coordinates(self, altaz=False):
        '''Give the coordinates as pairs of RA, DEC'''

        if altaz:
            return np.array(np.dstack([self['az'],
                                       self['el']]))
        else:
            return np.array(np.dstack([self['ra'],
                                       self['dec']]))

    def create_wcs(self, altaz=False):
        '''Create a wcs object from the pointing information'''
        if altaz:
            hor, ver = 'az', 'el'
        else:
            hor, ver = 'ra', 'dec'
        npix = np.array(self.meta['npix'])
        self.wcs = wcs.WCS(naxis=2)

        self.wcs.wcs.crpix = npix / 2
        delta_hor = self.meta['max_' + hor] - self.meta['min_' + hor]
        delta_ver = self.meta['max_' + ver] - self.meta['min_' + ver]

        if not hasattr(self.meta, 'reference_' + hor):
            self.meta['reference_' + hor] = self.meta['mean_' + hor]
        if not hasattr(self.meta, 'reference_' + ver):
            self.meta['reference_' + ver] = self.meta['mean_' + ver]

        # TODO: check consistency of units
        # Here I'm assuming all angles are radians
        crval = np.array([self.meta['reference_' + hor],
                          self.meta['reference_' + ver]])
        self.wcs.wcs.crval = np.degrees(crval)

        cdelt = np.array([-delta_hor / npix[0],
                          delta_ver / npix[1]])
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

    def convert_coordinates(self, altaz=False):
        '''Convert the coordinates from sky to pixel.'''
        if altaz:
            hor, ver = 'az', 'el'
        else:
            hor, ver = 'ra', 'dec'
        self.create_wcs(altaz)

        self['x'] = np.zeros_like(self[hor])
        self['y'] = np.zeros_like(self[ver])
        coords = np.degrees(self.get_coordinates())
        for f in range(len(self[hor][0, :])):
            pixcrd = self.wcs.wcs_world2pix(coords[:, f], 0)

            self['x'][:, f] = pixcrd[:, 0]
            self['y'][:, f] = pixcrd[:, 1]

    def calculate_images(self, scrunch=False, no_offsets=False, altaz=False):
        '''Obtain image from all scans.

        scrunch:         sum all channels
        no_offsets:      use positions from feed 0 for all feeds'''
        images = {}
        xbins = np.linspace(np.min(self['x']),
                            np.max(self['x']),
                            self.meta['npix'][0] + 1)
        ybins = np.linspace(np.min(self['y']),
                            np.max(self['y']),
                            self.meta['npix'][1] + 1)

        total_expo = 0
        total_img = 0
        total_sdev = 0
        for ch in self.chan_columns:
            feeds = self[ch+'_feed']
            allfeeds = list(set(feeds))
            assert len(allfeeds) == 1, 'Feeds are mixed up in channels'
            if no_offsets:
                feed = 0
            else:
                feed = feeds[0]
            print(feed)
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
            total_img += mean.T
            mean[good] /= expomap[good]
            # For Numpy vs FITS image conventions...
            images[ch] = mean.T
            img_sdev = img_sq
            total_sdev += img_sdev.T
            img_sdev[good] = img_sdev[good] / expomap[good] - mean[good] ** 2

            images['{}-Sdev'.format(ch)] = img_sdev.T
            total_expo += expomap.T

        self.images = images
        if scrunch:
            # Filter the part of the image whose value of exposure is higher
            # than the 10% percentile (avoid underexposed parts)
            good = total_expo > np.percentile(total_expo, 10)
            bad = np.logical_not(good)
            total_img[bad] = 0
            total_sdev[bad] = 0
            total_img[good] /= total_expo[good]
            total_sdev[good] = total_sdev[good] / total_expo[good] - \
                total_img[good] ** 2

            images = {self.chan_columns[0]: total_img,
                      '{}-Sdev'.format(self.chan_columns[0]): total_sdev,
                      '{}-EXPO'.format(self.chan_columns[0]): total_expo}

        return images

    def interactive_display(self, ch = None, recreate=False):
        '''Modify original scans from the image display'''
        from .interactive_filter import ImageSelector

        if not hasattr(self, 'images') or recreate:
            self.calculate_images()

        if ch is None:
            chs = self.chan_columns
        else:
            chs = [ch]
        for ch in chs:
            fig = plt.figure('Imageactive Display')
            ax = fig.add_subplot(111)
            img = self.images['{}-Sdev'.format(ch)]
            self.current = ch
            imagesel = ImageSelector(img, ax, fun=self.rerun_scan_analysis)
            plt.show()

    def rerun_scan_analysis(self, x, y, key):
        print(x, y, key)
        if key == 'a':
            ra_xs = {}
            ra_ys = {}
            dec_xs = {}
            dec_ys = {}
            scan_ids = {}

            ch = self.current
            feed = list(set(self[ch+'_feed']))[0]

            # Select data inside the pixel +- 1
            good_entries = \
                np.logical_and(
                    np.abs(self['x'][:, feed].astype(int) - int(x)) <= 1,
                    np.abs(self['y'][:, feed].astype(int) - int(y)) <= 1)
            sids = list(set(self['Scan_id'][good_entries]))
            vars_to_filter = {}
            ra_masks = {}
            dec_masks = {}
            for sid in sids:
                sname = self.meta['scan_list'][sid].decode()
                s = Scan(sname)

                try:
                    chan_mask = s['{}-filt'.format(ch)]
                except:
                    chan_mask = np.zeros_like(s[ch])

                scan_ids[sname] = sid
                ras = s['ra'][:, feed]
                decs = s['dec'][:, feed]

                z = s[ch]

                ravar = np.max(ras) - np.min(ras)
                decvar = np.max(decs) - np.min(decs)
                if ravar > decvar:
                    vars_to_filter[sname] = 'ra'
                    ra_xs[sname] = ras
                    ra_ys[sname] = z
                    ra_masks[sname] = chan_mask
                else:
                    vars_to_filter[sname] = 'dec'
                    dec_xs[sname] = decs
                    dec_ys[sname] = z
                    dec_masks[sname] = chan_mask

            info = select_data(ra_xs, ra_ys, masks=ra_masks,
                               xlabel='RA', title='RA')
            plt.show()
            for sname in info.keys():
                mask = self['Scan_id'] == scan_ids[sname]
                s = Scan(sname)
                dim = vars_to_filter[sname]
                if len(info[sname]['zap'].xs) > 0:

                    xs = info[sname]['zap'].xs
                    good = np.ones(len(s[dim]), dtype=bool)
                    if len(xs) >= 2:
                        intervals = list(zip(xs[:-1:2], xs[1::2]))
                        for i in intervals:
                            good[np.logical_and(s[dim][:, feed] >= i[0],
                                                s[dim][:, feed] <= i[1])] = False
                    s['{}-filt'.format(ch)] = good
                    print(s['{}-filt'.format(ch)])
                    self['{}-filt'.format(ch)][mask] = good

                if len(info[sname]['fitpars']) > 1:
                    s[ch] -= linear_fun(s[dim][:, feed],
                                        *info[sname]['fitpars'])
                # TODO: make it channel-independent
                    s.meta['backsub'] = True
                    self[ch][mask][:] = s[ch]

                # TODO: make it channel-independent
                if info[sname]['FLAG']:
                    s.meta['FLAG'] = True
                    self['{}-filt'.format(ch)][mask] = np.zeros(len(s[dim]),
                                                                dtype=bool)

                s.save()

            info = select_data(dec_xs, dec_ys, masks=dec_masks,
                               xlabel='Dec', title='Dec')
            plt.show()

            for sname in info.keys():
                mask = self['Scan_id'] == scan_ids[sname]
                s = Scan(sname)
                dim = vars_to_filter[sname]
                if len(info[sname]['zap'].xs) > 0:

                    xs = info[sname]['zap'].xs
                    good = np.ones(len(s[dim]), dtype=bool)
                    if len(xs) >= 2:
                        intervals = list(zip(xs[:-1:2], xs[1::2]))
                        for i in intervals:
                            good[np.logical_and(s[dim][:, feed] >= i[0],
                                                s[dim][:, feed] <= i[1])] = False
                    s['{}-filt'.format(ch)] = good
                    print(s['{}-filt'.format(ch)])
                    self['{}-filt'.format(ch)][mask] = good

                if len(info[sname]['fitpars']) > 1:
                    s[ch] -= linear_fun(s[dim][:, feed],
                                        *info[sname]['fitpars'])
                # TODO: make it channel-independent
                    s.meta['backsub'] = True
                    self[ch][mask][:] = s[ch]

                # TODO: make it channel-independent
                if info[sname]['FLAG']:
                    s.meta['FLAG'] = True
                    self['{}-filt'.format(ch)][mask] = np.zeros(len(s[dim]),
                                                                dtype=bool)

                s.save()

            self.interactive_display(ch=ch, recreate=True)


        elif key == 'h':
            pass
        elif key == 'v':
            pass

    def write(self, fname, **kwargs):
        '''Set default path and call Table.write'''
        t = Table(self)
        t.write(fname, path='scanset', **kwargs)

    def save_ds9_images(self, fname=None, save_sdev=False, scrunch=False,
                        no_offsets=False, altaz=False):
        '''Save a ds9-compatible file with one image per extension.'''
        if fname is None:
            fname = 'img.fits'
        images = self.calculate_images(scrunch=scrunch, no_offsets=no_offsets,
                                       altaz=altaz)
        self.create_wcs(altaz)

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
