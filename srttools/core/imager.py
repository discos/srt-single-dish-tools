"""Produce calibrated images.

``SDTimage`` is a script that, given a list of cross scans composing an
on-the-fly map, is able to calculate the map and save it in FITS format after
cleaning the data.
"""
from __future__ import (absolute_import, division,
                        print_function)

from .scan import Scan, chan_re, list_scans
from .read_config import read_config, sample_config_file
import numpy as np
from astropy import wcs
from astropy.table import Table, vstack, Column
import astropy.io.fits as fits
import astropy.units as u
try:
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

from .fit import linear_fun
from .interactive_filter import select_data
from .calibration import CalibratorTable

import sys
import warnings
import logging
import traceback
from .global_fit import fit_full_image
import six
import functools


__all__ = ["ScanSet"]


def _load_calibration(calibration, map_unit):
    caltable = CalibratorTable().read(calibration, path='table')
    caltable.update()
    caltable.compute_conversion_function(map_unit)

    if map_unit == "Jy/beam":
        conversion_units = u.Jy / u.ct
    elif map_unit in ["Jy/pixel", "Jy/sr"]:
        conversion_units = u.Jy / u.ct / u.steradian
    else:
        raise ValueError("Unit for calibration not recognized")
    return caltable, conversion_units


class ScanSet(Table):
    def __init__(self, data=None, norefilt=True, config_file=None,
                 freqsplat=None, nofilt=False, nosub=False, **kwargs):
        """Class obtained by a set of scans.

        Once the scans are loaded, this class contains all functionality that
        will be used to produce (calibrated or uncalibrated) maps with WCS
        information.

        Parameters
        ----------
        data : str or None
            data can be one of the following:
            + a config file, containing the information on the scans to load
            + an HDF5 archive, containing a former scanset
            + another ScanSet or an Astropy Table
        config_file : str
            Config file containing the parameters for the images and the
            directories containing the image and calibration data
        norefilt : bool
            See :class:`srttools.core.scan.Scan`
        freqsplat : str
            See :class:`srttools.core.scan.interpret_frequency_range`
        nofilt : bool
            See :class:`srttools.core.scan.clean_scan_using_variability`
        nosub : bool
            See :class:`srttools.core.scan.Scan`

        Other Parameters
        ----------------
        kwargs : additional arguments
            These will be passed to Scan initializers

        Examples
        --------
        >>> scanset = ScanSet()  # An empty scanset
        >>> isinstance(scanset, ScanSet)
        True
        """
        if data is None and config_file is None:
            Table.__init__(self, data, **kwargs)
            return
        self.norefilt = norefilt
        self.freqsplat = freqsplat

        if isinstance(data, six.string_types) and data.endswith('hdf5'):
            data = Table.read(data, path='scanset')

            txtfile = data.meta['scan_list_file']

            with open(txtfile, 'r') as fobj:
                self.scan_list = []
                for i in fobj.readlines():
                    self.scan_list.append(i.strip())

        if isinstance(data, Table):
            Table.__init__(self, data, **kwargs)
            if config_file is not None:
                config = read_config(config_file)
                self.meta.update(config)

            self.create_wcs()

        else:  # data is a config file
            config_file = data
            config = read_config(config_file)
            self.meta.update(config)
            self.meta['config_file'] = config_file

            scan_list = \
                self.list_scans()

            scan_list.sort()

            tables = []

            for i_s, s in self.load_scans(scan_list,
                                          freqsplat=freqsplat, nofilt=nofilt,
                                          nosub=nosub, **kwargs):

                if 'FLAG' in s.meta.keys() and s.meta['FLAG']:
                    continue
                s['Scan_id'] = i_s + np.zeros(len(s['time']), dtype=np.long)

                del s.meta['filename']
                del s.meta['calibrator_directories']
                del s.meta['list_of_directories']
                tables.append(s)

            scan_table = Table(vstack(tables))

            Table.__init__(self, scan_table)
            self.scan_list = scan_list

            self.meta['scan_list_file'] = None

            self.analyze_coordinates(altaz=False)

            self.convert_coordinates()

        self.chan_columns = np.array([i for i in self.columns
                                      if chan_re.match(i)])
        self.current = None

    def analyze_coordinates(self, altaz=False):
        """Save statistical information on coordinates."""
        if altaz:
            hor, ver = 'delta_az', 'delta_el'
        else:
            hor, ver = 'ra', 'dec'

        if 'delta_az' not in self.columns and altaz:
            self.calculate_delta_altaz()

        allhor = self[hor]
        allver = self[ver]
        hor_unit = self[hor].unit
        ver_unit = self[ver].unit

        # These seemingly useless float() calls are needed for serialize_meta
        self.meta['mean_' + hor] = float(np.mean(allhor)) * hor_unit
        self.meta['mean_' + ver] = float(np.mean(allver)) * ver_unit
        self.meta['min_' + hor] = float(np.min(allhor)) * hor_unit
        self.meta['min_' + ver] = float(np.min(allver)) * ver_unit
        self.meta['max_' + hor] = float(np.max(allhor)) * hor_unit
        self.meta['max_' + ver] = float(np.max(allver)) * ver_unit

        if 'reference_ra' not in self.meta:
            self.meta['reference_ra'] = self.meta['RA']
        if 'reference_dec' not in self.meta:
            self.meta['reference_dec'] = self.meta['Dec']

    def list_scans(self, datadir=None, dirlist=None):
        """List all scans contained in the directory listed in config."""
        if datadir is None:
            datadir = self.meta['datadir']
            dirlist = self.meta['list_of_directories']
        return list_scans(datadir, dirlist)

    def load_scans(self, scan_list, freqsplat=None, nofilt=False, **kwargs):
        """Load the scans in the list one by ones."""
        nscan = len(scan_list)
        for i, f in enumerate(scan_list):
            print("{}/{}".format(i + 1, nscan), end="\r")
            try:
                s = Scan(f, norefilt=self.norefilt, freqsplat=freqsplat,
                         nofilt=nofilt, **kwargs)
                yield i, s
            except KeyError as e:
                warnings.warn(
                    "Error while processing {}: Missing key: {}".format(f,
                                                                        str(e))
                )
            except Exception as e:
                warnings.warn(traceback.format_exc())
                warnings.warn("Error while processing {}: {}".format(f,
                                                                     str(e)))

    def get_coordinates(self, altaz=False):
        """Give the coordinates as pairs of RA, DEC."""
        if altaz:
            return np.array(np.dstack([self['delta_az'],
                                       self['delta_el']]))
        else:
            return np.array(np.dstack([self['ra'],
                                       self['dec']]))

    def get_obstimes(self):
        """Get `astropy.Time` object for time at the telescope location."""
        from astropy.time import Time
        from .io import locations
        return Time((self['time']) * u.day, format='mjd', scale='utc',
                    location=locations[self.meta['site']])

    def apply_user_filter(self, user_func=None, out_column=None):
        """Apply a user-supplied function as filter.

        Parameters
        ----------
        user_func : function
            This function needs to accept a `scanset` as only argument.
            `ScanSet` object. It has to return an array with the same length of
            a column of `scanset`
        out_column : str
            column where the results will be stored

        Returns
        -------
        retval : array
            the result of user_func
        """
        if user_func is None:
            raise ValueError('user_func needs to be specified')
        retval = user_func(self)
        if out_column is not None:
            self[out_column] = retval
        return retval

    def calculate_delta_altaz(self):
        """Construction of delta altaz coordinates.

        Calculate the delta of altazimutal coordinates wrt the position
        of the source
        """
        from astropy.coordinates import SkyCoord

        from .io import locations

        ref_coords = SkyCoord(ra=self.meta['reference_ra'],
                              dec=self.meta['reference_dec'],
                              obstime=self.get_obstimes(),
                              location=locations[self.meta['site']])
        ref_altaz_coords = ref_coords.altaz
        ref_az = ref_altaz_coords.az.to(u.rad)
        ref_el = ref_altaz_coords.alt.to(u.rad)

        self.meta['reference_delta_az'] = 0*u.rad
        self.meta['reference_delta_el'] = 0*u.rad
        self['delta_az'] = np.zeros_like(self['az'])
        self['delta_el'] = np.zeros_like(self['el'])
        for f in range(len(self['el'][0, :])):
            self['delta_az'][:, f] = \
                (self['az'][:, f] - ref_az) * np.cos(ref_el)
            self['delta_el'][:, f] = self['el'][:, f] - ref_el

        if HAS_MPL:
            fig1 = plt.figure("adsfasdfasd")
            plt.plot(self['delta_az'], self['delta_el'])
            plt.savefig('delta_altaz.png')
            plt.close(fig1)

            fig2 = plt.figure("adsfasdf")
            plt.plot(self['az'], self['el'])
            plt.plot(ref_az, ref_el)
            plt.savefig('altaz_with_src.png')
            plt.close(fig2)

    def create_wcs(self, altaz=False):
        """Create a wcs object from the pointing information."""
        if altaz:
            hor, ver = 'delta_az', 'delta_el'
        else:
            hor, ver = 'ra', 'dec'
        pixel_size = self.meta['pixel_size']
        self.wcs = wcs.WCS(naxis=2)

        if 'max_' + hor not in self.meta:
            self.analyze_coordinates(altaz)

        delta_hor = self.meta['max_' + hor] - self.meta['min_' + hor]
        delta_ver = self.meta['max_' + ver] - self.meta['min_' + ver]

        npix_hor = np.ceil(delta_hor / pixel_size)
        npix_ver = np.ceil(delta_ver / pixel_size)

        self.meta['npix'] = np.array([npix_hor, npix_ver])

        self.wcs.wcs.crpix = self.meta['npix'] / 2

        # TODO: check consistency of units
        # Here I'm assuming all angles are radians
        crval = np.array([self.meta['reference_' + hor].to(u.rad).value,
                          self.meta['reference_' + ver].to(u.rad).value])

        self.wcs.wcs.crval = np.degrees(crval)

        cdelt = np.array([-pixel_size.to(u.rad).value,
                          pixel_size.to(u.rad).value])
        self.wcs.wcs.cdelt = np.degrees(cdelt)

        self.wcs.wcs.ctype = \
            ["RA---{}".format(self.meta['projection']),
             "DEC--{}".format(self.meta['projection'])]

    def convert_coordinates(self, altaz=False):
        """Convert the coordinates from sky to pixel."""
        if altaz:
            hor, ver = 'delta_az', 'delta_el'
        else:
            hor, ver = 'ra', 'dec'
        self.create_wcs(altaz)

        self['x'] = np.zeros_like(self[hor])
        self['y'] = np.zeros_like(self[ver])
        coords = np.degrees(self.get_coordinates(altaz=altaz))
        for f in range(len(self[hor][0, :])):
            pixcrd = self.wcs.wcs_world2pix(coords[:, f], 0)

            self['x'][:, f] = pixcrd[:, 0]
            self['y'][:, f] = pixcrd[:, 1]
        self['x'].meta['altaz'] = altaz
        self['y'].meta['altaz'] = altaz

    def calculate_images(self, scrunch=False, no_offsets=False, altaz=False,
                         calibration=None, elevation=None, map_unit="Jy/beam",
                         calibrate_scans=False):
        """Obtain image from all scans.

        scrunch:         sum all channels
        no_offsets:      use positions from feed 0 for all feeds.
        """
        if altaz != self['x'].meta['altaz']:
            self.convert_coordinates(altaz)

        images = {}

        xbins = np.linspace(0,
                            self.meta['npix'][0],
                            self.meta['npix'][0] + 1)
        ybins = np.linspace(0,
                            self.meta['npix'][1],
                            self.meta['npix'][1] + 1)

        total_expo = 0
        total_img = 0
        total_sdev = 0
        for ch in self.chan_columns:
            feeds = self[ch+'_feed']
            allfeeds = list(set(feeds))
            if not len(allfeeds) == 1:
                raise ValueError('Feeds are mixed up in channels')
            if no_offsets:
                feed = 0
            else:
                feed = feeds[0]

            if elevation is None:
                elevation = np.mean(self['el'][:, feed])

            if '{}-filt'.format(ch) in self.keys():
                good = self['{}-filt'.format(ch)]
            else:
                good = np.ones(len(self[ch]), dtype=bool)

            expomap, _, _ = np.histogram2d(self['x'][:, feed][good],
                                           self['y'][:, feed][good],
                                           bins=[xbins, ybins])

            counts = np.array(self[ch][good])

            if calibration is not None and calibrate_scans:
                caltable, conversion_units = _load_calibration(calibration,
                                                               map_unit)
                area_conversion, final_unit = \
                    self._calculate_calibration_factors(map_unit)

                Jy_over_counts, Jy_over_counts_err = conversion_units * \
                    caltable.Jy_over_counts(channel=ch, map_unit=map_unit,
                                            elevation=self['el'][:, feed][good]
                                            )

                counts = counts * u.ct * area_conversion * Jy_over_counts
                counts = counts.to(final_unit).value

            img, _, _ = np.histogram2d(self['x'][:, feed][good],
                                       self['y'][:, feed][good],
                                       bins=[xbins, ybins],
                                       weights=counts)
            img_sq, _, _ = np.histogram2d(self['x'][:, feed][good],
                                          self['y'][:, feed][good],
                                          bins=[xbins, ybins],
                                          weights=counts ** 2)

            good = expomap > 0
            mean = img.copy()
            total_img += mean.T
            mean[good] /= expomap[good]
            # For Numpy vs FITS image conventions...
            images[ch] = mean.T
            img_sdev = img_sq
            total_sdev += img_sdev.T
            img_sdev[good] = img_sdev[good] / expomap[good] - mean[good] ** 2

            img_sdev = np.sqrt(img_sdev)
            if calibration is not None and calibrate_scans:
                cal_rel_err = \
                    np.mean(Jy_over_counts_err / Jy_over_counts).value
                img_sdev += mean * cal_rel_err

            images['{}-Sdev'.format(ch)] = img_sdev.T
            total_expo += expomap.T

        self.images = images
        if calibration is not None and not calibrate_scans:
            self.calibrate_images(calibration, elevation=elevation,
                                  map_unit=map_unit)

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
                      '{}-Sdev'.format(self.chan_columns[0]): np.sqrt(
                          total_sdev),
                      '{}-EXPO'.format(self.chan_columns[0]): total_expo}

        return images

    def fit_full_images(self, chans=None, fname=None, save_sdev=False,
                        scrunch=False, no_offsets=False, altaz=False,
                        calibration=None, excluded=None, par=None,
                        map_unit="Jy/beam"):
        """Flatten the baseline with a global fit.

        Fit a linear trend to each scan to minimize the scatter in an image
        """

        if not hasattr(self, 'images'):
            self.calculate_images(scrunch=scrunch, no_offsets=no_offsets,
                                  altaz=altaz, calibration=calibration,
                                  map_unit=map_unit)

        if chans is not None:
            chans = chans.split(',')
        else:
            chans = self.chan_columns

        for ch in chans:
            print("Fitting channel {}".format(ch))
            feeds = self[ch + '_feed']
            allfeeds = list(set(feeds))
            if not len(allfeeds) == 1:
                raise ValueError('Feeds are mixed up in channels')
            if no_offsets:
                feed = 0
            else:
                feed = feeds[0]
            self[ch + "_save"] = self[ch].copy()
            self[ch] = Column(fit_full_image(self, chan=ch, feed=feed,
                                             excluded=excluded, par=par))

        self.calculate_images(scrunch=scrunch, no_offsets=no_offsets,
                              altaz=altaz, calibration=calibration,
                              map_unit=map_unit)

    def _calculate_calibration_factors(self, map_unit):
        if map_unit == "Jy/beam":
            area_conversion = 1
            final_unit = u.Jy
        elif map_unit == "Jy/sr":
            area_conversion = 1
            final_unit = u.Jy / u.sr
        elif map_unit == "Jy/pixel":
            area_conversion = self.meta['pixel_size'] ** 2
            final_unit = u.Jy
        return area_conversion, final_unit

    def calibrate_images(self, calibration, elevation=np.pi/4,
                         map_unit="Jy/beam"):
        """Calibrate the images."""
        if not hasattr(self, 'images'):
            self.calculate_images()

        caltable, conversion_units = _load_calibration(calibration, map_unit)

        for ch in self.chan_columns:
            Jy_over_counts, Jy_over_counts_err = \
                caltable.Jy_over_counts(channel=ch, map_unit=map_unit,
                                        elevation=elevation) * \
                conversion_units

            if np.isnan(Jy_over_counts):
                warnings.warn("The Jy/counts factor is nan")
                continue
            A = self.images[ch].copy() * u.ct
            eA = self.images['{}-Sdev'.format(ch)].copy() * u.ct

            self.images['{}-RAW'.format(ch)] = \
                self.images['{}'.format(ch)].copy()
            self.images['{}-Sdev-RAW'.format(ch)] = \
                self.images['{}-Sdev'.format(ch)].copy()
            bad = eA != eA
            A[bad] = 1 * u.ct
            eA[bad] = 0 * u.ct

            bad = np.logical_or(A == 0, A != A)
            A[bad] = 1 * u.ct
            eA[bad] = 0 * u.ct

            B = Jy_over_counts
            eB = Jy_over_counts_err

            area_conversion, final_unit = \
                self._calculate_calibration_factors(map_unit)

            C = A * area_conversion * Jy_over_counts
            C[bad] = 0

            self.images[ch] = C.to(final_unit).value

            eC = C * (eA / A + eB / B)

            self.images['{}-Sdev'.format(ch)] = eC.to(final_unit).value

    def interactive_display(self, ch=None, recreate=False, test=False):
        """Modify original scans from the image display."""
        from .interactive_filter import ImageSelector
        if not HAS_MPL:
            raise ImportError('interactive_display: '
                              'matplotlib is not installed')

        if not hasattr(self, 'images') or recreate:
            self.calculate_images()

        if ch is None:
            chs = self.chan_columns
        else:
            chs = [ch]
        if test:
            chs = ['Ch0']
        for ch in chs:
            fig = plt.figure('Imageactive Display')
            gs = GridSpec(1, 2, width_ratios=(3, 2))
            ax = fig.add_subplot(gs[0])
            ax2 = fig.add_subplot(gs[1])
            imgch = ch
            sdevch = '{}-Sdev'.format(ch)
            if '{}-RAW'.format(ch) in self.images.keys():
                imgch = '{}-RAW'.format(ch)
                sdevch = '{}-Sdev-RAW'.format(ch)
            img = self.images[imgch]
            ax2.imshow(img, origin='lower',
                       vmin=np.percentile(img, 20), cmap="gnuplot2",
                       interpolation="nearest")

            img = self.images[sdevch].copy()
            self.current = ch
            bad = np.logical_or(img == 0, img != img)
            img[bad] = np.mean(img[np.logical_not(bad)])
            fun = functools.partial(self.rerun_scan_analysis, test=test)
            imgsel = ImageSelector(img, ax, fun=fun,
                                   test=test)
        return imgsel

    def rerun_scan_analysis(self, x, y, key, test=False):
        """Rerun the analysis of single scans."""
        logging.debug("{} {} {}".format(x, y, key))
        if key == 'a':
            self.reprocess_scans_through_pixel(x, y, test=test)
        elif key == 'h':
            pass
        elif key == 'v':
            pass

    def reprocess_scans_through_pixel(self, x, y, test=False):
        """Given a pixel in the image, find all scans passing through it."""
        ch = self.current

        ra_xs, ra_ys, dec_xs, dec_ys, scan_ids, ra_masks, dec_masks, \
            vars_to_filter = \
            self.find_scans_through_pixel(x, y, test=test)

        info = select_data(ra_xs, ra_ys, masks=ra_masks,
                           xlabel="RA", title="RA", test=test)

        for sname in info.keys():
            self.update_scan(sname, scan_ids[sname], vars_to_filter[sname],
                             info[sname]['zap'],
                             info[sname]['fitpars'], info[sname]['FLAG'])

        info = select_data(dec_xs, dec_ys, masks=dec_masks, xlabel="Dec",
                           title="Dec", test=test)

        for sname in info.keys():
            self.update_scan(sname, scan_ids[sname], vars_to_filter[sname],
                             info[sname]['zap'],
                             info[sname]['fitpars'], info[sname]['FLAG'])

        display = self.interactive_display(ch=ch, recreate=True, test=test)
        return display

    def find_scans_through_pixel(self, x, y, test=False):
        """Find scans passing through a pixel."""
        ra_xs = {}
        ra_ys = {}
        dec_xs = {}
        dec_ys = {}
        scan_ids = {}
        ra_masks = {}
        dec_masks = {}
        vars_to_filter = {}

        if not test:
            ch = self.current
        else:
            ch = 'Ch0'

        feed = list(set(self[ch+'_feed']))[0]

        # Select data inside the pixel +- 1

        good_entries = \
            np.logical_and(
                np.abs(self['x'][:, feed] - x) < 1,
                np.abs(self['y'][:, feed] - y) < 1)

        sids = list(set(self['Scan_id'][good_entries]))

        for sid in sids:
            sname = self.scan_list[sid]
            try:
                s = Scan(sname)
            except Exception:
                warnings.warn("Errors while opening scan {}".format(sname))
                continue
            try:
                chan_mask = s['{}-filt'.format(ch)]
            except Exception:
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

        return ra_xs, ra_ys, dec_xs, dec_ys, scan_ids, ra_masks, dec_masks, \
            vars_to_filter

    def update_scan(self, sname, sid, dim, zap_info, fit_info, flag_info,
                    test=False):
        """Update a scan in the scanset after filtering."""
        ch = self.current
        if test:
            ch = 'Ch0'
        feed = list(set(self[ch+'_feed']))[0]
        mask = self['Scan_id'] == sid
        try:
            s = Scan(sname)
        except Exception:
            return

        if len(zap_info.xs) > 0:

            xs = zap_info.xs
            good = np.ones(len(s[dim]), dtype=bool)
            if len(xs) >= 2:
                intervals = list(zip(xs[:-1:2], xs[1::2]))
                for i in intervals:
                    good[np.logical_and(s[dim][:, feed] >= i[0],
                                        s[dim][:, feed] <= i[1])] = False
            s['{}-filt'.format(ch)] = good
            self['{}-filt'.format(ch)][mask] = good

        if len(fit_info) > 1:
            s[ch] -= linear_fun(s[dim][:, feed],
                                *fit_info)
        # TODO: make it channel-independent
            s.meta['backsub'] = True
            self[ch][mask] = s[ch]

        # TODO: make it channel-independent
        if flag_info:
            s.meta['FLAG'] = True
            self['{}-filt'.format(ch)][mask] = np.zeros(len(s[dim]),
                                                        dtype=bool)
            s['{}-filt'.format(ch)] = np.zeros(len(s[dim]), dtype=bool)

        s.save()

    def barycenter_times(self):
        """Create barytime column with observing times converted to TDB."""
        obstimes_tdb = self.get_obstimes().tdb.mjd
        self['barytime'] = obstimes_tdb
        return obstimes_tdb

    def write(self, fname, **kwargs):
        """Same as Table.write, but adds path information for HDF5.

        Moreover, saves the scan list to a txt file, that will be read when
        data are reloaded. This is a *temporary solution*
        """
        import os
        f, _ = os.path.splitext(fname)
        txtfile = f + '_scan_list.txt'
        self.meta['scan_list_file'] = txtfile
        with open(txtfile, 'w') as fobj:
            for i in self.scan_list:
                print(i, file=fobj)

        Table.write(self, fname, path='scanset', serialize_meta=True, **kwargs)

    def save_ds9_images(self, fname=None, save_sdev=False, scrunch=False,
                        no_offsets=False, altaz=False, calibration=None,
                        map_unit="Jy/beam", calibrate_scans=False):
        """Save a ds9-compatible file with one image per extension."""
        if fname is None:
            tail = '.fits'
            if altaz:
                tail = '_altaz.fits'
            fname = self.meta['config_file'].replace('.ini', tail)
        images = self.calculate_images(scrunch=scrunch, no_offsets=no_offsets,
                                       altaz=altaz, calibration=calibration,
                                       map_unit=map_unit,
                                       calibrate_scans=calibrate_scans)

        self.create_wcs(altaz)

        hdulist = fits.HDUList()

        header = self.wcs.to_header()
        if map_unit == "Jy/beam" and calibration is not None:
            caltable = CalibratorTable.read(calibration)
            beam, _ = caltable.beam_width()
            std_to_fwhm = np.sqrt(8 * np.log(2))
            header['bmaj'] = np.degrees(beam) * std_to_fwhm
            header['bmin'] = np.degrees(beam) * std_to_fwhm
            header['bpa'] = 0

        if calibration is not None:
            header['bunit'] = map_unit

        hdu = fits.PrimaryHDU(header=header)
        hdulist.append(hdu)

        keys = list(images.keys())
        keys.sort()
        for ch in keys:
            is_sdev = ch.endswith('Sdev')

            if is_sdev and not save_sdev:
                continue

            hdu = fits.ImageHDU(images[ch], header=header, name='IMG' + ch)
            hdulist.append(hdu)

        hdulist.writeto(fname, overwrite=True)


def main_imager(args=None):
    """Main function."""
    import argparse

    description = ('Load a series of scans from a config file '
                   'and produce a map.')
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("file", nargs='?',
                        help="Load intermediate scanset from this file",
                        default=None, type=str)

    parser.add_argument("--sample-config", action='store_true', default=False,
                        help='Produce sample config file')

    parser.add_argument("-c", "--config", type=str, default=None,
                        help='Config file')

    parser.add_argument("--refilt", default=False,
                        action='store_true',
                        help='Re-run the scan filtering')

    parser.add_argument("--altaz", default=False,
                        action='store_true',
                        help='Do images in Az-El coordinates')

    parser.add_argument("--sub", default=False,
                        action='store_true',
                        help='Subtract the baseline from single scans')

    parser.add_argument("--interactive", default=False,
                        action='store_true',
                        help='Open the interactive display')

    parser.add_argument("--calibrate", type=str, default=None,
                        help='Calibration file')

    parser.add_argument("--nofilt", action='store_true', default=False,
                        help='Do not filter noisy channels')

    parser.add_argument("-g", "--global-fit", action='store_true',
                        default=False,
                        help='Perform global fitting of baseline')

    parser.add_argument("-e", "--exclude", nargs='+', default=None,
                        help='Exclude region from global fitting of baseline')

    parser.add_argument("--chans", type=str, default=None,
                        help=('Comma-separated channels to include in global '
                              'fitting (Ch0, Ch1, ...)'))

    parser.add_argument("-o", "--outfile", type=str, default=None,
                        help='Save intermediate scanset to this file.')

    parser.add_argument("-u", "--unit", type=str, default="Jy/beam",
                        help='Unit of the calibrated image. Jy/beam or '
                             'Jy/pixel')

    parser.add_argument("--debug", action='store_true', default=False,
                        help='Plot stuff and be verbose')

    parser.add_argument("--quick", action='store_true', default=False,
                        help='Calibrate after image creation, for speed '
                             '(bad when calibration depends on elevation)')

    parser.add_argument("--splat", type=str, default=None,
                        help=("Spectral scans will be scrunched into a single "
                              "channel containing data in the given frequency "
                              "range, starting from the frequency of the first"
                              " bin. E.g. '0:1000' indicates 'from the first "
                              "bin of the spectrum up to 1000 MHz above'. ':' "
                              "or 'all' for all the channels."))

    args = parser.parse_args(args)

    if args.sample_config:
        sample_config_file()
        sys.exit()

    outfile = args.outfile

    excluded = None
    if args.exclude is not None:
        nexc = len(args.exclude)
        if nexc % 3 != 0:
            raise ValueError("Exclusion region has to be specified as "
                             "centerX0, centerY0, radius0, centerX1, "
                             "centerY1, radius1, ... (in X,Y coordinates)")
        excluded = \
            np.array([np.float(e)
                      for e in args.exclude]).reshape((nexc // 3, 3))

    if args.file is not None:
        scanset = ScanSet(args.file, config_file=args.config)
        infile = args.file
        if outfile is None:
            outfile = infile
    else:
        if args.config is None:
            raise ValueError("Please specify the config file!")
        scanset = ScanSet(args.config, norefilt=not args.refilt,
                          freqsplat=args.splat, nosub=not args.sub,
                          nofilt=args.nofilt, debug=args.debug)
        infile = args.config

        if outfile is None:
            outfile = infile.replace('.ini', '_dump.hdf5')

    scanset.write(outfile, overwrite=True)

    if args.interactive:
        scanset.interactive_display()

    if args.global_fit:
        scanset.fit_full_images(excluded=excluded, chans=args.chans,
                                altaz=args.altaz)
        scanset.write(outfile.replace('.hdf5', '_baselinesub.hdf5'),
                      overwrite=True)

    scanset.save_ds9_images(save_sdev=True, calibration=args.calibrate,
                            map_unit=args.unit,
                            altaz=args.altaz, calibrate_scans=not args.quick)


def main_preprocess(args=None):
    """Preprocess the data."""
    import argparse

    description = ('Load a series of scans from a config file '
                   'and preprocess them, or preprocess a single scan.')
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("files", nargs='*',
                        help="Single files to preprocess",
                        default=None, type=str)

    parser.add_argument("-c", "--config", type=str, default=None,
                        help='Config file')

    parser.add_argument("--sub", default=False,
                        action='store_true',
                        help='Subtract the baseline from single scans')

    parser.add_argument("--interactive", default=False,
                        action='store_true',
                        help='Open the interactive display for each scan')

    parser.add_argument("--nofilt", action='store_true', default=False,
                        help='Do not filter noisy channels')

    parser.add_argument("--debug", action='store_true', default=False,
                        help='Plot stuff and be verbose')

    parser.add_argument("--splat", type=str, default=None,
                        help=("Spectral scans will be scrunched into a single "
                              "channel containing data in the given frequency "
                              "range, starting from the frequency of the first"
                              " bin. E.g. '0:1000' indicates 'from the first "
                              "bin of the spectrum up to 1000 MHz above'. ':' "
                              "or 'all' for all the channels."))

    args = parser.parse_args(args)

    if args.files is not None and args.files:
        for f in args.files:
            Scan(f, freqsplat=args.splat, nosub=not args.sub, norefilt=False,
                 debug=args.debug, interactive=args.interactive)
    else:
        if args.config is None:
            raise ValueError("Please specify the config file!")
        ScanSet(args.config, norefilt=False, freqsplat=args.splat,
                nosub=not args.sub, nofilt=args.nofilt, debug=args.debug,
                interactive=args.interactive)
