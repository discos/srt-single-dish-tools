"""Produce calibrated images.

``SDTimage`` is a script that, given a list of cross scans composing an
on-the-fly map, is able to calculate the map and save it in FITS format after
cleaning the data.
"""

import os
import sys
import warnings
import traceback
import copy
import functools
from datetime import datetime
from collections.abc import Iterable

from scipy.stats import binned_statistic_2d
import logging
import numpy as np
import astropy
from astropy import wcs
from astropy.table import Table, vstack, Column
from astropy.utils.metadata import MergeConflictWarning
import astropy.io.fits as fits
import astropy.units as u
from astropy.table.np_utils import TableMergeError
from astropy.time import Time

from .scan import Scan, list_scans
from .read_config import read_config, sample_config_file
from .utils import calculate_zernike_moments, calculate_beam_fom, HAS_MAHO
from .utils import compare_anything, ds9_like_log_scale, njit

from .io import chan_re, get_channel_feed, detect_data_kind
from .fit import linear_fun
from .interactive_filter import select_data
from .calibration import CalibratorTable
from .opacity import calculate_opacity
from .global_fit import fit_full_image
from .interactive_filter import create_empty_info

try:
    from tqdm import tqdm as show_progress
except ImportError:

    def show_progress(a, **kwargs):
        return a


try:
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    HAS_MPL = True
except ImportError:
    HAS_MPL = False


IMG_STR = "__img_dump_"
IMG_HOR_STR = "__img_hor_dump_"
IMG_VER_STR = "__img_ver_dump_"


__all__ = ["ScanSet"]


def _load_and_merge_subscans(indices_and_subscans):
    """

    Examples
    --------
    >>> t0 = Table({"time": [0, 1], "ra": [[0], [0.4]], "dec": [[0], [0.1]]})
    >>> t1 = Table({"time": [2, 3], "ra": [[0], [0.2]], "dec": [[0], [0.3]]})
    >>> t0.meta = {"FLAG": True, "filename": "puff.fits"}
    >>> t1.meta = {"FLAG": False, "filename": "puff2.fits"}
    >>> s = _load_and_merge_subscans([(1, t0), (3, t1)])
    >>> np.allclose(s["ra"], [[0], [0.2]])
    True
    """
    tables = []
    for i_s, s in indices_and_subscans:
        if "FLAG" in s.meta and s.meta["FLAG"]:
            logging.info("{} is flagged".format(s.meta["filename"]))
            continue
        s["Scan_id"] = i_s + np.zeros(len(s["time"]), dtype=int)

        ras = s["ra"][:, 0]
        decs = s["dec"][:, 0]

        ravar = (np.max(ras) - np.min(ras)) / np.cos(np.mean(decs))
        decvar = np.max(decs) - np.min(decs)
        s["direction"] = np.array(ravar > decvar, dtype=bool)

        # Remove conflicting keys
        for key in [
            "filename",
            "calibrator_directories",
            "skydip_directories",
            "list_of_directories",
            "mean_dec",
            "mean_ra",
            "max_dec",
            "max_ra",
            "min_dec",
            "min_ra",
            "dec_offset",
            "ra_offset",
            "RightAscension Offset",
            "Declination Offset",
        ]:
            if key in s.meta:
                s.meta.pop(key)

        tables.append(s)
    return merge_tables(tables)


def merge_tables(tables):
    """

    Examples
    --------
    >>> t0 = Table({"Ch0": [0, 1]})
    >>> t1 = Table({"Ch0": [[0, 1]]})
    >>> t0.meta["filename"] = "sss"
    >>> t1.meta["filename"] = "asd"
    >>> s = merge_tables([t0, t1])
    Traceback (most recent call last):
        ...
    astropy.table.np_utils.TableMergeError:...
        ...
    """
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", MergeConflictWarning)
            scan_table = vstack(tables)
    except TableMergeError as e:
        raise TableMergeError(
            "ERROR while merging tables. "
            "Tables:\n" + "\n".join([t.meta["filename"] for t in tables])
        )

    return scan_table


def all_lower(list_of_strings):
    return [s.lower() for s in list_of_strings]


def _coord_names(frame):
    """Name the coordinates based on the frame.

    Examples
    --------
    >>> _coord_names('icrs')[0]
    'ra'
    >>> _coord_names('altaz')[0]
    'delta_az'
    >>> _coord_names('sun')[0]
    'hpln'
    >>> _coord_names('galactic')[0]
    'glon'
    >>> _coord_names('ecliptic')[0]
    'elon'
    >>> _coord_names('turuturu')[0]
    Traceback (most recent call last):
       ...
    ValueError: turuturu: Unknown frame

    """
    if frame in ["icrs", "fk5"]:
        hor, ver = "ra", "dec"
    elif frame == "altaz":
        hor, ver = "delta_az", "delta_el"
    elif frame == "sun":
        hor, ver = "hpln", "hplt"
    elif frame == "galactic":
        hor, ver = "glon", "glat"
    elif frame == "ecliptic":
        hor, ver = "elon", "elat"
    else:
        raise ValueError(f"{frame}: Unknown frame")

    return hor, ver


def _wcs_ctype_names(frame, projection):
    """WCS ctype names

    Example
    -------
    >>> hor, ver = _wcs_ctype_names('icrs', 'ARC')
    >>> hor
    'RA---ARC'
    >>> ver
    'DEC--ARC'
    """
    hor, ver = _coord_names(frame)
    nproj = len(projection)
    nhor = len(hor)
    nver = len(ver)

    hor_str = f"{hor.upper()}{'-' * (8 - nhor - nproj)}{projection.upper()}"
    ver_str = f"{ver.upper()}{'-' * (8 - nver - nproj)}{projection.upper()}"

    return hor_str, ver_str


def _load_calibration(calibration, map_unit):
    caltable = CalibratorTable().read(calibration)
    caltable.update()
    caltable.compute_conversion_function(map_unit)

    if map_unit == "Jy/beam":
        conversion_units = u.Jy / u.ct
    elif map_unit in ["Jy/pixel", "Jy/sr"]:
        conversion_units = u.Jy / u.ct / u.steradian
    else:
        raise ValueError("Unit for calibration not recognized")
    return caltable, conversion_units


@njit
def _outlier_score(x):
    """Give a score to data series, larger if higher chance of outliers.

    Inspired by https://stackoverflow.com/questions/22354094/
    pythonic-way-of-detecting-outliers-in-one-dimensional-observation-data
    """
    xdiff = np.diff(x)
    good = xdiff != 0
    if not np.any(good):
        return 0
    xdiff = xdiff[good]
    if len(xdiff) < 2:
        return 0
    ref_dev = np.std(xdiff - np.median(xdiff))
    if ref_dev == 0.0:
        return 0

    median = np.median(x)
    diff = np.abs(x - median)
    return np.max(0.6745 * diff / ref_dev)


def outlier_score(x):
    """Give a score to data series, larger if higher chance of outliers.

    Inspired by https://stackoverflow.com/questions/22354094/
    pythonic-way-of-detecting-outliers-in-one-dimensional-observation-data
    """
    if len(x) == 0:
        return 0

    return _outlier_score(np.asarray(x))


class ScanSet(Table):
    def __init__(
        self,
        data=None,
        norefilt=True,
        config_file=None,
        freqsplat=None,
        nofilt=False,
        nosub=False,
        plot=False,
        debug=False,
        **kwargs,
    ):
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
            See :class:`srttools.scan.Scan`
        freqsplat : str
            See :class:`srttools.scan.interpret_frequency_range`
        nofilt : bool
            See :class:`srttools.scan.clean_scan_using_variability`
        nosub : bool
            See :class:`srttools.scan.Scan`

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
        self.norefilt = norefilt
        self.freqsplat = freqsplat
        self.images = None
        self.images_hor = None
        self.images_ver = None
        self.nofilt = nofilt
        self.nosub = nosub
        self.plot = plot
        self.debug = debug

        if data is None and config_file is None:
            pass
        elif isinstance(data, Iterable) and not isinstance(data, str):
            data = self._merge_input_data(data, config_file, **kwargs)
        elif isinstance(data, str) and data.endswith("hdf5"):
            data = Table.read(data)
        elif isinstance(data, str) and data.endswith("ini"):  # data is a config file
            data = self._read_data_from_config(data, **kwargs)
        elif not isinstance(data, Table):  # data needs to be a Table object
            raise ValueError(f"Invalid data: \n{data}")

        if config_file is not None:
            data.meta["config_file"] = config_file
            config = read_config(config_file)
            self.meta.update(config)

        super().__init__(data)

        if data is not None and not "x" in self.colnames:
            self.convert_coordinates()
        self.current = None
        self._scan_list = None
        self._chan_columns = None

    @property
    def chan_columns(self):
        if self._chan_columns is None:
            self._chan_columns = np.array([i for i in self.columns if chan_re.match(i)])
        return self._chan_columns

    @property
    def scan_list(self):
        if self._scan_list is None:
            self._scan_list = self.meta["scan_list"]
        return self._scan_list

    def _read_data_from_config(self, data, **kwargs):
        config_file = data
        config = read_config(config_file)

        self.meta.update(config)
        self.meta["config_file"] = config_file

        scan_list = self.list_scans()

        scan_list.sort()

        indices_and_subscans = self.load_scans(
            scan_list,
            debug=self.debug,
            freqsplat=self.freqsplat,
            nofilt=self.nofilt,
            nosub=self.nosub,
            plot=self.plot,
            **kwargs,
        )
        tables = _load_and_merge_subscans(indices_and_subscans)

        scan_table = merge_tables(tables)

        scan_table.meta["scan_list"] = scan_list
        return scan_table

    def _merge_input_data(self, data, config_file, **kwargs):
        alldata = [
            ScanSet(
                d,
                norefilt=self.norefilt,
                config_file=config_file,
                freqsplat=self.freqsplat,
                nofilt=self.nofilt,
                nosub=self.nosub,
                **kwargs,
            )
            for d in data
        ]

        scan_list = []
        max_scan_id = 0
        for d in alldata:
            scan_list += d.scan_list
            d["Scan_id"] += max_scan_id

            max_scan_id += len(d.scan_list)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", MergeConflictWarning)
            data = vstack(alldata)

        data.meta["scan_list"] = scan_list
        data.meta = data.meta
        return data

    def analyze_coordinates(self, frame="icrs"):
        """Save statistical information on coordinates."""
        hor, ver = _coord_names(frame)

        if "delta_az" not in self.colnames and frame == "altaz":
            self.calculate_delta_altaz()

        allhor = self[hor]
        allver = self[ver]
        hor_unit = self[hor].unit
        ver_unit = self[ver].unit

        # These seemingly useless float() calls are needed for serialize_meta
        self.meta["min_" + hor] = float(np.min(allhor)) * hor_unit
        self.meta["min_" + ver] = float(np.min(allver)) * ver_unit
        self.meta["max_" + hor] = float(np.max(allhor)) * hor_unit
        self.meta["max_" + ver] = float(np.max(allver)) * ver_unit

        self.meta["mean_" + hor] = (self.meta["max_" + hor] + self.meta["min_" + hor]) / 2
        self.meta["mean_" + ver] = (self.meta["max_" + ver] + self.meta["min_" + ver]) / 2

        if "reference_ra" not in self.meta:
            self.meta["reference_ra"] = self.meta["RA"]
        if "reference_dec" not in self.meta:
            self.meta["reference_dec"] = self.meta["Dec"]

    def list_scans(self, datadir=None, dirlist=None):
        """List all scans contained in the directory listed in config."""
        if datadir is None:
            datadir = self.meta["datadir"]
            dirlist = self.meta["list_of_directories"]
        return list_scans(datadir, dirlist)

    def get_opacity(self, datadir=None, dirlist=None):
        """List all scans contained in the directory listed in config."""
        self.opacities = {}
        if "skydip_directories" not in self.meta:
            return

        if datadir is None:
            datadir = self.meta["datadir"]
            dirlist = self.meta["skydip_directories"]
        scans = list_scans(datadir, dirlist)
        if len(scans) == 0:
            return

        for s in scans:
            if "summary.fits" in s:
                continue
            try:
                results = calculate_opacity(s, plot=False)
                self.opacities[results["time"]] = np.mean([results["Ch0"], results["Ch1"]])
            except KeyError as e:
                logging.warning("Error while processing {}: Missing key: {}".format(s, str(e)))

    def load_scans(self, scan_list, freqsplat=None, nofilt=False, debug=False, **kwargs):
        """Load the scans in the list one by ones."""
        for i, f in enumerate(show_progress(scan_list)):
            try:
                s = Scan(
                    f,
                    norefilt=self.norefilt,
                    debug=debug,
                    freqsplat=freqsplat,
                    nofilt=nofilt,
                    **kwargs,
                )
                yield i, s
            except KeyError as e:
                logging.warning("Error while processing {}: Missing key: {}".format(f, str(e)))
            except Exception as e:
                logging.warning(traceback.format_exc())
                logging.warning("Error while processing {}: {}".format(f, str(e)))

    def get_coordinates(self, frame="icrs"):
        """Give the coordinates as pairs of RA, DEC."""
        hor, ver = _coord_names(frame)
        return np.array(np.dstack([self[hor], self[ver]]))

    def get_obstimes(self):
        """Get :class:`astropy.time.Time` at the telescope location."""
        from .io import locations

        return Time(
            (self["time"]) * u.day,
            format="mjd",
            scale="utc",
            location=locations[self.meta["site"]],
        )

    def apply_user_filter(self, user_func=None, out_column=None):
        """Apply a user-supplied function as filter.

        Parameters
        ----------
        user_func : function
            This function needs to accept a :class:`srttools.imager.ScanSet`
            object as only argument.
            It has to return an array with the same length of
            a column of the input dataset
        out_column : str
            column where the results will be stored

        Returns
        -------
        retval : array
            the result of user_func
        """
        if user_func is None:
            raise ValueError("user_func needs to be specified")
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

        ref_coords = SkyCoord(
            ra=self.meta["reference_ra"],
            dec=self.meta["reference_dec"],
            obstime=self.get_obstimes(),
            location=locations[self.meta["site"]],
        )
        ref_altaz_coords = ref_coords.altaz
        ref_az = ref_altaz_coords.az.to(u.rad)
        ref_el = ref_altaz_coords.alt.to(u.rad)

        self.meta["reference_delta_az"] = 0 * u.rad
        self.meta["reference_delta_el"] = 0 * u.rad
        self["delta_az"] = np.zeros_like(self["az"])
        self["delta_el"] = np.zeros_like(self["el"])
        for f in range(len(self["el"][0, :])):
            self["delta_az"][:, f] = (self["az"][:, f] - ref_az) * np.cos(ref_el)
            self["delta_el"][:, f] = self["el"][:, f] - ref_el

        if HAS_MPL:
            fig1 = plt.figure("adsfasdfasd")
            plt.plot(np.degrees(self["delta_az"]), np.degrees(self["delta_el"]))
            plt.xlabel("Delta Azimuth (deg)")
            plt.ylabel("Delta Elevation (deg)")

            plt.savefig("delta_altaz.png")
            plt.close(fig1)

            fig2 = plt.figure("adsfasdf")
            plt.plot(np.degrees(self["az"]), np.degrees(self["el"]))
            plt.plot(np.degrees(ref_az), np.degrees(ref_el))
            plt.xlabel("Azimuth (deg)")
            plt.ylabel("Elevation (deg)")
            plt.savefig("altaz_with_src.png")
            plt.close(fig2)

    def create_wcs(self, frame="icrs"):
        """Create a wcs object from the pointing information."""
        hor, ver = _coord_names(frame)
        pixel_size = self.meta["pixel_size"]
        self.wcs = wcs.WCS(naxis=2)

        if "max_" + hor not in self.meta:
            self.analyze_coordinates(frame)

        delta_hor = self.meta["max_" + hor] - self.meta["min_" + hor]
        if frame != "sun":
            delta_hor *= np.cos(self.meta["reference_" + ver])
        delta_ver = self.meta["max_" + ver] - self.meta["min_" + ver]

        # npix >= 1!
        npix_hor = np.ceil(delta_hor / pixel_size)
        npix_ver = np.ceil(delta_ver / pixel_size)

        self.meta["npix"] = np.array([npix_hor, npix_ver])

        # the first pixel is starts from 1, 1!
        self.wcs.wcs.crpix = self.meta["npix"] / 2 + 1

        # TODO: check consistency of units
        # Here I'm assuming all angles are radians
        # crval = np.array([self.meta['reference_' + hor].to(u.rad).value,
        #                   self.meta['reference_' + ver].to(u.rad).value])
        crhor = np.mean([self.meta["max_" + hor].value, self.meta["min_" + hor].value])
        crver = np.mean([self.meta["max_" + ver].value, self.meta["min_" + ver].value])
        crval = np.array([crhor, crver])

        self.wcs.wcs.crval = np.degrees(crval)

        cdelt = np.array([-pixel_size.to(u.rad).value, pixel_size.to(u.rad).value])
        self.wcs.wcs.cdelt = np.degrees(cdelt)

        hor_str, ver_str = _wcs_ctype_names(frame, self.meta["projection"])
        self.wcs.wcs.ctype = [hor_str, ver_str]

    def convert_coordinates(self, frame="icrs"):
        """Convert the coordinates from sky to pixel."""
        hor, ver = _coord_names(frame)

        self.create_wcs(frame)

        self["x"] = np.zeros_like(self[hor])
        self["y"] = np.zeros_like(self[ver])
        coords = np.degrees(self.get_coordinates(frame=frame))
        for f in range(len(self[hor][0, :])):
            pixcrd = self.wcs.all_world2pix(coords[:, f], 0.5)

            self["x"][:, f] = pixcrd[:, 0] + 0.5
            self["y"][:, f] = pixcrd[:, 1] + 0.5
        self["x"].meta["frame"] = frame
        self["y"].meta["frame"] = frame

    def calculate_images(
        self,
        no_offsets=False,
        frame="icrs",
        calibration=None,
        elevation=None,
        map_unit="Jy/beam",
        calibrate_scans=False,
        direction=None,
        onlychans=None,
    ):
        """Obtain image from all scans.

        Other parameters
        ----------------
        no_offsets: bool
            use positions from feed 0 for all feeds.
        frame: str
            One of ``icrs``, ``altaz``, ``sun``, ``galactic``, default ``icrs``
        calibration: CalibratorTable
            Optional Calibrator table for calibration, default ``None``
        elevation: Angle
            Optional elevation angle. Defaults to mean elevation
        map_unit: str
            Only used if ``calibration`` is not ``None``. One of ``Jy/beam``
            or ``Jy/pixel``
        calibrate_scans: bool
            Calibrate from subscans instead of from the binned image. Slower
            but more precise
        direction: int
            Optional: only select horizontal or vertical scans for this image.
            0 if horizontal, 1 if vertical, defaults to ``None`` that means
            that all scans will be used.
        onlychans: List
            List of channels for which images are calculated. If defaults to
            all channels
        """
        if frame != self["x"].meta["frame"]:
            self.convert_coordinates(frame)

        images = {}

        xbins = np.linspace(0, self.meta["npix"][0], int(self.meta["npix"][0] + 1))
        ybins = np.linspace(0, self.meta["npix"][1], int(self.meta["npix"][1] + 1))

        for ch in self.chan_columns:
            if direction is None:
                logging.info("Calculating image in channel {}".format(ch))
            else:
                dir_string = "horizontal" if direction == 1 else "vertical"
                logging.info("Calculating image in channel {}, {}".format(ch, dir_string))
            if (
                onlychans is not None
                and ch not in onlychans
                and self.images is not None
                and ch in self.images.keys()
            ):
                images[ch] = self.images[ch]
                images["{}-Sdev".format(ch)] = self.images["{}-Sdev".format(ch)]
                images["{}-EXPO".format(ch)] = self.images["{}-EXPO".format(ch)]
                images["{}-Outliers".format(ch)] = self.images["{}-Outliers".format(ch)]
                continue

            feed = get_channel_feed(ch)

            if elevation is None:
                elevation = np.mean(self["el"][:, feed])

            if "{}-filt".format(ch) in self.keys():
                good = self["{}-filt".format(ch)]
            else:
                good = np.ones(len(self[ch]), dtype=bool)

            if direction == 0:
                good = good & self["direction"]
            elif direction == 1:
                good = good & np.logical_not(self["direction"])

            expomap, _, _ = np.histogram2d(
                self["x"][:, feed][good],
                self["y"][:, feed][good],
                bins=[xbins, ybins],
            )

            counts = np.array(self[ch][good])

            if calibration is not None and calibrate_scans:
                caltable, conversion_units = _load_calibration(calibration, map_unit)
                (
                    area_conversion,
                    final_unit,
                ) = self._calculate_calibration_factors(map_unit)

                (
                    Jy_over_counts,
                    Jy_over_counts_err,
                ) = conversion_units * caltable.Jy_over_counts(
                    channel=ch,
                    map_unit=map_unit,
                    elevation=self["el"][:, feed][good],
                )

                counts = counts * u.ct * area_conversion * Jy_over_counts
                counts = counts.to(final_unit).value

            filtered_x = self["x"][:, feed][good]
            filtered_y = self["y"][:, feed][good]

            img, _, _ = np.histogram2d(filtered_x, filtered_y, bins=[xbins, ybins], weights=counts)

            img_sq, _, _ = np.histogram2d(
                filtered_x,
                filtered_y,
                bins=[xbins, ybins],
                weights=counts**2,
            )

            img_outliers, _, _, _ = binned_statistic_2d(
                filtered_x,
                filtered_y,
                counts,
                statistic=outlier_score,
                bins=[xbins, ybins],
            )

            good = expomap > 0
            mean = img.copy()
            mean[good] /= expomap[good]
            # For Numpy vs FITS image conventions...
            images[ch] = mean.T
            img_sdev = img_sq
            img_sdev[good] = img_sdev[good] / expomap[good] - mean[good] ** 2

            img_sdev[good] = np.sqrt(img_sdev[good])
            if calibration is not None and calibrate_scans:
                cal_rel_err = np.mean(Jy_over_counts_err / Jy_over_counts).value
                img_sdev += mean * cal_rel_err

            images["{}-Sdev".format(ch)] = img_sdev.T
            images["{}-EXPO".format(ch)] = expomap.T
            images["{}-Outliers".format(ch)] = img_outliers.T

        if direction is None:
            self.images = images
        elif direction == 0:
            self.images_hor = images
        elif direction == 1:
            self.images_ver = images

        if calibration is not None and not calibrate_scans:
            self.calibrate_images(
                calibration,
                elevation=elevation,
                map_unit=map_unit,
                direction=direction,
            )

        return images

    def destripe_images(self, niter=10, npix_tol=None, **kwargs):
        from .destripe import destripe_wrapper

        if self.images is None:
            images = self.calculate_images(**kwargs)
        else:
            images = self.images

        destriped = {}
        for ch in self.chan_columns:
            if ch in images:
                destriped[ch + "_dirty"] = images[ch]

        if self.images_hor is None:
            self.calculate_images(direction=0, **kwargs)
        if self.images_ver is None:
            self.calculate_images(direction=1, **kwargs)

        images_hor, images_ver = self.images_hor, self.images_ver

        for ch in images_hor:
            if "Sdev" in ch:
                destriped[ch] = (images_hor[ch] ** 2 + images_ver[ch] ** 2) ** 0.5
                continue
            if "EXPO" in ch:
                destriped[ch] = images_hor[ch] + images_ver[ch]
                continue
            if "Outlier" in ch:
                destriped[ch] = images_hor[ch] + images_ver[ch]
                continue

            destriped[ch] = destripe_wrapper(
                images_hor[ch],
                images_ver[ch],
                niter=niter,
                npix_tol=npix_tol,
                expo_hor=images_hor[ch + "-EXPO"],
                expo_ver=images_ver[ch + "-EXPO"],
                label=ch,
            )

        for ch in destriped:
            self.images[ch] = destriped[ch]

        return self.images

    def scrunch_images(self, bad_chans=[]):
        """Sum the images from all channels."""
        total_expo = 0
        total_img = 0
        total_sdev = 0
        count = 0
        lower_bad_chans = all_lower(bad_chans)
        for ch in self.chan_columns:
            if ch.lower() in lower_bad_chans:
                logging.info("Discarding {}".format(ch))
                continue
            total_expo += self.images["{}-EXPO".format(ch)]
            total_sdev += self.images["{}-Sdev".format(ch)] ** 2
            total_img += self.images[ch]
            count += 1
        total_sdev = total_sdev**0.5 / count
        total_img /= count

        total_images = {
            "TOTAL": total_img,
            "TOTAL-Sdev": np.sqrt(total_sdev),
            "TOTAL-EXPO": total_expo,
        }
        self.images.update(total_images)
        return total_images

    def fit_full_images(
        self,
        chans=None,
        fname=None,
        save_sdev=False,
        no_offsets=False,
        frame="icrs",
        calibration=None,
        excluded=None,
        par=None,
        map_unit="Jy/beam",
    ):
        """Flatten the baseline with a global fit.

        Fit a linear trend to each scan to minimize the scatter in an image
        """

        if self.images is None:
            self.calculate_images(
                no_offsets=no_offsets,
                frame=frame,
                calibration=calibration,
                map_unit=map_unit,
            )

        if chans is not None:
            chans = chans.split(",")
        else:
            chans = self.chan_columns

        for ch in chans:
            logging.info("Fitting channel {}".format(ch))
            feed = get_channel_feed(ch)
            self[ch + "_save"] = self[ch].copy()
            self[ch] = Column(fit_full_image(self, chan=ch, feed=feed, excluded=excluded, par=par))
            self[ch].meta = self[ch + "_save"].meta

        self.calculate_images(
            no_offsets=no_offsets,
            frame=frame,
            calibration=calibration,
            map_unit=map_unit,
        )

    def _calculate_calibration_factors(self, map_unit):
        if map_unit == "Jy/beam":
            area_conversion = 1
            final_unit = u.Jy
        elif map_unit == "Jy/sr":
            area_conversion = 1
            final_unit = u.Jy / u.sr
        elif map_unit == "Jy/pixel":
            area_conversion = self.meta["pixel_size"] ** 2
            final_unit = u.Jy
        return area_conversion, final_unit

    def calibrate_images(
        self,
        calibration,
        elevation=np.pi / 4,
        map_unit="Jy/beam",
        direction=None,
    ):
        """Calibrate the images."""
        if self.images is None:
            self.calculate_images(direction=direction)

        if direction == 0:
            images = self.images_hor
        elif direction == 1:
            images = self.images_ver
        else:
            images = self.images

        caltable, conversion_units = _load_calibration(calibration, map_unit)

        for ch in self.chan_columns:
            Jy_over_counts, Jy_over_counts_err = (
                caltable.Jy_over_counts(channel=ch, map_unit=map_unit, elevation=elevation)
                * conversion_units
            )

            if np.isnan(Jy_over_counts):
                warnings.warn("The Jy/counts factor is nan")
                continue
            A = images[ch].copy() * u.ct
            eA = images["{}-Sdev".format(ch)].copy() * u.ct

            images["{}-RAW".format(ch)] = images["{}".format(ch)].copy()
            images["{}-RAW-Sdev".format(ch)] = images["{}-Sdev".format(ch)].copy()
            images["{}-RAW-EXPO".format(ch)] = images["{}-EXPO".format(ch)].copy()
            bad = eA != eA
            A[bad] = 1 * u.ct
            eA[bad] = 0 * u.ct

            bad = np.logical_or(A == 0, A != A)
            A[bad] = 1 * u.ct
            eA[bad] = 0 * u.ct

            B = Jy_over_counts
            eB = Jy_over_counts_err

            area_conversion, final_unit = self._calculate_calibration_factors(map_unit)

            C = A * area_conversion * Jy_over_counts
            C[bad] = 0

            images[ch] = C.to(final_unit).value

            eC = C * (eA / A + eB / B)

            images["{}-Sdev".format(ch)] = eC.to(final_unit).value

        if direction == 0:
            self.images_hor = images
        elif direction == 1:
            self.images_ver = images
        else:
            self.images = images

    def interactive_display(self, ch=None, recreate=False, test=False):
        """Modify original scans from the image display."""
        from .interactive_filter import ImageSelector

        if not HAS_MPL:
            raise ImportError("interactive_display: " "matplotlib is not installed")

        if self.images is None:
            recreate = True

        self.display_instructions = """
        -------------------------------------------------------------

        Imageactive display.

        You see here two images. The left one gives, for each bin, a number
        measuring the probability of outliers (based on the median absolute
        deviation if there are >10 scans per bin, and on the standard deviation
        otherwise), The right one is the output image of the processing.
        The right image is normalized with a ds9-like log scale.

        -------------------------------------------------------------

        Point the mouse on a pixel in the Outlier image and press a key:

        a    open a window to filter all scans passing through this pixel
        h    print help
        q    quit

        -------------------------------------------------------------
        """
        print(self.display_instructions)

        if ch is None:
            chs = self.chan_columns
        else:
            chs = [ch]
        if test:
            chs = ["Feed0_RCP"]

        for ch in chs:
            if recreate:
                self.calculate_images(onlychans=ch)
            fig = plt.figure("Imageactive Display - " + ch)
            gs = GridSpec(1, 2)
            ax = fig.add_subplot(gs[0])
            ax.set_title("Outlier plot")
            ax2 = fig.add_subplot(gs[1])
            ax2.set_title("Draft image")
            imgch = ch

            expo = np.mean(self.images["{}-EXPO".format(ch)])
            mean_expo = np.mean(expo[expo > 0])

            stats_for_outliers = "Outliers" if mean_expo > 6 else "Sdev"
            sdevch = "{}-{}".format(ch, stats_for_outliers)
            if "{}-RAW".format(ch) in self.images.keys():
                imgch = "{}-RAW".format(ch)
                if stats_for_outliers == "Sdev":
                    sdevch = "{}-RAW-{}".format(ch, stats_for_outliers)
            img = ds9_like_log_scale(self.images[imgch])
            ax2.imshow(
                img,
                origin="lower",
                vmin=np.percentile(img, 20),
                cmap="gnuplot2",
                interpolation="nearest",
            )

            img = self.images[sdevch].copy()
            self.current = ch
            bad = np.logical_or(img == 0, img != img)
            img[bad] = np.mean(img[np.logical_not(bad)])
            fun = functools.partial(self.rerun_scan_analysis, test=test)

            imgsel = ImageSelector(img, ax, fun=fun, test=test)
        return imgsel

    def rerun_scan_analysis(self, x, y, key, test=False):
        """Rerun the analysis of single scans."""
        logging.debug("{} {} {}".format(x, y, key))
        if key == "a":
            self.reprocess_scans_through_pixel(x, y, test=test)
        elif key == "h":
            print(self.display_instructions)
        elif key == "v":
            pass

    def reprocess_scans_through_pixel(self, x, y, test=False):
        """Given a pixel in the image, find all scans passing through it."""
        ch = self.current

        (
            ra_xs,
            ra_ys,
            dec_xs,
            dec_ys,
            scan_ids,
            ra_masks,
            dec_masks,
            vars_to_filter,
        ) = self.find_scans_through_pixel(x, y, test=test)

        if ra_xs != {}:
            empty = create_empty_info(ra_xs.keys())
            info = select_data(
                ra_xs,
                ra_ys,
                masks=ra_masks,
                xlabel="RA",
                title="RA",
                test=test,
            )

            if not compare_anything(empty, info):
                for sname in info.keys():
                    self.update_scan(
                        sname,
                        scan_ids[sname],
                        vars_to_filter[sname],
                        info[sname]["zap"],
                        info[sname]["fitpars"],
                        info[sname]["FLAG"],
                    )

        if dec_xs != {}:
            empty = create_empty_info(ra_xs.keys())
            info = select_data(
                dec_xs,
                dec_ys,
                masks=dec_masks,
                xlabel="Dec",
                title="Dec",
                test=test,
            )

            if not compare_anything(empty, info):
                for sname in info.keys():
                    self.update_scan(
                        sname,
                        scan_ids[sname],
                        vars_to_filter[sname],
                        info[sname]["zap"],
                        info[sname]["fitpars"],
                        info[sname]["FLAG"],
                    )

        # Only recreate images if there were changes!
        display = self.interactive_display(ch=ch, recreate=(dec_xs != {} or ra_xs != {}), test=test)
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
            ch = "Feed0_RCP"

        feed = get_channel_feed(ch)

        # Select data inside the pixel +- 1

        good_entries = np.logical_and(
            np.abs(self["x"][:, feed] - x) < 1,
            np.abs(self["y"][:, feed] - y) < 1,
        )

        sids = list(set(self["Scan_id"][good_entries]))

        for sid in sids:
            sname = self.scan_list[sid]
            try:
                s = Scan(sname)
            except Exception:
                logging.warning("Errors while opening scan {}".format(sname))
                continue
            try:
                chan_mask = s["{}-filt".format(ch)]
            except Exception:
                chan_mask = np.zeros_like(s[ch])

            scan_ids[sname] = sid
            ras = s["ra"][:, feed]
            decs = s["dec"][:, feed]

            z = s[ch]

            ravar = np.max(ras) - np.min(ras)
            decvar = np.max(decs) - np.min(decs)
            if ravar > decvar:
                vars_to_filter[sname] = "ra"
                ra_xs[sname] = ras
                ra_ys[sname] = z
                ra_masks[sname] = chan_mask
            else:
                vars_to_filter[sname] = "dec"
                dec_xs[sname] = decs
                dec_ys[sname] = z
                dec_masks[sname] = chan_mask

        return (
            ra_xs,
            ra_ys,
            dec_xs,
            dec_ys,
            scan_ids,
            ra_masks,
            dec_masks,
            vars_to_filter,
        )

    def update_scan(self, sname, sid, dim, zap_info, fit_info, flag_info, test=False):
        """Update a scan in the scanset after filtering."""
        ch = self.current
        if test:
            ch = "Feed0_RCP"
        feed = get_channel_feed(ch)
        mask = self["Scan_id"] == sid
        try:
            logging.info("Updating scan {}".format(sname))
            s = Scan(sname)
        except Exception as e:
            warnings.warn("Impossible to write to scan {}".format(sname))
            warnings.warn(str(e))
            return

        resave = False
        if len(zap_info.xs) > 0:
            resave = True
            xs = zap_info.xs
            good = np.ones(len(s[dim]), dtype=bool)
            if len(xs) >= 2:
                intervals = list(zip(xs[:-1:2], xs[1::2]))
                for i in intervals:
                    i = sorted(i)
                    good[np.logical_and(s[dim][:, feed] >= i[0], s[dim][:, feed] <= i[1])] = False
            s["{}-filt".format(ch)] = good
            self["{}-filt".format(ch)][mask] = good

        if len(fit_info) > 1:
            resave = True
            sub = linear_fun(s[dim][:, feed], *fit_info)
            s[ch][:] = np.array(s[ch]) - sub
            # TODO: make it channel-independent
            s.meta["backsub"] = True
            self[ch][mask] = s[ch]

        # TODO: make it channel-independent
        if flag_info is not None:
            resave = True
            s.meta["FLAG"] = flag_info
            flag_array = np.zeros(len(s[dim]), dtype=bool) + flag_info
            for c in self.chan_columns:
                self["{}-filt".format(c)][mask] = np.logical_not(flag_array)
                s["{}-filt".format(c)] = np.logical_not(flag_array)

        if resave:
            s.save()

    def barycenter_times(self):
        """Create barytime column with observing times converted to TDB."""
        obstimes_tdb = self.get_obstimes().tdb.mjd
        self["barytime"] = obstimes_tdb
        return obstimes_tdb

    def write(self, fname, *args, **kwargs):
        """Same as Table.write, but adds path information for HDF5."""

        self.update_meta_with_images()

        try:
            kwargs["serialize_meta"] = kwargs.pop("serialize_meta", True)
            super(ScanSet, self).write(fname, *args, **kwargs)
        except astropy.io.registry.IORegistryError as e:
            raise astropy.io.registry.IORegistryError(fname + ": " + str(e))

    def update_meta_with_images(self):
        if self.images is not None:
            for key in self.images.keys():
                self.meta[IMG_STR + key] = self.images[key]
        if self.images_hor is not None:
            for key in self.images_hor.keys():
                self.meta[IMG_HOR_STR + key] = self.images_hor[key]
        if self.images_ver is not None:
            for key in self.images_ver.keys():
                self.meta[IMG_VER_STR + key] = self.images_ver[key]

    def read_images_from_meta(self):
        for key in self.meta.keys():
            logging.info("Caught key press: {}".format(key))
            if IMG_STR in key:
                self.images = {}
                self.images[key.replace(IMG_STR, "")] = self.meta[key]
            elif IMG_HOR_STR in key:
                self.images_hor = {}
                self.images_hor[key.replace(IMG_HOR_STR, "")] = self.meta[key]
            elif IMG_VER_STR in key:
                self.images_ver = {}
                self.images_ver[key.replace(IMG_VER_STR, "")] = self.meta[key]

    def calculate_zernike_moments(
        self, im, cm=None, radius=0.3, norder=8, label=None, use_log=False
    ):
        """Calculate the Zernike moments of the image.

        These moments are useful to single out asymmetries in the image:
        for example, when characterizing the beam of the radio telescope using
        a map of a calibrator, it is useful to calculate these moments to
        understand if the beam is radially symmetric or has distorted side
        lobes.

        Parameters
        ----------
        im : 2-d array
            The image to be analyzed
        cm : [int, int]
            'Center of mass' of the image
        radius : float
            The radius around the center of mass, in percentage of the image
            size (0 <= radius <= 0.5)
        norder : int
            Maximum order of the moments to calculate

        Returns
        -------
        moments_dict : dict
            Dictionary containing the order, the sub-index and the moment, e.g.
            {0: {0: 0.3}, 1: {0: 1e-16}, 2: {0: 0.95, 2: 6e-19}, ...}
            Moments are symmetrical, so only the unique values are reported.
        """
        if isinstance(im, str):
            im = self.images[im]

        return calculate_zernike_moments(
            im,
            cm=cm,
            radius=radius,
            norder=norder,
            label=label,
            use_log=use_log,
        )

    def calculate_beam_fom(
        self,
        im,
        cm=None,
        radius=0.3,
        label=None,
        use_log=False,
        show_plot=False,
    ):
        """Calculate various figures of merit (FOMs) in an image.

        These FOMs are useful to single out asymmetries in a beam shape:
        for example, when characterizing the beam of the radio telescope using
        a map of a calibrator, it is useful to understand if there are lobes
        appearing only in one direction.

        Parameters
        ----------
        im : 2-d array
            The image to be analyzed

        Other parameters
        ----------------
        cm : [int, int]
            'Center of mass' of the image
        radius : float
            The radius around the center of mass, in percentage of the image
            size (0 <= radius <= 0.5)
        use_log: bool
            Rescale the image to a log scale before calculating the
            coefficients.
            The scale is the same documented in the ds9 docs, for consistency.
            After normalizing the image from 0 to 1, the log-rescaled image is
            log(ax + 1) / log a, with ``x`` the normalized image and ``a`` a
            constant fixed here at 1000
        show_plot : bool, default False
            show the plots immediately

        Returns
        -------
        results_dict : dict
            Dictionary containing the results
        """
        if isinstance(im, str):
            im = self.images[im]

        return calculate_beam_fom(
            im,
            cm=cm,
            radius=radius,
            label=label,
            use_log=use_log,
            show_plot=show_plot,
        )

    def save_ds9_images(
        self,
        fname=None,
        save_sdev=False,
        scrunch=False,
        no_offsets=False,
        frame="icrs",
        calibration=None,
        map_unit="Jy/beam",
        calibrate_scans=False,
        destripe=False,
        npix_tol=None,
        bad_chans=[],
    ):
        """Save a ds9-compatible file with one image per extension."""
        if fname is None:
            tail = ".fits"
            if frame != "icrs":
                tail = f"_{frame}.fits"
            if scrunch:
                tail = tail.replace(".fits", "_scrunch.fits")
            if calibration is not None:
                tail = tail.replace(".fits", "_cal.fits")
            if destripe:
                tail = tail.replace(".fits", "_destripe.fits")
            fname = self.meta["config_file"].replace(".ini", tail)

        if destripe:
            logging.info("Destriping....")
            images = self.destripe_images(
                no_offsets=no_offsets,
                frame=frame,
                calibration=calibration,
                map_unit=map_unit,
                npix_tol=npix_tol,
                calibrate_scans=calibrate_scans,
            )
        else:
            images = self.calculate_images(
                no_offsets=no_offsets,
                frame=frame,
                calibration=calibration,
                map_unit=map_unit,
                calibrate_scans=calibrate_scans,
            )

        if scrunch:
            self.scrunch_images(bad_chans=bad_chans)

        self.create_wcs(frame)

        hdulist = fits.HDUList()

        header = self.wcs.to_header()
        if map_unit == "Jy/beam" and calibration is not None:
            caltable = CalibratorTable.read(calibration)
            beam, _ = caltable.beam_width()
            std_to_fwhm = np.sqrt(8 * np.log(2))
            header["bmaj"] = np.degrees(beam) * std_to_fwhm
            header["bmin"] = np.degrees(beam) * std_to_fwhm
            header["bpa"] = 0

        if calibration is not None:
            header["bunit"] = map_unit

        if "dsun" in self.colnames:
            print("SAVING DSUN")
            header["dsun_obs"] = np.mean(self["dsun"])
            header["dsun_ref"] = 149597870700.0

        for (
            key
        ) in "ANTENNA,site,RightAscension,Declination,backend,receiver,DATE,Project_Name,SiteLongitude,SiteLatitude,SiteHeight,ScheduleName".split(
            ","
        ):
            if key not in self.meta:
                warnings.warn(f"Key {key} not found in metadata")
                continue
            headkey = key
            if len(key) > 8:
                headkey = f"HIERARCH {key}"
            header[headkey] = self.meta[key]

        firstchan = self.chan_columns[0]
        for key in "frequency,bandwidth,local_oscillator".split(","):
            if key not in self[firstchan].meta[key]:
                warnings.warn(f"Key {key} not in the {firstchan} metadata")

            val = self[firstchan].meta[key].to(u.GHz)
            headkey = key
            if len(key) > 8:
                headkey = f"HIERARCH {key}"
            header[headkey] = (val.value, val.unit)

        header["SOURCE"] = self.meta["SOURCE"].replace("_RA", "").replace("_DEC", "")

        header["CREATOR"] = "SDT"
        ut = Time(datetime.utcnow(), scale="utc")
        header["COMMENT"] = f"Made with the SRT Single-Dish Tools on UT {ut.fits}"
        hdu = fits.PrimaryHDU(header=header)
        hdulist.append(hdu)

        keys = list(images.keys())
        keys.sort()

        header_mod = copy.deepcopy(header)
        for ch in keys:
            is_sdev = ch.endswith("Sdev")
            is_expo = "EXPO" in ch
            is_outl = "Outliers" in ch
            is_stokes = ("Q" in ch) or ("U" in ch)

            do_moments = not (is_sdev or is_expo or is_stokes or is_outl)
            do_moments = do_moments and frame and HAS_MAHO

            if is_sdev and not save_sdev:
                continue

            if do_moments:
                moments_dict = self.calculate_zernike_moments(
                    images[ch],
                    cm=None,
                    radius=0.3,
                    norder=8,
                    label=ch,
                    use_log=True,
                )
                for k in moments_dict.keys():
                    if k == "Description":
                        continue
                    for k1 in moments_dict[k].keys():
                        header_mod["ZK_{:02d}_{:02d}".format(k, k1)] = moments_dict[k][k1]
                moments_dict = self.calculate_beam_fom(
                    images[ch], cm=None, radius=0.3, label=ch, use_log=True
                )
                for k in moments_dict.keys():
                    if k == "Description":
                        continue
                    logging.info("FOM_{}: {}".format(k, moments_dict[k]))
                    # header_mod['FOM_{}'.format(k)] = moments_dict[k]

            hdu = fits.ImageHDU(images[ch], header=header_mod, name="IMG" + ch)
            hdulist.append(hdu)

        hdulist.writeto(fname, overwrite=True)


def _excluded_regions_from_args(args_exclude):
    excluded_xy = None
    excluded_radec = None
    if args_exclude is not None and not len(args_exclude) == 1:
        nexc = len(args_exclude)
        if nexc % 3 != 0:
            raise ValueError(
                "Exclusion region has to be specified as "
                "centerX0, centerY0, radius0, centerX1, "
                "centerY1, radius1, ... (in X,Y coordinates)"
            )
        excluded_xy = np.array([float(e) for e in args_exclude]).reshape((nexc // 3, 3))
        excluded_radec = None
    elif args_exclude is not None:
        import pyregion

        regions = pyregion.open(args_exclude[0])
        nregs = len(regions)
        excluded_xy = []
        excluded_radec = []
        for i in range(nregs):
            region = regions[i]
            if region.name != "circle":
                logging.warning("Only circular regions are allowed!")
                continue
            if region.coord_format == "fk5":
                excluded_radec.append(np.radians(region.coord_list))
            elif region.coord_format == "image":
                excluded_xy.append(region.coord_list)
            else:
                logging.warning("Only regions in fk5 or image coordinates are allowed!")
                continue
    return excluded_xy, excluded_radec


def main_imager(args=None):
    """Main function."""
    import argparse

    description = "Load a series of scans from a config file " "and produce a map."
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "file",
        nargs="?",
        help="Load intermediate scanset from this file",
        default=None,
        type=str,
    )

    parser.add_argument(
        "--sample-config",
        action="store_true",
        default=False,
        help="Produce sample config file",
    )

    parser.add_argument("-c", "--config", type=str, default=None, help="Config file")

    parser.add_argument(
        "--refilt",
        default=False,
        action="store_true",
        help="Re-run the scan filtering",
    )

    parser.add_argument(
        "--altaz",
        default=False,
        action="store_true",
        help="Do images in Az-El coordinates (deprecated in favor of --frame " "altaz)",
    )

    parser.add_argument(
        "--sub",
        default=False,
        action="store_true",
        help="Subtract the baseline from single scans",
    )

    parser.add_argument(
        "--interactive",
        default=False,
        action="store_true",
        help="Open the interactive display",
    )

    parser.add_argument("--calibrate", type=str, default=None, help="Calibration file")

    parser.add_argument(
        "--nofilt",
        action="store_true",
        default=False,
        help="Do not filter noisy channels",
    )

    parser.add_argument(
        "-g",
        "--global-fit",
        action="store_true",
        default=False,
        help="Perform global fitting of baseline",
    )

    parser.add_argument(
        "-e",
        "--exclude",
        nargs="+",
        default=None,
        help="Exclude region from global fitting of baseline "
        "and baseline subtraction. It can be specified "
        "as X1, Y1, radius1, X2, Y2, radius2 in image "
        "coordinates or as a ds9-compatible region file "
        "in image or fk5 coordinates containing circular "
        "regions to be excluded. Currently, baseline "
        "subtraction only takes into account fk5 "
        "coordinates and global fitting image coordinates"
        ". This will change in the future.",
    )

    parser.add_argument(
        "--chans",
        type=str,
        default=None,
        help=(
            "Comma-separated channels to include in global " "fitting (Feed0_RCP, Feed0_LCP, ...)"
        ),
    )

    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        default=None,
        help="Save intermediate scanset to this file.",
    )

    parser.add_argument(
        "-u",
        "--unit",
        type=str,
        default="Jy/beam",
        help="Unit of the calibrated image. Jy/beam or " "Jy/pixel",
    )

    parser.add_argument(
        "--frame",
        type=str,
        default="icrs",
        choices=["icrs", "altaz", "sun"],
        help="Reference frame for the image. One of icrs, altaz, sun",
    )

    parser.add_argument(
        "--destripe",
        action="store_true",
        default=False,
        help="Destripe the image",
    )

    parser.add_argument(
        "--npix-tol",
        type=int,
        default=None,
        help="Number of pixels with zero exposure to tolerate"
        " when destriping the image, or the full row or "
        "column is discarded."
        " Default None, meaning that the image will be"
        " destriped as a whole",
    )

    parser.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="Plot stuff and be verbose",
    )

    parser.add_argument(
        "--quick",
        action="store_true",
        default=False,
        help="Calibrate after image creation, for speed "
        "(bad when calibration depends on elevation)",
    )

    parser.add_argument(
        "--scrunch-channels",
        action="store_true",
        default=False,
        help="Sum all the images from the single channels into" " one.",
    )

    parser.add_argument(
        "--nosave",
        action="store_true",
        default=False,
        help="Do not save the hdf5 intermediate files when" "loading subscans.",
    )

    parser.add_argument(
        "--noplot",
        action="store_true",
        default=False,
        help="Do not produce diagnostic plots for data " "processing",
    )

    parser.add_argument(
        "--bad-chans",
        default="",
        type=str,
        help="Channels to be discarded when scrunching, "
        "separated by a comma (e.g. "
        "--bad-chans Feed2_RCP,Feed3_RCP )",
    )

    parser.add_argument(
        "--splat",
        type=str,
        default=None,
        help=(
            "Spectral scans will be scrunched into a single "
            "channel containing data in the given frequency "
            "range, starting from the frequency of the first"
            " bin. E.g. '0:1000' indicates 'from the first "
            "bin of the spectrum up to 1000 MHz above'. ':' "
            "or 'all' for all the channels."
        ),
    )

    args = parser.parse_args(args)

    if args.sample_config:
        sample_config_file()
        sys.exit()

    if args.bad_chans == "":
        bad_chans = []
    else:
        bad_chans = args.bad_chans.split(",")

    outfile = args.outfile

    excluded_xy, excluded_radec = _excluded_regions_from_args(args.exclude)

    if args.altaz:
        args.frame = "altaz"

    if args.file is not None:
        scanset = ScanSet(args.file, config_file=args.config, plot=not args.noplot)
        infile = args.file
        if outfile is None:
            outfile = infile
    else:
        if args.config is None:
            raise ValueError("Please specify the config file!")
        scanset = ScanSet(
            args.config,
            norefilt=not args.refilt,
            freqsplat=args.splat,
            nosub=not args.sub,
            nofilt=args.nofilt,
            debug=args.debug,
            avoid_regions=excluded_radec,
            nosave=args.nosave,
            plot=not args.noplot,
        )
        infile = args.config

        if outfile is None:
            outfile = infile.replace(".ini", "_dump.hdf5")

    if args.interactive:
        scanset.interactive_display()

    if args.global_fit:
        scanset.fit_full_images(chans=args.chans, frame=args.frame, excluded=excluded_xy)
    scanset.write(outfile, overwrite=True)

    scanset.save_ds9_images(
        save_sdev=True,
        scrunch=args.scrunch_channels,
        frame=args.frame,
        calibration=args.calibrate,
        map_unit=args.unit,
        calibrate_scans=not args.quick,
        destripe=args.destripe,
        npix_tol=args.npix_tol,
        bad_chans=bad_chans,
    )

    scanset.write(outfile, overwrite=True)


def main_preprocess(args=None):
    """Preprocess the data."""
    import argparse

    description = (
        "Load a series of scans from a config file "
        "and preprocess them, or preprocess a single scan."
    )
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "files",
        nargs="*",
        help="Single files to preprocess",
        default=None,
        type=str,
    )

    parser.add_argument("-c", "--config", type=str, default=None, help="Config file")

    parser.add_argument(
        "--sub",
        default=False,
        action="store_true",
        help="Subtract the baseline from single scans",
    )

    parser.add_argument(
        "--interactive",
        default=False,
        action="store_true",
        help="Open the interactive display for each scan",
    )

    parser.add_argument(
        "--nofilt",
        action="store_true",
        default=False,
        help="Do not filter noisy channels",
    )

    parser.add_argument("--debug", action="store_true", default=False, help="Be verbose")

    parser.add_argument("--plot", action="store_true", default=False, help="Plot stuff")

    parser.add_argument(
        "--nosave",
        action="store_true",
        default=False,
        help="Do not save the hdf5 intermediate files when" "loading subscans.",
    )

    parser.add_argument(
        "--splat",
        type=str,
        default=None,
        help=(
            "Spectral scans will be scrunched into a single "
            "channel containing data in the given frequency "
            "range, starting from the frequency of the first"
            " bin. E.g. '0:1000' indicates 'from the first "
            "bin of the spectrum up to 1000 MHz above'. ':' "
            "or 'all' for all the channels."
        ),
    )

    parser.add_argument(
        "-e",
        "--exclude",
        nargs="+",
        default=None,
        help="Exclude region from global fitting of baseline "
        "and baseline subtraction. It can be specified "
        "as X1, Y1, radius1, X2, Y2, radius2 in image "
        "coordinates or as a ds9-compatible region file "
        "in image or fk5 coordinates containing circular "
        "regions to be excluded. Currently, baseline "
        "subtraction only takes into account fk5 "
        "coordinates and global fitting image coordinates"
        ". This will change in the future.",
    )

    args = parser.parse_args(args)

    excluded_xy, excluded_radec = _excluded_regions_from_args(args.exclude)

    if args.files is not None and args.files:
        for f in args.files:
            if not os.path.exists(f):
                warnings.warn("File {} does not exist".format(f))
                continue

            kind = detect_data_kind(f)
            if kind is None:
                continue

            Scan(
                f,
                freqsplat=args.splat,
                nosub=not args.sub,
                norefilt=False,
                debug=args.debug,
                plot=args.plot,
                interactive=args.interactive,
                avoid_regions=excluded_radec,
                config_file=args.config,
                nosave=args.nosave,
            )
    else:
        if args.config is None:
            raise ValueError("Please specify the config file!")
        ScanSet(
            args.config,
            norefilt=False,
            freqsplat=args.splat,
            nosub=not args.sub,
            nofilt=args.nofilt,
            debug=args.debug,
            plot=args.plot,
            interactive=args.interactive,
            avoid_regions=excluded_radec,
            nosave=args.nosave,
        )
    return 0
