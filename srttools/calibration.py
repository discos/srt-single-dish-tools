"""Produce calibrated light curves.

``SDTlcurve`` is a script that, given a list of cross scans from different
sources, is able to recognize calibrators and use them to convert the observed
counts into a density flux value in Jy.
"""

import configparser
import copy
import glob
import json
import logging
import os
import re
import traceback
import warnings

import numpy as np
from scipy.optimize import curve_fit

import astropy.units as u
from astropy.stats import mad_std
from astropy.table import Column, Table

from .crossscan import beam_fwhm_from_frequency
from .fit import fit_baseline_plus_bell
from .io import get_channel_feed, mkdir_p
from .read_config import (
    add_directory_arguments,
    apply_directory_arguments,
    get_config_file,
    output_file_name,
    read_config,
    sample_config_file,
)
from .scan import Scan, list_scans
from .utils import (
    HAS_STATSM,
    TWOPI,
    calculate_moments,
    info_once,
    new_figure,
    scantype,
    standard_byte,
)

try:
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    HAS_MPL = True
except ImportError:
    HAS_MPL = False

CALIBRATOR_CONFIG = None


__all__ = ["CalibratorTable", "read_calibrator_config"]


# Version of the calibrator fits (2: amplitudes corrected for the offset across
# the subscan, RA err on the sky). Cached fits of other versions are recalculated.
CALIBRATOR_FIT_VERSION = 2


def _constant(x, p):
    return p


FLUX_QUANTITIES = {
    "Jy/beam": "Flux",
    "Jy/pixel": "Flux Integral",
    "Jy/sr": "Flux Integral",
}


def _get_flux_quantity(map_unit):
    try:
        return FLUX_QUANTITIES[map_unit]
    except Exception:
        raise ValueError(
            "Incorrect map_unit for flux conversion. Use one " f"of {list(FLUX_QUANTITIES.keys())}"
        )


def read_calibrator_config():
    """Read the configuration of calibrators in data/calibrators.

    Returns
    -------
    configs : dict
        Dictionary containing the configuration for each calibrator. Each key
        is the name of a calibrator. Each entry is another dictionary, in one
        of the following formats:
        1) {'Kind' : 'FreqList', 'Frequencies' : [...], 'Bandwidths' : [...],
        'Fluxes' : [...], 'Flux Errors' : [...]}
        where 'Frequencies' is the list of observing frequencies in GHz,
        'Bandwidths' is the list of bandwidths in GHz, 'Fluxes' is the list of
        flux densities in Jy from the literature and 'Flux Errors' are the
        uncertainties on those fluxes.
        2) {'Kind' : 'CoeffTable', 'CoeffTable':
        {'coeffs' : 'time, a0, a0e, a1, a1e, a2, a2e, a3, a3e\n2010.0,0 ...}}
        where the 'coeffs' key contains a dictionary with the table of
        coefficients a la Perley & Butler ApJS 204, 19 (2013), as a
        comma-separated string.

    Examples
    --------
    >>> calibs = read_calibrator_config() # doctest: +ELLIPSIS
    >>> calibs['DummyCal']['Kind']
    'CoeffTable'
    >>> 'coeffs' in calibs['DummyCal']['CoeffTable']
    True
    """
    flux_re = re.compile(r"^Flux")
    curdir = os.path.dirname(__file__)
    calibdir = os.path.join(curdir, "data", "calibrators")
    calibrator_file_list = glob.glob(os.path.join(calibdir, "*.ini"))

    configs = {}
    for cfile in calibrator_file_list:
        cparser = configparser.ConfigParser()
        cparser.read(cfile)

        info_once(logging, f"Reading {cfile}")
        if "CoeffTable" not in list(cparser.sections()):
            configs[cparser.get("Info", "Name")] = {
                "Kind": "FreqList",
                "Frequencies": [],
                "Bandwidths": [],
                "Fluxes": [],
                "Flux Errors": [],
            }

            for section in cparser.sections():
                if not flux_re.match(section):
                    continue
                configs[cparser.get("Info", "Name")]["Frequencies"].append(
                    float(cparser.get(section, "freq"))
                )
                configs[cparser.get("Info", "Name")]["Bandwidths"].append(
                    float(cparser.get(section, "bwidth"))
                )
                configs[cparser.get("Info", "Name")]["Fluxes"].append(
                    float(cparser.get(section, "flux"))
                )
                configs[cparser.get("Info", "Name")]["Flux Errors"].append(
                    float(cparser.get(section, "eflux"))
                )
        else:
            configs[cparser.get("Info", "Name")] = {
                "CoeffTable": dict(cparser.items("CoeffTable")),
                "Kind": "CoeffTable",
            }
    return configs


def _match_calibrator_name(calibrator, calibrators, relax=False):
    """Match a calibrator name to the ones in the config file.

    Parameters
    ----------
    calibrator : str
        Name of the calibrator to match.
    calibrators : list of str
        List of calibrators in the config file.

    Returns
    -------
    calibrator : str
        Name of the calibrator in the config file.

    Examples
    --------
    >>> calibrators = ['3C48', '3C286', '3C147', '3C138', '3C295', '3C196', '3C9']
    >>> calibrator = '3C48'
    >>> _match_calibrator_name(calibrator, calibrators)
    '3C48'
    >>> calibrator = '3C286a'
    >>> cal = _match_calibrator_name(calibrator, calibrators)  # This should fail
    >>> cal is None
    True
    >>> cal = _match_calibrator_name(calibrator, calibrators, relax=True)  # This should work
    >>> cal
    '3C286'
    """
    if relax:
        for cal in calibrators:
            if calibrator == cal or calibrator in cal or cal in calibrator:
                return cal
    else:
        for cal in calibrators:
            if calibrator == cal:
                return cal
    return None


def find_calibrator_in_list(calibrator, calibrators):
    """Find a calibrator in the config file.

    Start by searching for the exact name, and if not present, look for
    similar names.

    Parameters
    ----------
    calibrator : str
        Name of the calibrator to match.
    calibrators : list of str
        List of calibrators in the config file.

    Returns
    -------
    calibrator : str
        Name of the calibrator in the config file.

    Examples
    --------
    >>> calibrators = ['3C48', '3C286', '3C147', '3C138', '3C295', '3C196', '3C9']
    >>> calibrator = '3C48'
    >>> find_calibrator_in_list(calibrator, calibrators)
    '3C48'
    >>> calibrator = '3C286a'
    >>> cal = find_calibrator_in_list(calibrator, calibrators)
    >>> cal
    '3C286'
    >>> calibrator = 'bla'
    >>> cal = find_calibrator_in_list(calibrator, calibrators)
    >>> cal is None
    True
    """
    new_calibrator = _match_calibrator_name(calibrator, calibrators)

    if new_calibrator is None:
        logging.warning(
            f"Calibrator {calibrator} not found with exact name in config file. Trying to relax the search."
        )
        new_calibrator = _match_calibrator_name(calibrator, calibrators, relax=True)

        if new_calibrator is None:
            logging.warning(f"Calibrator {calibrator} not found in config file.")
        else:
            logging.warning(f"Found similarly-named {new_calibrator} in config file.")

    return new_calibrator


def _get_calibrator_flux(calibrator, frequency, bandwidth=1, time=0):
    global CALIBRATOR_CONFIG

    info_once(logging, f"Getting calibrator flux from {calibrator}")

    if CALIBRATOR_CONFIG is None:
        CALIBRATOR_CONFIG = read_calibrator_config()

    calibrators = CALIBRATOR_CONFIG.keys()
    new_calibrator_name = find_calibrator_in_list(calibrator, calibrators)
    if new_calibrator_name is None:
        return None, None
    conf = CALIBRATOR_CONFIG[new_calibrator_name]

    # find closest value among frequencies
    if conf["Kind"] == "FreqList":
        idx = (np.abs(np.array(conf["Frequencies"]) - frequency)).argmin()
        return conf["Fluxes"][idx], conf["Flux Errors"][idx]
    elif conf["Kind"] == "CoeffTable":
        return _calc_flux_from_coeffs(conf, frequency, bandwidth, time)


def _safe_fit_baseline_plus_bell(x, y, label, **kwargs):
    """Run ``fit_baseline_plus_bell``, returning ``(None, None)`` if it raises.

    The fit can raise (e.g. ``NonFiniteValueError`` when the Gaussian width
    collapses on a single spike). One bad subscan should not stop the
    processing of all the others, so this warns and lets the caller skip it.
    """
    try:
        return fit_baseline_plus_bell(x, y, kind="gauss", **kwargs)
    except Exception as e:
        warnings.warn(f"Fit failed in {label}: {type(e).__name__}: {e}")
        return None, None


def _closest_approach(ras, decs, ra0, dec0):
    """Index and distance (deg) of the point of the track (ras, decs) closest to (ra0, dec0).

    All coordinates in degrees.
    """
    dra = (ras - ra0 + 180) % 360 - 180
    dist = np.hypot(dra * np.cos(np.radians(dec0)), decs - dec0)
    i_min = np.argmin(dist)
    return i_min, dist[i_min]


def _track_length(ras, decs):
    """Angular length (deg) of the track (ras, decs), in degrees."""
    dra = (np.diff(ras) + 180) % 360 - 180
    return np.sum(np.hypot(dra * np.cos(np.radians(decs[:-1])), np.diff(decs)))


def _treat_scan(scan_path, plot=False, max_offset=0.25, **kwargs):
    """Fit a Gaussian plus baseline to each channel of a calibrator subscan.

    Channels whose feed never passes within ``max_offset`` beam FWHMs of the
    source are skipped (``None`` fits them all): in a calibrator map most
    subscans miss the source, and their fits give meaningless Flux/Counts.
    """
    scandir, sname = os.path.split(scan_path)

    try:
        # For now, use nosave. HDF5 doesn't store meta, essential for
        # this
        scan = Scan(scan_path, norefilt=True, nosave=True, plot=plot, **kwargs)
    except KeyError as e:
        logging.warning(f"Missing key. Bad file? {sname}: {str(e)}")
        return False, None
    except Exception as e:
        logging.warning(f"Error while processing {sname}: {str(e)}")
        logging.warning(traceback.format_exc())
        return False, None

    if plot and HAS_MPL:
        # Diagnostic plots go next to the processed scan, in intermediatedir
        outdir = scan.root_name(scan_path) + "_scanfit"
        mkdir_p(outdir)

    chans = scan.chan_columns()

    rows = []
    for nch, channel in enumerate(chans):
        if channel[-1:] in ["Q", "U"]:
            continue
        # Each channel only sees the sky through its own feed
        feed = get_channel_feed(channel)
        if feed is None:
            feed = 0

        ras = np.degrees(scan["ra"][:, feed])
        decs = np.degrees(scan["dec"][:, feed])
        els = np.degrees(scan["el"][:, feed])
        azs = np.degrees(scan["az"][:, feed])
        time = np.mean(scan["time"][:])
        el = np.mean(els)
        az = np.mean(azs)
        source = scan.meta["SOURCE"]
        pnt_ra = np.degrees(scan.meta["RA"])
        pnt_dec = np.degrees(scan.meta["Dec"])
        i_closest, offset = _closest_approach(
            ras, decs, pnt_ra.to_value(u.deg), pnt_dec.to_value(u.deg)
        )
        beam_fwhm = beam_fwhm_from_frequency(scan, channel)
        if max_offset is not None:
            if offset > max_offset * beam_fwhm:
                logging.debug(f"{sname}, {channel}: misses {source} by {offset * 60:.2f}'")
                continue
        frequency = scan[channel].meta["frequency"]
        bandwidth = scan[channel].meta["bandwidth"]
        temperature = scan[channel + "-Temp"]

        y = scan[channel]

        # Fit for gain curves
        x, _ = scantype(ras, decs, els, azs)
        temperature_model, _ = _safe_fit_baseline_plus_bell(
            x, temperature, f"{sname}, {channel}-Temp"
        )
        source_temperature = np.nan
        if temperature_model is not None:
            source_temperature = temperature_model["Bell"].amplitude.value

        # Fit RA and/or Dec
        x, scan_type = scantype(ras, decs)
        # Start from where the feed passes closest to the source, with the beam width
        # (in the units of x). The default start, the middle of the subscan, misses
        # sources off centre, as in most subscans of a map.
        beam_sigma = beam_fwhm / (2 * np.sqrt(2 * np.log(2)))
        x_per_deg = np.sum(np.abs(np.diff(x))) / _track_length(ras, decs)
        model, fit_info = _safe_fit_baseline_plus_bell(
            x, y, f"{sname}, {channel}", mean=x[i_closest], stddev=beam_sigma * x_per_deg
        )
        if model is None:
            continue

        std = np.std(np.diff(y)) / np.sqrt(2)
        try:
            uncert = fit_info["param_cov"].diagonal() ** 0.5
        except Exception:
            message = fit_info["message"]
            warnings.warn(f"Fit failed in {sname}, {channel}: {message}")
            continue
        bell = model["Bell"]
        baseline = model["Baseline"]
        # pars = model.parameters
        pnames = model.param_names
        counts = model.amplitude_1.value

        backsub = y - baseline(x)
        resid_std = np.std(y - model(x))
        moments = calculate_moments(backsub)
        skewness = moments["skewness"]
        kurtosis = moments["kurtosis"]

        if scan_type.startswith("RA"):
            fit_ra = bell.mean.value
            fit_width = bell.stddev.value * np.cos(np.radians(pnt_dec))
            fit_dec = None
            # Offset on the sky, like Dec err
            ra_err = (fit_ra * u.degree - pnt_ra) * np.cos(np.radians(pnt_dec))
            dec_err = None
            fit_mean = fit_ra
            fit_label = "RA"
            # Where the track passes, across the scan, from the pointed position
            track_offset = decs[i_closest] - pnt_dec.to_value(u.deg)
            pnt = pnt_ra
        elif scan_type.startswith("Dec"):
            fit_ra = None
            fit_dec = bell.mean.value
            fit_width = bell.stddev.value
            dec_err = fit_dec * u.degree - pnt_dec
            ra_err = None
            fit_mean = fit_dec
            fit_label = "Dec"
            track_offset = (ras[i_closest] - pnt_ra.to_value(u.deg)) * np.cos(np.radians(pnt_dec))
            pnt = pnt_dec
        else:
            raise ValueError("Unknown scan type")

        index = pnames.index("amplitude_1")

        counts_err = uncert[index]
        index = pnames.index("stddev_1")
        width_err = uncert[index]

        flux_density, flux_density_err = 0, 0
        flux_over_counts, flux_over_counts_err = 0, 0
        calculated_flux, calculated_flux_err = 0, 0

        new_row = [
            scandir,
            sname,
            scan_type,
            source,
            channel,
            feed,
            time,
            frequency,
            bandwidth,
            std,
            counts,
            counts_err,
            fit_width,
            width_err,
            flux_density,
            flux_density_err,
            el,
            az,
            source_temperature,
            flux_over_counts,
            flux_over_counts_err,
            flux_over_counts,
            flux_over_counts_err,
            calculated_flux,
            calculated_flux_err,
            pnt_ra,
            pnt_dec,
            fit_ra,
            fit_dec,
            ra_err,
            dec_err,
            skewness,
            kurtosis,
            resid_std,
            False,
            track_offset,
            1.0,
        ]
        rows.append(new_row)

        if plot and HAS_MPL:
            fig = new_figure()
            gs = GridSpec(2, 1, figure=fig, height_ratios=(3, 1))
            ax0 = fig.add_subplot(gs[0])
            ax1 = fig.add_subplot(gs[1], sharex=ax0)

            ax0.plot(x, y, label="Data")
            ax0.plot(
                x,
                bell(x),
                label=f"Fit: Amp: {counts}, Wid: {fit_width}",
            )
            ax1.plot(x, y - bell(x))

            ax0.axvline(fit_mean, label=fit_label + " Fit", ls="-")
            ax0.axvline(pnt.to(u.deg).value, label=fit_label + " Pnt", ls="--")
            ax0.set_xlim([min(x), max(x)])
            ax1.set_xlabel(fit_label)
            ax0.set_ylabel("Counts")
            ax1.set_ylabel("Residual (cts)")

            ax0.legend()

            fig.savefig(os.path.join(outdir, f"Feed{feed}_chan{nch:03d}.png"))
            fig = new_figure()
            gs = GridSpec(2, 1, figure=fig, height_ratios=(3, 1))
            ax0 = fig.add_subplot(gs[0])
            ax1 = fig.add_subplot(gs[1], sharex=ax0)

            ax0.plot(x, temperature, label="Data")
            if temperature_model is not None:
                ax0.plot(x, temperature_model(x), label="Fit")
                ax1.plot(x, temperature - temperature_model(x))

            ax0.axvline(pnt.to(u.deg).value, label=fit_label + " Pnt", ls="--")
            ax0.set_xlim([min(x), max(x)])
            ax1.set_xlabel(fit_label)
            ax0.set_ylabel("Counts")
            ax1.set_ylabel("Residual (cts)")

            ax0.legend()
            fig.savefig(os.path.join(outdir, f"Feed{feed}_chan{nch:03d}_temp.png"))

    return True, rows


class CalibratorTable(Table):
    """Table composed of fitted and tabulated fluxes."""

    def __init__(self, *args, **kwargs):
        """Initialize the object."""
        Table.__init__(self, *args, **kwargs)
        self.calibration_coeffs = {}
        self.calibration_uncerts = {}
        self.calibration = {}
        self.valid_elevation = {}

        names = [
            "Dir",
            "File",
            "Scan Type",
            "Source",
            "Chan",
            "Feed",
            "Time",
            "Frequency",
            "Bandwidth",
            "Data Std",
            "Counts",
            "Counts Err",
            "Width",
            "Width Err",
            "Flux",
            "Flux Err",
            "Elevation",
            "Azimuth",
            "Source_temperature",
            "Flux/Counts",
            "Flux/Counts Err",
            "Flux Integral/Counts",
            "Flux Integral/Counts Err",
            "Calculated Flux",
            "Calculated Flux Err",
            "RA",
            "Dec",
            "Fit RA",
            "Fit Dec",
            "RA err",
            "Dec err",
            "Skewness",
            "Kurtosis",
            "Resid Std",
            "Outlier",
            "Track Offset",
            "Offset Corr",
        ]

        dtype = [
            "S200",
            "S200",
            "S200",
            "S200",
            "S200",
            int,
            np.double,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            float,
            bool,
            float,
            float,
        ]

        defaults = {"S200": "", int: 0, np.double: np.nan, float: np.nan, bool: False}
        for n, d in zip(names, dtype):
            if n not in self.keys():
                # Tables written by older versions lack the newest columns
                default = 1.0 if n == "Offset Corr" else defaults[d]
                self.add_column(Column(name=n, dtype=d, data=[default] * len(self)))

    def from_scans(
        self,
        scan_list=None,
        debug=False,
        freqsplat=None,
        config_file=None,
        nofilt=False,
        plot=False,
    ):
        """Load source table from a list of scans.

        For each scan, a fit is performed. Since we are assuming point-like
        sources here, the fit is a Gaussian plus a slope. The centroid, width
        and amplitude of the fit fill out new rows of the CalibratorTable
        ('Fit RA' or 'Fit Dec', 'Width' and 'Counts' respectively).
        Channels whose feed passes farther than 1/4 of the beam FWHM from the
        source are not fitted: in a calibrator map, only the subscans crossing
        the source are used.

        Parameters
        ----------
        scan_list : list of str
            List of files containing cross scans to be fitted
        config_file : str
            File containing the configuration (list of directories etc.)

        Other Parameters
        ----------------
        debug : bool
            Throw debug information
        freqsplat : str
            List of frequencies to be merged into one. See
            :func:`srttools.scan.interpret_frequency_range`
        nofilt : bool
            Do not filter the noisy channels of the scan. See
            :class:`srttools.scan.clean_scan_using_variability`
        plot : bool
            Plot diagnostic plots? Default False, True if debug is True.

        Returns
        -------
        retval : bool
            True if at least one scan was correctly processed

        See Also
        --------
        srttools.scan.interpret_frequency_range
        """
        if debug is True:
            plot = True

        if scan_list is None:
            if config_file is None:
                config_file = get_config_file()
            config = read_config(config_file)
            scan_list = list_scans(config["datadir"], config["list_of_directories"]) + list_scans(
                config["datadir"], config["calibrator_directories"]
            )
            scan_list.sort()
        nscan = len(scan_list)

        out_retval = False
        for i_s, s in enumerate(scan_list):
            logging.info(f"{i_s + 1}/{nscan}: Loading {s}")

            retval, rows = _treat_scan(
                s,
                plot=plot,
                debug=debug,
                freqsplat=freqsplat,
                nofilt=nofilt,
                config_file=config_file,
            )

            if retval:
                out_retval = True

                for r in rows:
                    self.add_row(r)

        self.correct_pointing_offsets()
        return out_retval

    def correct_pointing_offsets(self, max_corr=2.0):
        """Correct the amplitudes for the offset of the source across each subscan.

        The fit of a subscan measures the source position along it, not
        across it: if the source is off by ``d`` across the subscan, the
        track passes beside the peak and the amplitude is lower by
        ``exp(-d**2 / (2 sigma**2))``, with ``sigma`` the beam width in that
        direction. The offset across an RA subscan is measured by the Dec
        subscans of the same observation and channel, and vice versa: we use
        the mean ``Dec err`` (or ``RA err``) and ``Width`` of the (up to) two
        perpendicular subscans closest in time, and ``Track Offset``, where
        the track passes across the scan with respect to the pointed
        position. ``Counts`` and ``Counts Err`` are multiplied by the
        correction factor, stored in ``Offset Corr``.

        Subscans without perpendicular subscans (e.g. a map scanned in one
        direction) are not corrected (``Offset Corr`` = 1); nor, with a
        warning, those whose correction would be larger than ``max_corr``:
        the source is too far off the track for a reliable amplitude.
        Calling it again recomputes the correction from the uncorrected
        amplitudes.
        """
        if len(self) == 0:
            return
        old_corr = np.array(self["Offset Corr"], dtype=float)
        old_corr[~np.isfinite(old_corr)] = 1.0
        new_corr = np.ones(len(self))
        is_ra = np.array([str(s).startswith("RA") for s in self["Scan Type"]])
        times = np.array(self["Time"], dtype=float)
        width = np.array(self["Width"], dtype=float)
        track = np.array(self["Track Offset"], dtype=float)
        # Source offset measured along each subscan, on the sky
        along = np.where(
            is_ra, np.array(self["RA err"], dtype=float), np.array(self["Dec err"], dtype=float)
        )
        keys = np.array([f"{d}::{c}" for d, c in zip(self["Dir"], self["Chan"])])
        for key in np.unique(keys):
            group = np.flatnonzero(keys == key)
            for i in group:
                perp = group[(is_ra[group] != is_ra[i]) & np.isfinite(along[group])]
                perp = perp[np.isfinite(width[perp]) & (width[perp] > 0)]
                if len(perp) == 0 or not np.isfinite(track[i]):
                    continue
                closest = perp[np.argsort(np.abs(times[perp] - times[i]))[:2]]
                across = track[i] - np.mean(along[closest])
                corr = np.exp(across**2 / (2 * np.mean(width[closest]) ** 2))
                if corr > max_corr:
                    warnings.warn(
                        f"{self['File'][i]}, {self['Chan'][i]}: the source is "
                        f"{across * 60:.2f}' across the subscan; amplitude not corrected"
                    )
                    continue
                new_corr[i] = corr
        factor = new_corr / old_corr
        self["Counts"] *= factor
        self["Counts Err"] *= factor
        self["Offset Corr"] = new_corr

    def write(self, fname, *args, **kwargs):
        """Same as Table.write, but adds path information for HDF5."""
        if fname.endswith(".hdf5"):
            super().write(fname, *args, **kwargs)
        else:
            super().write(fname, *args, **kwargs)

    def check_not_empty(self):
        """Check that table is not empty.

        Returns
        -------
        good : bool
            True if all checks pass, False otherwise.
        """
        if len(self["Flux/Counts"]) == 0:
            warnings.warn("The calibrator table is empty!")
            return False
        return True

    def check_up_to_date(self):
        """Check that the calibration information is up to date.

        Returns
        -------
        good : bool
            True if all checks pass, False otherwise.
        """
        if not self.check_not_empty():
            return False

        if np.any(self["Flux/Counts"] == 0):
            warnings.warn("The calibrator table needs an update!")
            self.update()

        return True

    def update(self, outlier_nsigma=5.0):
        """Update the calibration information.

        Execute ``get_fluxes``, ``calibrate``, ``flag_outliers`` (unless
        ``outlier_nsigma`` is None) and ``compute_conversion_function``
        """
        if not self.check_not_empty():
            return

        self.get_fluxes()
        self.calibrate()
        if outlier_nsigma is not None:
            self.flag_outliers(nsigma=outlier_nsigma)
        self.compute_conversion_function()

    def _valid_mask(self, good_mask=None):
        """Rows usable for calibration: ``good_mask`` (default: calibrators) and not outliers."""
        if good_mask is None:
            good_mask = self["Flux"] > 0
        return np.asarray(good_mask) & ~np.asarray(self["Outlier"], dtype=bool)

    def flag_outliers(self, nsigma=5.0, min_fits=6):
        """Flag fits that are outliers with respect to similar fits, e.g. because of RFI.

        Fits are compared within groups with the same observation (``Dir``),
        source, channel and arm (RA or Dec): within a single observation the
        source is not expected to change. A fit is flagged in the
        ``Outlier`` column if any of the following is more than ``nsigma``
        robust standard deviations (1.4826 times the median absolute
        deviation) away from the group median: the rms of the fit residuals
        ("Resid Std", only excesses), the amplitude ("Counts") or the width.
        The last two catch fits that locked on to an RFI burst instead of the
        source. Groups with fewer than ``min_fits`` fits are not checked.

        Other Parameters
        ----------------
        nsigma : float
            Rejection threshold
        min_fits : int
            Minimum number of fits in a group

        Returns
        -------
        outliers : array of bool
            The new ``Outlier`` column
        """
        outliers = np.zeros(len(self), dtype=bool)
        arms = np.array([str(s).strip()[:2].upper() for s in self["Scan Type"]])
        keys = [np.asarray(self[c]).astype(str) for c in ["Dir", "Source", "Chan"]]
        groups = {}
        for i, key in enumerate(zip(*keys, arms)):
            groups.setdefault(key, []).append(i)

        def robust_z(values):
            values = np.asarray(values, dtype=float)
            if not np.any(np.isfinite(values)):
                # E.g. "Resid Std" in tables written by older versions
                return np.zeros_like(values)
            med = np.nanmedian(values)
            spread = 1.4826 * np.nanmedian(np.abs(values - med))
            spread = max(spread, 1e-6 * np.abs(med), np.finfo(float).tiny)
            return (values - med) / spread

        for key, idx in groups.items():
            if len(idx) < min_fits:
                continue
            idx = np.array(idx)
            z_resid = robust_z(self["Resid Std"][idx])
            bad = z_resid > nsigma
            for col in ["Counts", "Width"]:
                bad |= np.abs(robust_z(self[col][idx])) > nsigma
            outliers[idx[bad]] = True
            if np.any(bad):
                files = (
                    list(np.asarray(self["File"][idx[bad]]).astype(str))
                    if "File" in self.colnames
                    else []
                )
                warnings.warn(
                    f"{np.count_nonzero(bad)}/{len(idx)} fits flagged as outliers "
                    f"(e.g. RFI) in {key[0]}, {key[1]}, {key[2]}, {key[3]} scans: "
                    f"{files}. They are excluded from calibration and fluxes. "
                    "Use outlier_nsigma=None (--outlier-nsigma 0) to disable."
                )
                if np.count_nonzero(bad) > 0.2 * len(idx):
                    warnings.warn(
                        f"More than 20% of the fits flagged in {key[0]}, {key[1]}, "
                        f"{key[2]}. Check the data."
                    )
        self["Outlier"][:] = outliers
        return outliers

    def get_fluxes(self):
        """Get the tabulated flux of the source, if listed as calibrators.

        Updates the table.
        """
        if not self.check_not_empty():
            return

        for it, t in enumerate(self["Time"]):
            source = self["Source"][it]
            frequency = self["Frequency"][it] / 1000
            bandwidth = self["Bandwidth"][it] / 1000
            flux, eflux = _get_calibrator_flux(source, frequency, bandwidth, time=t)

            self["Flux"][it] = flux
            self["Flux Err"][it] = eflux

    def calibrate(self):
        """Calculate the calibration constants.

        The following conversion functions are calculated for each tabulated
        cross scan belonging to a calibrator:

        + 'Flux/Counts' and 'Flux/Counts Err': Tabulated flux density divided
          by the _height_ of the fitted Gaussian. This is used, e.g. to
          calibrate images in Jy/beam, as it calibrates the local amplitude to
          the flux density

        + 'Flux Integral/Counts' and 'Flux Integral/Counts Err': Tabulated flux
          density divided by the _volume_ of the 2D Gaussian corresponding to
          the fitted cross scans, assuming a symmetrical beam (which is
          generally not the case, but a good approximation). This is used,
          e.g., to perform the calibration in Jy/pixel: Each pixel will be
          normalized to the expected total flux in the corresponding pixel
          area

        See Also
        --------
        srttools.calibration.CalibratorTable.from_scans
        """
        if not self.check_not_empty():
            return

        flux = self["Flux"] * u.Jy
        eflux = self["Flux Err"] * u.Jy
        counts = self["Counts"] * u.ct
        ecounts = self["Counts Err"] * u.ct
        width = np.radians(self["Width"]) * u.radian
        ewidth = np.radians(self["Width Err"]) * u.radian

        # Volume in a beam: For a 2-d Gaussian with amplitude A and sigmas sx
        # and sy, this is 2 pi A sx sy.
        total = TWOPI * counts * width**2
        etotal = TWOPI * ecounts * width**2

        flux_integral_over_counts = flux / total
        flux_integral_over_counts_err = (
            etotal / total + eflux / flux + 2 * ewidth / width
        ) * flux_integral_over_counts

        flux_over_counts = flux / counts
        flux_over_counts_err = (ecounts / counts + eflux / flux) * flux_over_counts

        self["Flux/Counts"][:] = flux_over_counts.to(u.Jy / u.ct).value
        self["Flux/Counts Err"][:] = flux_over_counts_err.to(u.Jy / u.ct).value

        self["Flux Integral/Counts"][:] = flux_integral_over_counts.to(
            u.Jy / u.ct / u.steradian
        ).value
        self["Flux Integral/Counts Err"][:] = flux_integral_over_counts_err.to(
            u.Jy / u.ct / u.steradian
        ).value

    def compute_conversion_function(self, map_unit="Jy/beam", good_mask=None):
        """Compute the conversion between Jy and counts.

        Try to get a meaningful second-degree polynomial fit over elevation.
        Revert to the rough function :func:`Jy_over_counts_rough` in case
        ``statsmodels`` is not installed. In this latter case, only the baseline
        value is given for flux conversion and error.
        These values are saved in the ``calibration_coeffs`` and
        ``calibration_uncerts`` attributes of ``CalibratorTable``, and a
        dictionary called ``calibration`` is also created. For each channel,
        this dictionary contains either None or an object. This object is the
        output of a ``fit`` procedure in ``statsmodels``. The method
        object.predict(X) returns the calibration corresponding to elevation X.
        """
        if not HAS_STATSM:
            channels = list(set(self["Chan"]))
            for channel in channels:
                fc, fce = self.Jy_over_counts_rough(
                    channel=channel, map_unit=map_unit, good_mask=None
                )
                self.calibration_coeffs[channel] = [fc, 0, 0]
                self.calibration_uncerts[channel] = [fce, 0, 0]
                self.calibration[channel] = None
            return
        else:
            import statsmodels.api as sm

        good_mask = self._valid_mask(good_mask)

        flux_quantity = _get_flux_quantity(map_unit)

        channels = list(set(self["Chan"]))
        for channel in channels:
            good_chans = (self["Chan"] == channel) & good_mask
            # N.B. doing this after checking the list of channels in the table, which
            # might be some numpy string type
            channel = str(channel)

            f_c_ratio = self[flux_quantity + "/Counts"][good_chans]
            f_c_ratio_err = self[flux_quantity + "/Counts Err"][good_chans]
            elvs = np.radians(self["Elevation"][good_chans])

            # isfinite: fits with zero amplitude give an infinite ratio
            good_fc = np.isfinite(f_c_ratio) & (f_c_ratio > 0)
            good_fce = np.isfinite(f_c_ratio_err) & (f_c_ratio_err >= 0)

            good = good_fc & good_fce

            x_to_fit = np.array(elvs[good])
            y_to_fit = np.array(f_c_ratio[good])
            ye_to_fit = np.array(f_c_ratio_err[good])

            order = np.argsort(x_to_fit)
            x_to_fit = x_to_fit[order]
            y_to_fit = y_to_fit[order]
            ye_to_fit = ye_to_fit[order]

            if len(x_to_fit) < 3:
                # Too few points for a line in elevation: constant calibration
                fc, fce = self.Jy_over_counts_rough(
                    channel=channel, map_unit=map_unit, good_mask=good_mask
                )
                self.calibration_coeffs[channel] = [fc, 0, 0]
                self.calibration_uncerts[channel] = [fce, 0, 0]
                self.calibration[channel] = None
                continue

            X = np.column_stack((np.ones(len(x_to_fit)), x_to_fit))
            # X = np.c_[np.ones(len(x_to_fit)), X]
            # X = sm.add_constant(X)
            model = sm.RLM(y_to_fit, X, missing="drop")
            results = model.fit()
            self.calibration_coeffs[channel] = results.params
            self.calibration_uncerts[channel] = results.cov_params().diagonal() ** 0.5
            self.calibration[channel] = results
            self.valid_elevation[channel] = [np.min(x_to_fit), np.max(x_to_fit)]

    def Jy_over_counts(self, channel=None, elevation=None, map_unit="Jy/beam", good_mask=None):
        """Compute the Jy/Counts conversion corresponding to a given map unit.

        Parameters
        ----------
        channel : str
            Channel name (e.g. 'Feed0_RCP', 'Feed0_LCP' etc.)
        elevation : float or array-like
            The elevation or a list of elevations
        map_unit : str
            A valid unit for the calibrated map (See the keys of
            FLUX_QUANTITIES)
        good_mask : array of bools, default None
            This mask can be used to specify the valid entries of the table.
            If None, the mask is set to an array of True values

        Returns
        -------
        fc : float or array-like
            One conversion value for each elevation
        fce : float or array-like
            the uncertainties corresponding to each ``fc``
        """
        channel = str(channel)
        use_rough = False
        if not HAS_STATSM:
            warnings.warn("No statsmodels found.")
            use_rough = True

        good_mask = self._valid_mask(good_mask)

        flux_quantity = _get_flux_quantity(map_unit)
        if channel not in self["Chan"]:
            warnings.warn(f"No calibration found for channel {channel}")
            return None, None

        if channel not in self.calibration.keys():
            self.compute_conversion_function(map_unit, good_mask=good_mask)
        if elevation is None or self.valid_elevation == {}:
            warnings.warn("No elevation given.")
            use_rough = True
        elif channel is None:
            warnings.warn("No channel given.")
            use_rough = True
        elif self.valid_elevation == {}:
            warnings.warn(f"No elevation info found in the calibration for channel {channel}")
            use_rough = True
        elif self.calibration.get(channel) is None:
            # Too few calibrator fits for an elevation dependence
            use_rough = True
        else:
            elevation_range = self.valid_elevation[channel]
            elevation_span = elevation_range[1] - elevation_range[0]
            tolerated_el_min, tolerated_el_max = (
                elevation_range[0] - 0.5 * elevation_span,
                elevation_range[1] + 0.5 * elevation_span,
            )

            elevation_out_of_range = np.any(elevation < tolerated_el_min) or np.any(
                elevation > tolerated_el_max
            )
            if elevation_out_of_range:
                warnings.warn("Some values of elevation are too far from the calibrated ones.")
                use_rough = True

        if use_rough:
            logging.info("Using rough calibration (with no elevation dependence).")
            elevation = np.array(elevation)
            fc, fce = self.Jy_over_counts_rough(
                channel=channel, map_unit=map_unit, good_mask=good_mask
            )
            if elevation.size > 1:
                fc = np.zeros_like(elevation) + fc
                fce = np.zeros_like(elevation) + fce
            return fc, fce

        X = np.column_stack((np.ones(np.array(elevation).size), np.array(elevation)))
        fc = self.calibration[channel].predict(X)
        goodch = self["Chan"] == channel
        good = good_mask & goodch

        fce = np.sqrt(np.mean(self[flux_quantity + "/Counts Err"][good] ** 2)) + np.zeros_like(fc)

        if len(fc) == 1:
            fc, fce = fc[0], fce[0]

        return fc, fce

    def Jy_over_counts_rough(self, channel=None, map_unit="Jy/beam", good_mask=None):
        """Get the conversion from counts to Jy.

        Other Parameters
        ----------------
        channel : str
            Name of the data channel
        map_unit : str
            A valid unit for the calibrated map (See the keys of
            FLUX_QUANTITIES)
        good_mask : array of bools, default None
            This mask can be used to specify the valid entries of the table.
            If None, the mask is set to an array of True values

        Returns
        -------
        fc : float
            flux density /count ratio
        fce : float
            uncertainty on ``fc``: never smaller than the robust standard
            deviation of all the fits (including those rejected as outliers)
            divided by the square root of their number
        """
        self.check_up_to_date()

        flux_quantity = _get_flux_quantity(map_unit)

        good_mask = self._valid_mask(good_mask)

        good_chans = np.ones(len(self["Time"]), dtype=bool)
        if channel is not None:
            good_chans = self["Chan"] == channel

        good_chans = good_chans & good_mask

        f_c_ratio = self[flux_quantity + "/Counts"][good_chans]
        f_c_ratio_err = self[flux_quantity + "/Counts Err"][good_chans]
        times = self["Time"][good_chans]

        # isfinite: fits with zero amplitude give an infinite ratio
        good_fc = np.isfinite(f_c_ratio) & (f_c_ratio > 0)
        good_fce = np.isfinite(f_c_ratio_err) & (f_c_ratio_err >= 0)

        good = good_fc & good_fce

        x_to_fit = np.array(times[good])
        y_to_fit = np.array(f_c_ratio[good])
        ye_to_fit = np.array(f_c_ratio_err[good])

        p = [np.median(y_to_fit)]
        pcov = np.array([[np.median(ye_to_fit) ** 2]])
        first = True
        # Error of the mean from the robust scatter of all the fits, before
        # rejecting outliers: the fits that survive the rejection can agree
        # among themselves much better than the calibration is known
        n_all = len(y_to_fit)
        scatter_err = mad_std(y_to_fit) / np.sqrt(n_all) if n_all > 1 else 0.0

        while 1:
            bad = np.abs((y_to_fit - _constant(x_to_fit, p)) / ye_to_fit) > 5

            if not np.any(bad) and not first:
                break

            if len(x_to_fit[bad]) > len(x_to_fit) - 5:
                warnings.warn("Calibration fit is shaky")
                # Too few points, or points scattered well beyond their
                # formal errors: median, with the error of the median from
                # the robust scatter of the points (never smaller than the
                # formal errors allow)
                n = len(y_to_fit)
                scatter = max(mad_std(y_to_fit), np.median(ye_to_fit))
                p = [np.median(y_to_fit)]
                pcov = np.array([[(1.2533 * scatter / np.sqrt(n)) ** 2]])
                break

            xbad = x_to_fit[bad]
            ybad = y_to_fit[bad]
            for xb, yb in zip(xbad, ybad):
                logging.warning(f"Outliers: {xb}, {yb}")

            good = np.logical_not(bad)
            x_to_fit = x_to_fit[good]
            y_to_fit = y_to_fit[good]
            ye_to_fit = ye_to_fit[good]

            p, pcov = curve_fit(_constant, x_to_fit, y_to_fit, sigma=ye_to_fit, p0=p)
            first = False
        fc = p[0]
        fce = max(np.sqrt(pcov[0, 0]), scatter_err)

        return fc, fce

    def calculate_src_flux(self, channel=None, map_unit="Jy/beam", source=None):
        """Calculate source flux and error, pointing by pointing.

        Uses the conversion factors calculated from the tabulated fluxes for
        all sources but the current, and the fitted Gaussian amplitude for the
        current source.
        Updates the calibrator table and returns the average flux

        Parameters
        ----------
        channel : str or list of str
            Data channel
        map_unit : str
            Units in the map (default Jy/beam)
        source : str
            Source name. Must match one of the sources in the table.
            Default

        Returns
        -------
        mean_flux : array of floats
            Array with as many channels as the input ones
        mean_flux_err : array of floats
            Uncertainties corresponding to mean_flux
        """
        if source is None:
            good_source = np.ones_like(self["Flux"], dtype=bool)
        else:
            good_source = self["Source"] == source

        non_source = np.logical_not(good_source)
        # Outliers are excluded from the calibration by Jy_over_counts
        good_source = good_source & ~np.asarray(self["Outlier"], dtype=bool)

        if channel is None:
            channels = list(set(self["Chan"]))
        else:
            channels = [channel]

        mean_flux = []
        mean_flux_err = []
        for ch in channels:
            good_chan = self["Chan"] == ch
            good = good_source & good_chan
            elevation = np.radians(self["Elevation"][good])
            fc, fce = self.Jy_over_counts(
                channel=ch,
                elevation=elevation,
                map_unit=map_unit,
                good_mask=non_source,
            )

            calculated_flux = copy.deepcopy(self["Calculated Flux"])
            calculated_flux_err = copy.deepcopy(self["Calculated Flux Err"])
            counts = np.array(self["Counts"])
            counts_err = np.array(self["Counts Err"])

            calculated_flux[good] = counts[good] * fc
            calculated_flux_err[good] = (
                counts_err[good] / counts[good] + fce / fc
            ) * calculated_flux[good]

            self["Calculated Flux"][:] = calculated_flux
            self["Calculated Flux Err"][:] = calculated_flux_err

            mean_flux.append(np.mean(calculated_flux[good]))
            mean_flux_err.append(np.sqrt(np.mean(calculated_flux_err[good] ** 2)))

        return mean_flux, mean_flux_err

    def check_consistency(self, channel=None, epsilon=0.05):
        """Check the consistency of calculated and fitted flux densities.

        For each source in the ``srttools``' calibrator list, use
        ``calculate_src_flux`` to calculate the source flux ignoring the
        tabulated value, and compare the calculated and tabulated values.

        Returns
        -------
        retval : bool
            True if, for all calibrators, the tabulated and calculated values
            of the flux are consistent. False otherwise.
        """
        is_cal = (~np.isnan(self["Flux"])) & (self["Flux"] > 0)
        calibrators = list(set(self["Source"][is_cal]))
        for cal in calibrators:
            self.calculate_src_flux(channel=channel, source=cal)

        if channel is None:
            good_chan = np.ones_like(self["Chan"], dtype=bool)
        else:
            good_chan = self["Chan"] == channel

        calc_fluxes = self["Calculated Flux"][is_cal & good_chan]
        biblio_fluxes = self["Flux"][is_cal & good_chan]
        names = self["Source"][is_cal & good_chan]
        times = self["Time"][is_cal & good_chan]

        consistent = np.abs(biblio_fluxes - calc_fluxes) < epsilon * biblio_fluxes

        for (
            n,
            t,
            b,
            c,
            cons,
        ) in zip(names, times, biblio_fluxes, calc_fluxes, consistent):
            if not cons:
                warnings.warn(f"{n}, MJD {t}: Expected {b}, " f"measured {c}")

        return consistent

    def beam_width(self, channel=None):
        """Calculate the (weighted) mean beam width, in radians.

        Checks for invalid (nan and such) values.
        """
        goodch = np.ones(len(self), dtype=bool)
        if channel is not None:
            goodch = self["Chan"] == channel
        allwidths = self[goodch]["Width"]
        allwidth_errs = self[goodch]["Width Err"]
        good = (allwidth_errs > 0) & (allwidth_errs == allwidth_errs)
        allwidths = allwidths[good]
        allwidth_errs = allwidth_errs[good]

        # Weighted mean
        width = np.sum(allwidths / allwidth_errs) / np.sum(1 / allwidth_errs)

        width_err = np.sqrt(np.sum(allwidth_errs**2))
        return np.radians(width), np.radians(width_err)

    def counts_over_Jy(self, channel=None, elevation=None):
        """Get the conversion from Jy to counts."""
        self.check_up_to_date()

        fc, fce = self.Jy_over_counts(channel=channel, elevation=elevation)
        cf = 1 / fc
        return cf, fce / fc * cf

    def plot_two_columns(
        self,
        xcol,
        ycol,
        xerrcol=None,
        yerrcol=None,
        ax=None,
        channel=None,
        xfactor=1,
        yfactor=1,
        color=None,
        test=False,
    ):
        """Plot the data corresponding to two given columns."""
        showit = False
        if ax is None:
            # Shown in a window, unless testing
            showit = not test
            ax = new_figure(interactive=showit).add_subplot()

        good = (self[xcol] == self[xcol]) & (self[ycol] == self[ycol])
        mask = np.ones_like(good)
        label = ""
        if channel is not None:
            mask = self["Chan"] == channel
            label = f"_{channel}"

        good = good & mask
        x_to_plot = np.array(self[xcol][good]) * xfactor
        order = np.argsort(x_to_plot)
        y_to_plot = np.array(self[ycol][good]) * yfactor
        y_to_plot = y_to_plot[order]
        yerr_to_plot = None
        xerr_to_plot = None
        if xerrcol is not None:
            xerr_to_plot = np.array(self[xerrcol][good]) * xfactor
            xerr_to_plot = xerr_to_plot[order]
        if yerrcol is not None:
            yerr_to_plot = np.array(self[yerrcol][good]) * yfactor
            yerr_to_plot = yerr_to_plot[order]

        if xerrcol is not None or yerrcol is not None:
            ax.errorbar(
                x_to_plot,
                y_to_plot,
                xerr=xerr_to_plot,
                yerr=yerr_to_plot,
                label=ycol + label,
                fmt="none",
                color=color,
                ecolor=color,
            )
        else:
            ax.scatter(x_to_plot, y_to_plot, label=ycol + label, color=color)

        if showit:
            plt.show()
        return x_to_plot, y_to_plot

    def show(self, filename=None):
        """Save a summary plot of the calibration (default ``calibration_summary.png``)."""
        from matplotlib import cm

        # TODO: this is meant to become interactive. I will make different
        # panels linked to each other.
        if filename is None:
            filename = "calibration_summary.png"
        fig = new_figure(figsize=(16, 16))
        fig.suptitle("Summary")
        gs = GridSpec(2, 2, figure=fig, hspace=0)
        ax00 = fig.add_subplot(gs[0, 0])
        ax01 = fig.add_subplot(gs[0, 1], sharey=ax00)
        ax10 = fig.add_subplot(gs[1, 0], sharex=ax00)
        ax11 = fig.add_subplot(gs[1, 1], sharex=ax01, sharey=ax10)

        channels = list(set(self["Chan"]))
        colors = cm.rainbow(np.linspace(0, 1, len(channels)))
        for ic, channel in enumerate(channels):
            # Ugly workaround for python 2-3 compatibility
            channel_str = channel

            color = colors[ic]
            self.plot_two_columns(
                "Elevation",
                "Flux/Counts",
                yerrcol="Flux/Counts Err",
                ax=ax00,
                channel=channel,
                color=color,
            )

            elevations = np.arange(np.min(self["Elevation"]), np.max(self["Elevation"]), 0.001)
            jy_over_cts, jy_over_cts_err = self.Jy_over_counts(channel_str, np.radians(elevations))
            ax00.plot(elevations, jy_over_cts, color=color)
            ax00.plot(elevations, jy_over_cts + jy_over_cts_err, color=color)
            ax00.plot(elevations, jy_over_cts - jy_over_cts_err, color=color)
            self.plot_two_columns(
                "Elevation",
                "RA err",
                ax=ax10,
                channel=channel,
                yfactor=60,
                color=color,
            )
            self.plot_two_columns(
                "Elevation",
                "Dec err",
                ax=ax10,
                channel=channel,
                yfactor=60,
                color=color,
            )
            self.plot_two_columns(
                "Azimuth",
                "Flux/Counts",
                yerrcol="Flux/Counts Err",
                ax=ax01,
                channel=channel,
                color=color,
            )

            jy_over_cts, jy_over_cts_err = self.Jy_over_counts(
                channel_str, np.radians(np.mean(elevations))
            )

            ax01.axhline(jy_over_cts, color=color)
            ax01.axhline(jy_over_cts + jy_over_cts_err, color=color)
            ax01.axhline(jy_over_cts - jy_over_cts_err, color=color)
            self.plot_two_columns(
                "Azimuth",
                "RA err",
                ax=ax11,
                channel=channel,
                yfactor=60,
                color=color,
            )
            self.plot_two_columns(
                "Azimuth",
                "Dec err",
                ax=ax11,
                channel=channel,
                yfactor=60,
                color=color,
            )

        for i in np.arange(-1, 1, 0.1):
            # Arcmin errors
            ax10.axhline(i, ls="--", color="gray")
            ax11.axhline(i, ls="--", color="gray")
        #            ax11.text(1, i, "{}".format())
        ax00.legend()
        ax01.legend()
        ax10.legend()
        ax11.legend()
        ax10.set_xlabel("Elevation")
        ax11.set_xlabel("Azimuth")
        ax00.set_ylabel("Flux / Counts")
        ax10.set_ylabel("Pointing error (arcmin)")
        fig.savefig(filename)


def flux_function(start_frequency, bandwidth, coeffs, ecoeffs):
    """Flux function from Perley & Butler ApJS 204, 19 (2013) (PB13).

    Parameters
    ----------
    start_frequency : float
        Starting frequency of the data, in GHz
    bandwidth : float
        Bandwidth, in GHz
    coeffs : list of floats
        Parameters of the PB13 interpolation
    ecoeffs : list of floats
        Uncertainties of the PB13 interpolation
    """
    a0, a1, a2, a3 = coeffs

    if np.all(ecoeffs < 1e10):
        # assume 5% error on calibration parameters!
        ecoeffs = coeffs * 0.05
    a0e, a1e, a2e, a3e = ecoeffs
    f0 = start_frequency
    f1 = start_frequency + bandwidth

    fs = np.linspace(f0, f1, 21)
    df = np.diff(fs)[0]
    fmean = (fs[:-1] + fs[1:]) / 2

    logf = np.log10(fmean)
    logS = a0 + a1 * logf + a2 * logf**2 + a3 * logf**3
    elogS = a0e + a1e * logf + a2e * logf**2 + a3e * logf**3

    S = 10**logS
    eS = S * elogS

    # Error is not random, should add linearly; divide by bandwidth
    return np.sum(S) * df / bandwidth, np.sum(eS) * df / bandwidth


def _calc_flux_from_coeffs(conf, frequency, bandwidth=1, time=0):
    """Return the flux of a calibrator at a given frequency.

    Uses Perley & Butler ApJS 204, 19 (2013).
    """
    import io

    coefftable = conf["CoeffTable"]["coeffs"]
    fobj = io.BytesIO(standard_byte(coefftable))
    table = Table.read(fobj, format="ascii.csv")

    idx = np.argmin(np.abs(np.longdouble(table["time"]) - time))

    a0, a0e = table["a0", "a0e"][idx]
    a1, a1e = table["a1", "a1e"][idx]
    a2, a2e = table["a2", "a2e"][idx]
    a3, a3e = table["a3", "a3e"][idx]
    coeffs = np.array([a0, a1, a2, a3], dtype=float)

    ecoeffs = np.array([a0e, a1e, a2e, a3e], dtype=float)

    return flux_function(frequency, bandwidth, coeffs, ecoeffs)


def main_cal(args=None):
    """Main function."""
    import argparse

    description = "Load a series of cross scans from a config file " "and use them as calibrators."
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "file",
        nargs="?",
        help="Input calibration file",
        default=None,
        type=str,
    )
    parser.add_argument(
        "--sample-config",
        action="store_true",
        default=False,
        help="Produce sample config file",
    )

    parser.add_argument(
        "--nofilt",
        action="store_true",
        default=False,
        help="Do not filter noisy channels",
    )

    parser.add_argument("-c", "--config", type=str, default=None, help="Config file")
    add_directory_arguments(parser)

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
        "-o",
        "--output",
        type=str,
        default=None,
        help="Output file containing the calibration",
    )

    parser.add_argument(
        "--snr-min",
        type=float,
        default=10,
        help="Minimum SNR for calibrator measurements " "to be considered valid",
    )

    parser.add_argument(
        "--show",
        action="store_true",
        default=False,
        help="Show calibration summary",
    )

    parser.add_argument(
        "--check",
        action="store_true",
        default=False,
        help="Check consistency of calibration",
    )

    parser.add_argument(
        "--outlier-nsigma",
        type=float,
        default=5.0,
        help="Flag fits whose residual rms, amplitude or width differ from those "
        "of the other subscans of the same observation, source, channel and arm "
        "by more than this many robust standard deviations (e.g. RFI). 0 "
        "disables the rejection",
    )

    args = parser.parse_args(args)
    # Before loading any scan: remembered by read_config for the whole run
    apply_directory_arguments(args)

    if args.sample_config:
        sample_config_file()
        return

    if args.file is not None:
        caltable = CalibratorTable().read(args.file)
        caltable.show(filename=args.file.replace(".hdf5", ".jpg"))
        return

    if args.config is None:
        raise ValueError("Please specify the config file!")

    config = read_config(args.config)

    calibrator_dirs = config["calibrator_directories"]
    if calibrator_dirs is None or not calibrator_dirs:
        raise ValueError("No calibrators specified in config file")

    scan_list = list_scans(config["datadir"], config["calibrator_directories"])

    scan_list.sort()

    outfile = args.output
    if outfile is None:
        outfile = output_file_name(args.config, "_cal.hdf5")
    # A cache, reused by later runs with the same calibrator scans and options
    outfile_unfilt = output_file_name(args.config, "_cal_unfilt.hdf5", intermediate=True)
    settings = json.dumps(
        {
            "scans": [os.path.abspath(s) for s in scan_list],
            "splat": args.splat,
            "nofilt": bool(args.nofilt),
            # Increase when the fits change, to recalculate older caches
            "fit_version": CALIBRATOR_FIT_VERSION,
        }
    )
    caltable = None
    if os.path.exists(outfile_unfilt):
        caltable = CalibratorTable.read(outfile_unfilt)
        if caltable.meta.get("sdtcal_settings") == settings:
            logging.info(
                f"Loading unfiltered calibration table from {outfile_unfilt} (delete to reprocess)"
            )
        else:
            warnings.warn(
                "The calibrator scans, the options (--splat, --nofilt) or the version of the "
                "fits changed since the "
                f"cached calibration table {outfile_unfilt} was made: recalculating it."
            )
            caltable = None
    if caltable is None:
        caltable = CalibratorTable()
        caltable.from_scans(
            scan_list,
            freqsplat=args.splat,
            nofilt=args.nofilt,
            plot=args.show,
            config_file=args.config,
        )
        caltable.meta["sdtcal_settings"] = settings
        caltable.write(outfile_unfilt, overwrite=True)

    snr = caltable["Counts"] / caltable["Data Std"]
    N = len(caltable)
    good = snr > args.snr_min
    chans = list(set(caltable["Chan"]))

    caltable = caltable[good]
    logging.info(
        f"{len(caltable)} good calibrator observations found above " f"SNR={args.snr_min} (of {N})"
    )
    for chan in chans:
        good_chan = caltable["Chan"] == chan
        if not np.any(good_chan):
            warnings.warn(
                f"No good data for channel {chan}. Try using the --snr-min option with some value lower than {args.snr_min}"
            )
    outlier_nsigma = args.outlier_nsigma if args.outlier_nsigma > 0 else None
    caltable.update(outlier_nsigma=outlier_nsigma)

    if args.check:
        for chan in list(set(caltable["Chan"])):
            caltable.check_consistency(chan)
    if args.show:
        caltable.show(filename=outfile.replace(".hdf5", ".jpg"))

    caltable.write(outfile, overwrite=True)
    caltable.write(outfile.replace(".hdf5", ".csv"), overwrite=True)


def main_lcurve(args=None):
    """Main function."""
    import argparse

    description = (
        "Load a series of cross scans from a config file " "and obtain a calibrated curve."
    )
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "file",
        nargs="?",
        help="Input calibration file",
        default=None,
        type=str,
    )
    parser.add_argument(
        "-s",
        "--source",
        nargs="+",
        type=str,
        default=None,
        help="Source or list of sources",
    )
    parser.add_argument(
        "--sample-config",
        action="store_true",
        default=False,
        help="Produce sample config file",
    )

    parser.add_argument(
        "--nofilt",
        action="store_true",
        default=False,
        help="Do not filter noisy channels",
    )

    parser.add_argument("-c", "--config", type=str, default=None, help="Config file")
    add_directory_arguments(parser)

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
        "-o",
        "--output",
        type=str,
        default=None,
        help="Output file containing the calibration",
    )

    parser.add_argument(
        "--outlier-nsigma",
        type=float,
        default=5.0,
        help="Flag fits whose residual rms, amplitude or width differ from those "
        "of the other subscans of the same observation, source, channel and arm "
        "by more than this many robust standard deviations (e.g. RFI). 0 "
        "disables the rejection",
    )

    args = parser.parse_args(args)
    # Before loading any scan: remembered by read_config for the whole run
    apply_directory_arguments(args)

    if args.sample_config:
        sample_config_file()
        return

    outlier_nsigma = args.outlier_nsigma if args.outlier_nsigma > 0 else None
    if args.file is not None:
        caltable = CalibratorTable.read(args.file)
        caltable.update(outlier_nsigma=outlier_nsigma)
    else:
        if args.config is None:
            raise ValueError("Please specify the config file!")
        caltable = CalibratorTable()
        caltable.from_scans(config_file=args.config)
        caltable.update(outlier_nsigma=outlier_nsigma)

        outfile = args.output
        if outfile is None:
            outfile = output_file_name(args.config, "_cal.hdf5")

        caltable.write(outfile, overwrite=True)

    sources = args.source
    if args.source is None:
        sources = list(set(caltable["Source"]))

    for s in sources:
        caltable.calculate_src_flux(source=s)
        good = (caltable["Source"] == s) & ~np.asarray(caltable["Outlier"], dtype=bool)
        lctable = Table()
        lctable.add_column(Column(name="Time", dtype=float))
        lctable["Time"] = caltable["Time"][good]
        lctable["Flux"] = caltable["Calculated Flux"][good]
        lctable["Flux Err"] = caltable["Calculated Flux Err"][good]
        lctable["Chan"] = caltable["Chan"][good]
        lctable.write(s.replace(" ", "_") + ".csv", overwrite=True)
