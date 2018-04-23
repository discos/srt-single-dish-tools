"""
Produce calibrated light curves.

``SDTlcurve`` is a script that, given a list of cross scans from different
sources, is able to recognize calibrators and use them to convert the observed
counts into a density flux value in Jy.


"""
from __future__ import (absolute_import, division,
                        print_function)

from .scan import Scan, list_scans
from .read_config import read_config, sample_config_file, get_config_file
from .fit import fit_baseline_plus_bell
from .io import mkdir_p
from .utils import standard_string, standard_byte, compare_strings
from .utils import HAS_STATSM, calculate_moments, scantype

import os
import sys
import glob
import re
import warnings
import traceback
from scipy.optimize import curve_fit
import logging
import astropy.units as u

import numpy as np
from astropy.table import Table, Column

# For Python 2 and 3 compatibility
try:
    import configparser
except ImportError:
    import ConfigParser as configparser

try:
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

CALIBRATOR_CONFIG = None


__all__ = ["CalibratorTable", "read_calibrator_config"]


def _constant(x, p):
    return p


FLUX_QUANTITIES = {"Jy/beam": "Flux",
                   "Jy/pixel": "Flux Integral",
                   "Jy/sr": "Flux Integral"}


def _get_flux_quantity(map_unit):
    try:
        return FLUX_QUANTITIES[map_unit]
    except Exception:
        raise ValueError("Incorrect map_unit for flux conversion. Use one "
                         "of {}".format(list(FLUX_QUANTITIES.keys())))


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

    See Also
    --------
    srttools.calibration.flux_function

    Examples
    --------
    >>> calibs = read_calibrator_config()
    >>> calibs['DummyCal']['Kind']
    'CoeffTable'
    >>> 'coeffs' in calibs['DummyCal']['CoeffTable']
    True
    """
    flux_re = re.compile(r'^Flux')
    curdir = os.path.dirname(__file__)
    calibdir = os.path.join(curdir, 'data', 'calibrators')
    calibrator_file_list = glob.glob(os.path.join(calibdir, '*.ini'))

    configs = {}
    for cfile in calibrator_file_list:
        cparser = configparser.ConfigParser()
        cparser.read(cfile)

        if 'CoeffTable' not in list(cparser.sections()):
            configs[cparser.get("Info", "Name")] = {"Kind": "FreqList",
                                                    "Frequencies": [],
                                                    "Bandwidths": [],
                                                    "Fluxes": [],
                                                    "Flux Errors": []}

            for section in cparser.sections():
                if not flux_re.match(section):
                    continue
                configs[cparser.get("Info", "Name")]["Frequencies"].append(
                    float(cparser.get(section, "freq")))
                configs[cparser.get("Info", "Name")]["Bandwidths"].append(
                    float(cparser.get(section, "bwidth")))
                configs[cparser.get("Info", "Name")]["Fluxes"].append(
                    float(cparser.get(section, "flux")))
                configs[cparser.get("Info", "Name")]["Flux Errors"].append(
                    float(cparser.get(section, "eflux")))
        else:
            configs[cparser.get("Info", "Name")] = \
                {"CoeffTable": dict(cparser.items("CoeffTable")),
                 "Kind": "CoeffTable"}

    return configs


def _get_calibrator_flux(calibrator, frequency, bandwidth=1, time=0):
    global CALIBRATOR_CONFIG

    if CALIBRATOR_CONFIG is None:
        CALIBRATOR_CONFIG = read_calibrator_config()

    calibrators = CALIBRATOR_CONFIG.keys()

    for cal in calibrators:
        if cal in calibrator:
            calibrator = cal
            break
    else:
        return None, None

    conf = CALIBRATOR_CONFIG[calibrator]
    # find closest value among frequencies
    if conf["Kind"] == "FreqList":
        idx = (np.abs(np.array(conf["Frequencies"]) - frequency)).argmin()
        return conf["Fluxes"][idx] * bandwidth, \
            conf["Flux Errors"][idx] * bandwidth
    elif conf["Kind"] == "CoeffTable":
        return _calc_flux_from_coeffs(conf, frequency, bandwidth, time)


def _treat_scan(scan_path, plot=False, **kwargs):
    scandir, sname = os.path.split(scan_path)
    if plot and HAS_MPL:
        outdir = os.path.splitext(sname)[0] + "_scanfit"
        outdir = os.path.join(scandir, outdir)
        mkdir_p(outdir)

    try:
        # For now, use nosave. HDF5 doesn't store meta, essential for
        # this
        scan = Scan(scan_path, norefilt=True, nosave=True, **kwargs)
    except KeyError as e:
        logging.warning("Missing key. Bad file? {}: {}".format(sname,
                                                               str(e)))
        return False, None
    except Exception as e:
        logging.warning("Error while processing {}: {}".format(sname,
                                                               str(e)))
        warnings.warn(traceback.format_exc())
        return False, None

    feeds = np.arange(scan['ra'].shape[1])
    chans = scan.chan_columns()

    chan_nums = np.arange(len(chans))
    F, N = np.meshgrid(feeds, chan_nums)
    F = F.flatten()
    N = N.flatten()
    rows = []
    for feed, nch in zip(F, N):
        channel = chans[nch]

        ras = np.degrees(scan['ra'][:, feed])
        decs = np.degrees(scan['dec'][:, feed])
        els = np.degrees(scan['el'][:, feed])
        azs = np.degrees(scan['az'][:, feed])
        time = np.mean(scan['time'][:])
        el = np.mean(els)
        az = np.mean(azs)
        source = scan.meta['SOURCE']
        pnt_ra = np.degrees(scan.meta['RA'])
        pnt_dec = np.degrees(scan.meta['Dec'])
        frequency = scan[channel].meta['frequency']
        bandwidth = scan[channel].meta['bandwidth']
        temperature = scan[channel + '-Temp']

        y = scan[channel]

        # Fit for gain curves
        x, _ = scantype(ras, decs, els, azs)
        temperature_model, _ = \
            fit_baseline_plus_bell(x, temperature, kind='gauss')
        source_temperature = temperature_model['Bell'].amplitude.value

        # Fit RA and/or Dec
        x, scan_type = scantype(ras, decs)
        model, fit_info = fit_baseline_plus_bell(x, y, kind='gauss')

        try:
            uncert = fit_info['param_cov'].diagonal() ** 0.5
        except Exception:
            message = fit_info['message']
            warnings.warn(
                "Fit failed in scan {s}: {m}".format(s=sname,
                                                     m=message))
            continue
        bell = model['Bell']
        baseline = model['Baseline']
        # pars = model.parameters
        pnames = model.param_names
        counts = model.amplitude_1.value

        backsub = y - baseline(x)
        moments = calculate_moments(backsub)
        skewness = moments['skewness']
        kurtosis = moments['kurtosis']

        if scan_type.startswith("RA"):
            fit_ra = bell.mean.value
            fit_width = bell.stddev.value * np.cos(np.radians(pnt_dec))
            fit_dec = None
            ra_err = fit_ra * u.degree - pnt_ra
            dec_err = None
            fit_mean = fit_ra
            fit_label = 'RA'
            pnt = pnt_ra
        elif scan_type.startswith("Dec"):
            fit_ra = None
            fit_dec = bell.mean.value
            fit_width = bell.stddev.value
            dec_err = fit_dec * u.degree - pnt_dec
            ra_err = None
            fit_mean = fit_dec
            fit_label = 'Dec'
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

        rows.append([scandir, sname, scan_type, source, channel, feed,
                     time, frequency, bandwidth, counts, counts_err,
                     fit_width, width_err,
                     flux_density, flux_density_err, el, az,
                     source_temperature,
                     flux_over_counts, flux_over_counts_err,
                     flux_over_counts, flux_over_counts_err,
                     calculated_flux, calculated_flux_err,
                     pnt_ra, pnt_dec, fit_ra, fit_dec, ra_err,
                     dec_err, skewness, kurtosis])

        if plot and HAS_MPL:
            fig = plt.figure("Fit information")
            gs = GridSpec(2, 1, height_ratios=(3, 1))
            ax0 = plt.subplot(gs[0])
            ax1 = plt.subplot(gs[1], sharex=ax0)

            ax0.plot(x, y, label="Data")
            ax0.plot(x, bell(x),
                     label="Fit: Amp: {}, Wid: {}".format(counts, fit_width))
            ax1.plot(x, y - bell(x))

            ax0.axvline(fit_mean, label=fit_label + " Fit", ls="-")
            ax0.axvline(pnt.to(u.deg).value, label=fit_label + " Pnt",
                        ls="--")
            ax0.set_xlim([min(x), max(x)])
            ax1.set_xlabel(fit_label)
            ax0.set_ylabel("Counts")
            ax1.set_ylabel("Residual (cts)")

            ax0.legend()
            ax1.legend()
            plt.savefig(os.path.join(outdir,
                                     "Feed{}_chan{}.png".format(feed,
                                                                nch)))
            plt.close(fig)
            fig = plt.figure("Fit information - temperature")
            gs = GridSpec(2, 1, height_ratios=(3, 1))
            ax0 = plt.subplot(gs[0])
            ax1 = plt.subplot(gs[1], sharex=ax0)

            ax0.plot(x, temperature, label="Data")
            ax0.plot(x, temperature_model['Bell'](x), label="Fit")
            ax1.plot(x, temperature - temperature_model['Bell'](x))

            ax0.axvline(pnt.to(u.deg).value, label=fit_label + " Pnt",
                        ls="--")
            ax0.set_xlim([min(x), max(x)])
            ax1.set_xlabel(fit_label)
            ax0.set_ylabel("Counts")
            ax1.set_ylabel("Residual (cts)")

            plt.legend()
            plt.savefig(os.path.join(outdir,
                                     "Feed{}_chan{}_temp.png".format(feed,
                                                                     nch)))
            plt.close(fig)

    return True, rows


class CalibratorTable(Table):
    """Table composed of fitted and tabulated fluxes."""

    def __init__(self, *args, **kwargs):
        """Initialize the object."""
        Table.__init__(self, *args, **kwargs)
        self.calibration_coeffs = {}
        self.calibration_uncerts = {}
        self.calibration = {}
        names = ["Dir", "File", "Scan Type", "Source",
                 "Chan", "Feed", "Time",
                 "Frequency", "Bandwidth",
                 "Counts", "Counts Err",
                 "Width", "Width Err",
                 "Flux", "Flux Err",
                 "Elevation", "Azimuth",
                 "Source_temperature",
                 "Flux/Counts", "Flux/Counts Err",
                 "Flux Integral/Counts", "Flux Integral/Counts Err",
                 "Calculated Flux", "Calculated Flux Err",
                 "RA", "Dec",
                 "Fit RA", "Fit Dec",
                 "RA err", "Dec err",
                 "Skewness", "Kurtosis"]

        dtype = ['S200', 'S200', 'S200', 'S200',
                 'S200', np.int, np.double,
                 np.float, np.float,
                 np.float, np.float,
                 np.float, np.float,
                 np.float, np.float, np.float,
                 np.float, np.float,
                 np.float, np.float,
                 np.float, np.float,
                 np.float, np.float,
                 np.float, np.float,
                 np.float, np.float,
                 np.float, np.float,
                 np.float, np.float]

        for n, d in zip(names, dtype):
            if n not in self.keys():
                self.add_column(Column(name=n, dtype=d))

    def from_scans(self, scan_list=None, debug=False, freqsplat=None,
                   config_file=None, nofilt=False, plot=False):
        """Load source table from a list of scans.

        For each scan, a fit is performed. Since we are assuming point-like
        sources here, the fit is a Gaussian plus a slope. The centroid, width
        and amplitude of the fit fill out new rows of the CalibratorTable
        ('Fit RA' or 'Fit Dec', 'Width' and 'Counts' respectively).

        Parameters
        ----------
        scan_list : list of str
            List of files containing cross scans to be fitted
        config_file : str
            File containing the configuration (list of directories etc.)

        Other parameters
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
            scan_list = \
                list_scans(config['datadir'],
                           config['list_of_directories']) + \
                list_scans(config['datadir'],
                           config['calibrator_directories'])
            scan_list.sort()
        nscan = len(scan_list)

        out_retval = False
        for i_s, s in enumerate(scan_list):
            logging.info('{}/{}: Loading {}'.format(i_s + 1, nscan, s))

            retval, rows = _treat_scan(s, plot=plot, debug=debug,
                                       freqsplat=freqsplat, nofilt=nofilt)

            if retval:
                out_retval = True

                for r in rows:
                    self.add_row(r)

        return out_retval

    def write(self, fname, *args, **kwargs):
        """Same as Table.write, but adds path information for HDF5."""
        if fname.endswith('.hdf5'):
            Table.write(self, fname, *args, path='table', **kwargs)
        else:
            Table.write(self, fname, *args, **kwargs)

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

    def update(self):
        """Update the calibration information.

        Execute ``get_fluxes``, ``calibrate`` and
        ``compute_conversion_function``
        """
        if not self.check_not_empty():
            return

        self.get_fluxes()
        self.calibrate()
        self.compute_conversion_function()

    def get_fluxes(self):
        """Get the tabulated flux of the source, if listed as calibrators.

        Updates the table.
        """
        if not self.check_not_empty():
            return

        for it, t in enumerate(self['Time']):
            source = standard_string(self['Source'][it])
            frequency = self['Frequency'][it] / 1000
            bandwidth = self['Bandwidth'][it] / 1000
            flux, eflux = \
                _get_calibrator_flux(source, frequency, bandwidth, time=t)

            self['Flux'][it] = flux
            self['Flux Err'][it] = eflux

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

        flux = self['Flux'] * u.Jy
        eflux = self['Flux Err'] * u.Jy
        counts = self['Counts'] * u.ct
        ecounts = self['Counts Err'] * u.ct
        width = np.radians(self['Width']) * u.radian
        ewidth = np.radians(self['Width Err']) * u.radian

        # Volume in a beam: For a 2-d Gaussian with amplitude A and sigmas sx
        # and sy, this is 2 pi A sx sy.
        total = 2 * np.pi * counts * width ** 2
        etotal = 2 * np.pi * ecounts * width ** 2

        flux_integral_over_counts = flux / total
        flux_integral_over_counts_err = \
            (etotal / total + eflux / flux +
             2 * ewidth / width) * flux_integral_over_counts

        flux_over_counts = flux / counts
        flux_over_counts_err = \
            (ecounts / counts + eflux / flux) * flux_over_counts

        self['Flux/Counts'][:] = \
            flux_over_counts.to(u.Jy / u.ct).value
        self['Flux/Counts Err'][:] = \
            flux_over_counts_err.to(u.Jy / u.ct).value

        self['Flux Integral/Counts'][:] = \
            flux_integral_over_counts.to(u.Jy / u.ct / u.steradian).value
        self['Flux Integral/Counts Err'][:] = \
            flux_integral_over_counts_err.to(u.Jy / u.ct / u.steradian).value

    def compute_conversion_function(self, map_unit="Jy/beam", good_mask=None):
        """Compute the conversion between Jy and counts.

        Try to get a meaningful second-degree polynomial fit over elevation.
        Revert to the rough function `Jy_over_counts_rough` in case
        `statsmodels` is not installed. In this latter case, only the baseline
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
                fc, fce = self.Jy_over_counts_rough(channel=channel,
                                                    map_unit=map_unit,
                                                    good_mask=None)
                self.calibration_coeffs[standard_string(channel)] = [fc, 0, 0]
                self.calibration_uncerts[standard_string(channel)] = \
                    [fce, 0, 0]
                self.calibration[standard_string(channel)] = None
            return
        else:
            import statsmodels.api as sm

        if good_mask is None:
            good_mask = self['Flux'] > 0

        flux_quantity = _get_flux_quantity(map_unit)

        channels = list(set(self["Chan"]))
        for channel in channels:
            good_chans = compare_strings(self["Chan"], channel) & good_mask

            f_c_ratio = self[flux_quantity + "/Counts"][good_chans]
            f_c_ratio_err = self[flux_quantity + "/Counts Err"][good_chans]
            elvs = np.radians(self["Elevation"][good_chans])

            good_fc = (f_c_ratio == f_c_ratio) & (f_c_ratio > 0)
            good_fce = (f_c_ratio_err == f_c_ratio_err) & (f_c_ratio_err >= 0)

            good = good_fc & good_fce

            x_to_fit = np.array(elvs[good])
            y_to_fit = np.array(f_c_ratio[good])
            ye_to_fit = np.array(f_c_ratio_err[good])

            order = np.argsort(x_to_fit)
            x_to_fit = x_to_fit[order]
            y_to_fit = y_to_fit[order]
            ye_to_fit = ye_to_fit[order]

            X = np.column_stack((np.ones(len(x_to_fit)), x_to_fit))
            # X = np.c_[np.ones(len(x_to_fit)), X]
            # X = sm.add_constant(X)
            model = sm.RLM(y_to_fit, X, missing='drop')
            results = model.fit()

            self.calibration_coeffs[standard_string(channel)] = results.params
            self.calibration_uncerts[standard_string(channel)] = \
                results.cov_params().diagonal()**0.5
            self.calibration[standard_string(channel)] = results

    def Jy_over_counts(self, channel=None, elevation=None,
                       map_unit="Jy/beam", good_mask=None):
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
        rough = False
        if not HAS_STATSM:
            rough = True

        if good_mask is None:
            good_mask = self['Flux'] > 0

        flux_quantity = _get_flux_quantity(map_unit)

        if standard_string(channel) not in self.calibration.keys():
            self.compute_conversion_function(map_unit, good_mask=good_mask)

        if elevation is None or rough is True or channel is None:
            elevation = np.array(elevation)
            fc, fce = self.Jy_over_counts_rough(channel=channel,
                                                map_unit=map_unit,
                                                good_mask=good_mask)
            if elevation.size > 1:
                fc = np.zeros_like(elevation) + fc
                fce = np.zeros_like(elevation) + fce
            return fc, fce

        X = np.column_stack((np.ones(np.array(elevation).size),
                             np.array(elevation)))

        fc = self.calibration[standard_string(channel)].predict(X)

        goodch = compare_strings(self["Chan"], channel)
        good = good_mask & goodch
        fce = np.sqrt(np.mean(
            self[flux_quantity + "/Counts Err"][good]**2)) + np.zeros_like(fc)

        if len(fc) == 1:
            fc, fce = fc[0], fce[0]

        return fc, fce

    def Jy_over_counts_rough(self, channel=None, map_unit="Jy/beam",
                             good_mask=None):
        """Get the conversion from counts to Jy.

        Other parameters
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
            uncertainty on `fc`
        """

        self.check_up_to_date()

        if good_mask is None:
            good_mask = np.ones(len(self["Time"]), dtype=bool)

        good_chans = np.ones(len(self["Time"]), dtype=bool)
        if channel is not None:
            good_chans = compare_strings(self['Chan'], channel)

        good_chans = good_chans & good_mask

        flux_quantity = _get_flux_quantity(map_unit)

        f_c_ratio = self[flux_quantity + "/Counts"][good_chans]
        f_c_ratio_err = self[flux_quantity + "/Counts Err"][good_chans]
        times = self["Time"][good_chans]

        good_fc = (f_c_ratio == f_c_ratio) & (f_c_ratio > 0)
        good_fce = (f_c_ratio_err == f_c_ratio_err) & (f_c_ratio_err >= 0)

        good = good_fc & good_fce

        x_to_fit = np.array(times[good])
        y_to_fit = np.array(f_c_ratio[good])
        ye_to_fit = np.array(f_c_ratio_err[good])

        p = [np.median(y_to_fit)]
        first = True

        while 1:
            bad = np.abs((y_to_fit - _constant(x_to_fit, p)) / ye_to_fit) > 5
            if not np.any(bad) and not first:
                break

            if len(x_to_fit[bad]) > len(x_to_fit) - 5:
                warnings.warn("Calibration fit is shaky")
                break

            xbad = x_to_fit[bad]
            ybad = y_to_fit[bad]
            for xb, yb in zip(xbad, ybad):
                logging.warning("Outliers: {}, {}".format(xb, yb))

            good = np.logical_not(bad)
            x_to_fit = x_to_fit[good]
            y_to_fit = y_to_fit[good]
            ye_to_fit = ye_to_fit[good]

            p, pcov = curve_fit(_constant, x_to_fit, y_to_fit, sigma=ye_to_fit,
                                p0=p)
            first = False
        fc = p[0]
        fce = np.sqrt(pcov[0, 0])

        return fc, fce

    def calculate_src_flux(self, channel=None,
                           map_unit="Jy/beam", source=None):
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
            good_source = np.ones_like(self['Flux'], dtype=bool)
        else:
            good_source = compare_strings(self['Source'], source)
        non_source = np.logical_not(good_source)

        if channel is None:
            channels = [standard_string(s) for s in set(self['Chan'])]
        else:
            channels = [standard_string(channel)]

        mean_flux = []
        mean_flux_err = []
        for channel in channels:
            good_chan = compare_strings(self['Chan'], channel)
            good = good_source & good_chan
            elevation = np.radians(self['Elevation'][good])
            fc, fce = self.Jy_over_counts(channel=channel, elevation=elevation,
                                          map_unit=map_unit,
                                          good_mask=non_source)

            calculated_flux = np.array(self['Calculated Flux'])
            calculated_flux_err = np.array(self['Calculated Flux Err'])
            counts = np.array(self['Counts'])
            counts_err = np.array(self['Counts Err'])

            calculated_flux[good] = counts[good] * fc
            calculated_flux_err[good] = \
                (counts_err[good] / counts[good] + fce / fc) * \
                calculated_flux[good]

            self['Calculated Flux'][:] = calculated_flux
            self['Calculated Flux Err'][:] = calculated_flux_err

            mean_flux.append(np.mean(calculated_flux[good]))
            mean_flux_err.append(
                np.sqrt(np.mean(calculated_flux_err[good] ** 2)))

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
        is_cal = self['Flux'] > 0
        calibrators = list(set(self['Source'][is_cal]))
        for cal in calibrators:
            self.calculate_src_flux(channel=channel, source=cal)
        if channel is None:
            good_chan = np.ones_like(self['Chan'], dtype=bool)
        else:
            good_chan = compare_strings(self['Chan'], channel)
        calc_fluxes = self['Calculated Flux'][is_cal & good_chan]
        biblio_fluxes = self['Flux'][is_cal & good_chan]
        names = self['Source'][is_cal & good_chan]
        times = self['Time'][is_cal & good_chan]

        consistent = \
            np.abs(biblio_fluxes - calc_fluxes) < epsilon * biblio_fluxes

        for n, t, b, c, in zip(names, times, biblio_fluxes, calc_fluxes):
            consistent = np.abs(b - c) < epsilon * b
            if not consistent:
                warnings.warn("{}, MJD {}: Expected {}, "
                              "measured {}".format(n, t, b, c))

        return consistent

    def beam_width(self, channel=None):
        """Calculate the (weighted) mean beam width, in radians.

        Checks for invalid (nan and such) values.
        """
        goodch = np.ones(len(self), dtype=bool)
        if channel is not None:
            goodch = compare_strings(self['Chan'], channel)
        allwidths = self[goodch]['Width']
        allwidth_errs = self[goodch]['Width Err']
        good = (allwidth_errs > 0) & (allwidth_errs == allwidth_errs)
        allwidths = allwidths[good]
        allwidth_errs = allwidth_errs[good]

        # Weighted mean
        width = np.sum(allwidths/allwidth_errs) / np.sum(1/allwidth_errs)

        width_err = np.sqrt(np.sum(allwidth_errs ** 2))
        return np.radians(width), np.radians(width_err)

    def counts_over_Jy(self, channel=None, elevation=None):
        """Get the conversion from Jy to counts."""
        self.check_up_to_date()

        fc, fce = self.Jy_over_counts(channel=channel, elevation=elevation)
        cf = 1 / fc
        return cf, fce / fc * cf

    def plot_two_columns(self, xcol, ycol, xerrcol=None, yerrcol=None, ax=None,
                         channel=None,
                         xfactor=1, yfactor=1, color=None, test=False):
        """Plot the data corresponding to two given columns."""
        showit = False
        if ax is None:
            plt.figure("{} vs {}".format(xcol, ycol))
            ax = plt.gca()
            showit = True

        good = (self[xcol] == self[xcol]) & (self[ycol] == self[ycol])
        mask = np.ones_like(good)
        label = ""
        if channel is not None:
            mask = self['Chan'] == channel
            label = "_{}".format(channel)

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
            ax.errorbar(x_to_plot, y_to_plot,
                        xerr=xerr_to_plot,
                        yerr=yerr_to_plot,
                        label=ycol + label,
                        fmt="none", color=color,
                        ecolor=color)
        else:
            ax.scatter(x_to_plot, y_to_plot, label=ycol + label,
                       color=color)

        if showit and not test:
            plt.show()
        return x_to_plot, y_to_plot

    def show(self):
        """Show a summary of the calibration."""

        from matplotlib import cm
        # TODO: this is meant to become interactive. I will make different
        # panels linked to each other.

        fig = plt.figure("Summary", figsize=(16, 16))
        plt.suptitle("Summary")
        gs = GridSpec(2, 2, hspace=0)
        ax00 = plt.subplot(gs[0, 0])
        ax01 = plt.subplot(gs[0, 1], sharey=ax00)
        ax10 = plt.subplot(gs[1, 0], sharex=ax00)
        ax11 = plt.subplot(gs[1, 1], sharex=ax01, sharey=ax10)

        channels = list(set(self['Chan']))
        colors = cm.rainbow(np.linspace(0, 1, len(channels)))
        for ic, channel in enumerate(channels):
            # Ugly workaround for python 2-3 compatibility
            channel_str = standard_string(channel)

            color = colors[ic]
            self.plot_two_columns('Elevation', "Flux/Counts",
                                  yerrcol="Flux/Counts Err", ax=ax00,
                                  channel=channel, color=color)

            elevations = np.arange(np.min(self['Elevation']),
                                   np.max(self['Elevation']), 0.001)
            jy_over_cts, jy_over_cts_err = \
                self.Jy_over_counts(channel_str, np.radians(elevations))
            ax00.plot(elevations, jy_over_cts, color=color)
            ax00.plot(elevations, jy_over_cts + jy_over_cts_err, color=color)
            ax00.plot(elevations, jy_over_cts - jy_over_cts_err, color=color)
            self.plot_two_columns('Elevation', "RA err", ax=ax10,
                                  channel=channel,
                                  yfactor=60, color=color)
            self.plot_two_columns('Elevation', "Dec err", ax=ax10,
                                  channel=channel,
                                  yfactor=60, color=color)
            self.plot_two_columns('Azimuth', "Flux/Counts",
                                  yerrcol="Flux/Counts Err", ax=ax01,
                                  channel=channel, color=color)

            jy_over_cts, jy_over_cts_err = \
                self.Jy_over_counts(channel_str,
                                    np.radians(np.mean(elevations)))

            ax01.axhline(jy_over_cts, color=color)
            ax01.axhline(jy_over_cts + jy_over_cts_err, color=color)
            ax01.axhline(jy_over_cts - jy_over_cts_err, color=color)
            self.plot_two_columns('Azimuth', "RA err", ax=ax11,
                                  channel=channel,
                                  yfactor=60, color=color)
            self.plot_two_columns('Azimuth', "Dec err", ax=ax11,
                                  channel=channel,
                                  yfactor=60, color=color)

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
        plt.savefig("calibration_summary.png")
        plt.close(fig)


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
    fmean = (fs[:-1] + fs[1:])/2

    logf = np.log10(fmean)
    logS = a0 + a1 * logf + a2 * logf**2 + a3 * logf**3
    elogS = a0e + a1e * logf + a2e * logf**2 + a3e * logf**3

    S = 10 ** logS
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
    table = Table.read(fobj, format='ascii.csv')

    idx = np.argmin(np.abs(np.longdouble(table["time"]) - time))

    a0, a0e = table['a0', 'a0e'][idx]
    a1, a1e = table['a1', 'a1e'][idx]
    a2, a2e = table['a2', 'a2e'][idx]
    a3, a3e = table['a3', 'a3e'][idx]
    coeffs = np.array([a0, a1, a2, a3], dtype=float)

    ecoeffs = np.array([a0e, a1e, a2e, a3e], dtype=float)

    return flux_function(frequency, bandwidth, coeffs, ecoeffs)


def main_cal(args=None):
    """Main function."""
    import argparse

    description = ('Load a series of cross scans from a config file '
                   'and use them as calibrators.')
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("file", nargs='?', help="Input calibration file",
                        default=None, type=str)
    parser.add_argument("--sample-config", action='store_true', default=False,
                        help='Produce sample config file')

    parser.add_argument("--nofilt", action='store_true', default=False,
                        help='Do not filter noisy channels')

    parser.add_argument("-c", "--config", type=str, default=None,
                        help='Config file')

    parser.add_argument("--splat", type=str, default=None,
                        help=("Spectral scans will be scrunched into a single "
                              "channel containing data in the given frequency "
                              "range, starting from the frequency of the first"
                              " bin. E.g. '0:1000' indicates 'from the first "
                              "bin of the spectrum up to 1000 MHz above'. ':' "
                              "or 'all' for all the channels."))

    parser.add_argument("-o", "--output", type=str, default=None,
                        help='Output file containing the calibration')

    parser.add_argument("--show", action='store_true', default=False,
                        help='Show calibration summary')

    parser.add_argument("--check", action='store_true', default=False,
                        help='Check consistency of calibration')

    args = parser.parse_args(args)

    if args.sample_config:
        sample_config_file()
        sys.exit()

    if args.file is not None:
        caltable = CalibratorTable().read(args.file)
        caltable.show()
        sys.exit()

    if args.config is None:
        raise ValueError("Please specify the config file!")

    config = read_config(args.config)

    calibrator_dirs = config['calibrator_directories']
    if calibrator_dirs is None or not calibrator_dirs:
        raise ValueError("No calibrators specified in config file")

    scan_list = \
        list_scans(config['datadir'],
                   config['calibrator_directories'])

    scan_list.sort()

    outfile = args.output
    if outfile is None:
        outfile = args.config.replace(".ini", "_cal.hdf5")
    caltable = CalibratorTable()
    caltable.from_scans(scan_list, freqsplat=args.splat, nofilt=args.nofilt,
                        plot=args.show)
    caltable.update()

    if args.check:
        for chan in list(set(caltable['Chan'])):
            caltable.check_consistency(chan)
    if args.show:
        caltable.show()

    caltable.write(outfile, overwrite=True)
    caltable.write(outfile.replace('.hdf5', '.csv'), overwrite=True)


def main_lcurve(args=None):
    """Main function."""
    import argparse

    description = ('Load a series of cross scans from a config file '
                   'and obtain a calibrated curve.')
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("file", nargs='?', help="Input calibration file",
                        default=None, type=str)
    parser.add_argument("-s", "--source", nargs='+', type=str, default=None,
                        help='Source or list of sources')
    parser.add_argument("--sample-config", action='store_true', default=False,
                        help='Produce sample config file')

    parser.add_argument("--nofilt", action='store_true', default=False,
                        help='Do not filter noisy channels')

    parser.add_argument("-c", "--config", type=str, default=None,
                        help='Config file')

    parser.add_argument("--splat", type=str, default=None,
                        help=("Spectral scans will be scrunched into a single "
                              "channel containing data in the given frequency "
                              "range, starting from the frequency of the first"
                              " bin. E.g. '0:1000' indicates 'from the first "
                              "bin of the spectrum up to 1000 MHz above'. ':' "
                              "or 'all' for all the channels."))

    parser.add_argument("-o", "--output", type=str, default=None,
                        help='Output file containing the calibration')

    args = parser.parse_args(args)

    if args.sample_config:
        sample_config_file()
        sys.exit()

    if args.file is not None:
        caltable = CalibratorTable.read(args.file)
        caltable.update()
    else:
        if args.config is None:
            raise ValueError("Please specify the config file!")
        caltable = CalibratorTable()
        caltable.from_scans(config_file=args.config)
        caltable.update()

        outfile = args.output
        if outfile is None:
            outfile = args.config.replace(".ini", "_cal.hdf5")

        caltable.write(outfile, overwrite=True)

    sources = args.source
    if args.source is None:
        sources = [standard_string(s) for s in set(caltable['Source'])]

    for s in sources:
        caltable.calculate_src_flux(source=s)
        good = compare_strings(caltable['Source'], s)
        lctable = Table()
        lctable.add_column(Column(name="Time", dtype=float))
        lctable['Time'] = caltable['Time'][good]
        lctable['Flux'] = caltable['Calculated Flux'][good]
        lctable['Flux Err'] = caltable['Calculated Flux Err'][good]
        lctable['Chan'] = caltable['Chan'][good]
        lctable.write(s.replace(' ', '_') + '.csv', overwrite=True)
