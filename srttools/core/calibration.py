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
import os
import sys
import glob
import re
import warnings
import traceback
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import logging
import astropy.units as u

try:
    import cPickle as pickle
except:
    import pickle

import numpy as np
from astropy.table import Table, vstack, Column
# For Python 2 and 3 compatibility
try:
    import configparser
except ImportError:
    import ConfigParser as configparser

CALIBRATOR_CONFIG = None


def _calibration_function(x, pars):
    return pars[0] + pars[1] * x + pars[2] * x**2


def _constant(x, p):
    return p


def _get_flux_quantity(map_unit):
    if map_unit == "Jy/beam":
        return "Flux"
    elif map_unit in ["Jy/pixel", "Jy/sr"]:
        return "Flux Integral"
    else:
        raise ValueError("Incorrect map_unit for flux conversion. Use one of"
                         "Jy/beam, Jy/pixel")


def _scantype(ras, decs):
    """Get if scan is along RA or Dec, and if forward or backward."""
    ravar = np.max(ras) - np.min(ras)
    decvar = np.max(decs) - np.min(decs)
    if ravar > decvar:
        x = ras
        xvariab = 'RA'
    else:
        x = decs
        xvariab = 'Dec'

    if x[-1] > x[0]:
        scan_direction = '>'
    else:
        scan_direction = '<'

    return x, xvariab + scan_direction


def read_calibrator_config():
    """Read the configuration of calibrators in data/calibrators."""
    flux_re = re.compile(r'^Flux')
    curdir = os.path.dirname(__file__)
    calibdir = os.path.join(curdir, '..', 'data', 'calibrators')
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


class SourceTable(Table):
    """Class containing all information and functions about sources."""

    def __init__(self, *args, **kwargs):
        """Initialize the object."""
        Table.__init__(self, *args, **kwargs)

        names = ["Dir", "File", "Scan Type", "Source",
                 "Chan", "Feed", "Time",
                 "Frequency", "Bandwidth",
                 "Counts", "Counts Err",
                 "Width", "Width Err",
                 "Flux", "Flux Err",
                 "Elevation", "Azimuth",
                 "Flux/Counts", "Flux/Counts Err",
                 "Flux Integral/Counts", "Flux Integral/Counts Err",
                 "RA", "Dec",
                 "Fit RA", "Fit Dec",
                 "RA err", "Dec err"]

        dtype = ['S200', 'S200', 'S200', 'S200',
                 'S200', np.int, np.double,
                 np.float, np.float,
                 np.float, np.float,
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
        """Load source table from a list of scans."""

        if debug is True:
            plot = True

        if scan_list is None:
            if config_file is None:
                config_file = get_config_file()
            config = read_config(config_file)
            scan_list = \
                list_scans(config['datadir'], config['list_of_directories'])
            scan_list.sort()
        nscan = len(scan_list)

        for i_s, s in enumerate(scan_list):
            logging.info('{}/{}: Loading {}'.format(i_s + 1, nscan, s))
            scandir, sname = os.path.split(s)
            if plot:
                outdir = os.path.splitext(sname)[0] + "_scanfit"
                outdir = os.path.join(scandir, outdir)
                mkdir_p(outdir)

            try:
                # For now, use nosave. HDF5 doesn't store meta, essential for
                # this
                # TODO: experiment with serialize_meta!
                scan = Scan(s, norefilt=True, nosave=True, debug=debug,
                            freqsplat=freqsplat, nofilt=nofilt)
            except KeyError as e:
                warnings.warn("Missing key. Bad file? {}: {}".format(s,
                                                                     str(e)))
                continue
            except Exception as e:
                warnings.warn("Error while processing {}: {}".format(s,
                                                                     str(e)))
                warnings.warn(traceback.format_exc())
                continue

            feeds = np.arange(scan['ra'].shape[1])
            chans = scan.chan_columns()

            chan_nums = np.arange(len(chans))
            F, N = np.meshgrid(feeds, chan_nums)
            F = F.flatten()
            N = N.flatten()
            for feed, nch in zip(F, N):
                channel = chans[nch]

                ras = np.degrees(scan['ra'][:, feed])
                decs = np.degrees(scan['dec'][:, feed])
                time = np.mean(scan['time'][:])
                el = np.degrees(np.mean(scan['el'][:, feed]))
                az = np.degrees(np.mean(scan['az'][:, feed]))
                source = scan.meta['SOURCE']
                pnt_ra = np.degrees(scan.meta['RA'])
                pnt_dec = np.degrees(scan.meta['Dec'])
                frequency = scan[channel].meta['frequency']
                bandwidth = scan[channel].meta['bandwidth']
                flux_density, flux_density_err = 0, 0
                flux_over_counts, flux_over_counts_err = 0, 0

                y = scan[channel]

                x, scan_type = _scantype(ras, decs)

                model, fit_info = fit_baseline_plus_bell(x, y, kind='gauss')

                try:
                    uncert = fit_info['param_cov'].diagonal() ** 0.5
                except:
                    message = fit_info['message']
                    warnings.warn(
                        "Fit failed in scan {s}: {m}".format(s=s,
                                                             m=message))
                    continue
                bell = model['Bell']
                # pars = model.parameters
                pnames = model.param_names
                counts = model.amplitude_1.value
                if plot:
                    fig = plt.figure("Fit information")
                    import matplotlib as mpl
                    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=(3, 1))
                    ax0 = plt.subplot(gs[0])
                    ax1 = plt.subplot(gs[1], sharex=ax0)

                    ax0.plot(x, y, label="Data")
                    ax0.plot(x, bell(x), label="Fit")
                    ax1.plot(x, y - bell(x))

                if scan_type.startswith("RA"):
                    fit_ra = bell.mean
                    fit_width = bell.stddev * np.cos(np.radians(pnt_dec))
                    fit_dec = None
                    ra_err = fit_ra * u.degree - pnt_ra
                    dec_err = None
                    if plot:
                        ax0.axvline(fit_ra, label="RA Fit", ls="-")
                    if plot:
                        ax0.axvline(pnt_ra.to(u.deg).value, label="RA Pnt",
                                    ls="--")
                        ax0.set_xlim([fit_ra - 2, fit_ra + 2])
                        ax1.set_xlabel('RA')

                elif scan_type.startswith("Dec"):
                    fit_ra = None
                    fit_dec = bell.mean
                    fit_width = bell.stddev
                    dec_err = fit_dec * u.degree - pnt_dec
                    ra_err = None
                    if plot:
                        ax0.axvline(fit_dec, label="Dec Fit", ls="-")
                    if plot:
                        ax0.axvline(pnt_dec.to(u.deg).value, label="Dec Pnt",
                                    ls="--")
                        ax0.set_xlim([fit_dec - 2, fit_dec + 2])
                        ax1.set_xlabel('Dec')
                        ax0.set_ylabel("Counts")
                        ax1.set_ylabel("Residual (cts)")
                index = pnames.index("amplitude_1")

                counts_err = uncert[index]
                index = pnames.index("stddev_1")
                width_err = uncert[index]

                self.add_row([scandir, sname, scan_type, source, channel, feed,
                              time, frequency, bandwidth, counts, counts_err,
                              fit_width, width_err,
                              flux_density, flux_density_err, el, az,
                              flux_over_counts, flux_over_counts_err,
                              flux_over_counts, flux_over_counts_err,
                              pnt_ra, pnt_dec, fit_ra, fit_dec, ra_err,
                              dec_err])

                if plot:
                    plt.legend()
                    plt.savefig(os.path.join(outdir,
                                             "Feed{}_chan{}.png".format(feed,
                                                                        nch)))
                    plt.close(fig)

    def write(self, fname, *args, **kwargs):
        if fname.endswith('.hdf5'):
            Table.write(self, fname, *args, path='table', **kwargs)
        else:
            Table.write(self, fname, *args, **kwargs)


class CalibratorTable(SourceTable):
    """Class containing all information and functions about calibrators."""

    def __init__(self, *args, **kwargs):
        """Initialize the object."""
        SourceTable.__init__(self, *args, **kwargs)
        self.calibration_coeffs = {}
        self.calibration_uncerts = {}
        self.calibration = {}

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
        """Update the calibration information."""
        if not self.check_not_empty():
            return

        self.get_fluxes()
        self.calibrate()
        self.compute_conversion_function()

    def get_fluxes(self):
        """Get the tabulated flux of the calibrator."""
        if not self.check_not_empty():
            return

        for it, t in enumerate(self['Time']):
            source = self['Source'][it].decode("utf-8")
            frequency = self['Frequency'][it] / 1000
            bandwidth = self['Bandwidth'][it] / 1000
            flux, eflux = \
                _get_calibrator_flux(source, frequency, bandwidth, time=t)

            self['Flux'][it] = flux
            self['Flux Err'][it] = eflux

    def calibrate(self):
        """Calculate the calibration constants."""
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

    def compute_conversion_function(self, map_unit="Jy/beam"):
        """Compute the conversion between Jy and counts.

        Try to get a meaningful fit over elevation. Revert to the rough
        function `Jy_over_counts_rough` in case `statsmodels` is not installed.
        """
        try:
            import statsmodels.api as sm
        except:
            channels = list(set(self["Chan"]))
            for channel in channels:
                fc, fce = self.Jy_over_counts_rough(channel=channel,
                                                    map_unit=map_unit)
                self.calibration_coeffs[channel] = [fc, 0, 0]
                self.calibration_uncerts[channel] = [fce, 0, 0]
                self.calibration[channel] = None
            return

        flux_quantity = _get_flux_quantity(map_unit)

        channels = list(set(self["Chan"]))
        for channel in channels:
            good_chans = self["Chan"] == channel

            f_c_ratio = self[flux_quantity + "/Counts"][good_chans]
            f_c_ratio_err = self[flux_quantity + "/Counts Err"][good_chans]
            elvs = self["Elevation"][good_chans]

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

            model = sm.RLM(y_to_fit, X)
            results = model.fit()

            self.calibration_coeffs[channel] = results.params
            self.calibration_uncerts[channel] = \
                results.cov_params().diagonal()**0.5
            self.calibration[channel] = results

    def Jy_over_counts(self, channel, elevation=None, map_unit="Jy/beam"):
        rough = False
        try:
            import statsmodels.api as sm
            from statsmodels.sandbox.regression.predstd import \
                wls_prediction_std
        except:
            warnings.warn("Statsmodels is not installed. "
                          "Reverting to rough mode.")
            rough = True

        flux_quantity = _get_flux_quantity(map_unit)

        if hasattr(channel, 'encode'):
            channel = channel.encode()

        if channel not in self.calibration.keys():
            self.compute_conversion_function(map_unit)

        if elevation is None or rough is True:
            elevation = np.array(elevation)
            fc, fce = self.Jy_over_counts_rough(channel=channel,
                                                map_unit=map_unit)
            if elevation.size > 1:
                fc = np.zeros_like(elevation) + fc
                fce = np.zeros_like(elevation) + fce
            return fc, fce

        X = np.column_stack((np.ones(np.array(elevation).size),
                             np.array(elevation)))

        fc = self.calibration[channel].predict(X)

        goodch = self["Chan"] == channel
        fce = np.mean(self[flux_quantity + " Err"][goodch]) + np.zeros_like(fc)

        if len(fc) == 1:
            fc, fce = fc[0], fce[0]

        return fc, fce

    def Jy_over_counts_rough(self, channel=None, map_unit="Jy/beam"):
        """Get the conversion from counts to Jy.

        Other parameters
        ----------------
        channel : str
            Name of the data channel

        Returns
        -------
        fc : float
            flux density /count ratio
        fce : float
            uncertainty on `fc`
        """

        if hasattr(channel, 'encode'):
            channel = channel.encode()

        self.check_up_to_date()

        good_chans = np.ones(len(self["Time"]), dtype=bool)
        if channel is not None:
            good_chans = self['Chan'] == channel

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

        p = [np.mean(y_to_fit)]
        while 1:
            p, pcov = curve_fit(_constant, x_to_fit, y_to_fit, sigma=ye_to_fit,
                                p0=p)

            bad = np.abs((y_to_fit - _constant(x_to_fit, p)) / ye_to_fit) > 5

            if not np.any(bad):
                break

            if len(x_to_fit[bad]) > len(x_to_fit) - 5:
                warnings.warn("Calibration fit is shaky")
                break

            xbad = x_to_fit[bad]
            ybad = y_to_fit[bad]
            for xb, yb in zip(xbad, ybad):
                logging.info("Outliers: {}, {}".format(xb, yb))
            good = np.logical_not(bad)
            x_to_fit = x_to_fit[good]
            y_to_fit = y_to_fit[good]
            ye_to_fit = ye_to_fit[good]

        fc = p[0]
        fce = np.sqrt(pcov[0, 0])

        return fc, fce

    def beam_width(self, channel=None):
        goodch = np.ones(len(self), dtype=bool)
        if channel is not None:
            goodch = self["Chan"] == channel
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
                         xfactor=1, yfactor=1, color=None):
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

        if showit:
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
            if type(channel) == bytes and not type(channel) == str:
                print("DEcoding")
                channel_str = channel.decode()
            else:
                channel_str = channel
            color = colors[ic]
            self.plot_two_columns('Elevation', "Flux/Counts",
                                  yerrcol="Flux/Counts Err", ax=ax00,
                                  channel=channel, color=color)

            elevations = np.arange(np.min(self['Elevation']),
                                   np.max(self['Elevation']), 0.001)
            jy_over_cts, jy_over_cts_err = self.Jy_over_counts(channel_str,
                                                               elevations)
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

            jy_over_cts, jy_over_cts_err = self.Jy_over_counts(channel_str, 45)
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
    """Flux function from Perley & Butler ApJS 204, 19 (2013)."""
    a0, a1, a2, a3 = coeffs

    if np.all(ecoeffs < 1e10):
        # assume 5% error on calibration parameters!
        ecoeffs = coeffs * 0.05
    a0e, a1e, a2e, a3e = ecoeffs
    f0 = start_frequency
    f1 = start_frequency + bandwidth

    fs = np.linspace(f0, f1, 21)
    df = np.diff(fs)[0]

    logf = np.log10(fs)
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
    fobj = io.BytesIO(coefftable.encode())
    table = Table.read(fobj, format='ascii.csv')

    idx = np.argmin(np.abs(np.longdouble(table["time"]) - time))

    a0, a0e = table['a0', 'a0e'][idx]
    a1, a1e = table['a1', 'a1e'][idx]
    a2, a2e = table['a2', 'a2e'][idx]
    a3, a3e = table['a3', 'a3e'][idx]
    coeffs = np.array([a0, a1, a2, a3], dtype=float)

    ecoeffs = np.array([a0e, a1e, a2e, a3e], dtype=float)

    return flux_function(frequency, bandwidth, coeffs, ecoeffs)


def main(args=None):
    """Main function."""
    import argparse
    import os

    description = ('Load a series of scans from a config file '
                   'and produce a map.')
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

    args = parser.parse_args(args)

    if args.sample_config:
        sample_config_file()
        sys.exit()

    if args.file is not None:
        caltable = CalibratorTable().read(args.file)
        caltable.show()
        sys.exit()
    assert args.config is not None, "Please specify the config file!"

    config = read_config(args.config)

    calibrator_dirs = config['calibrator_directories']
    if calibrator_dirs is None:
        warnings.warn("No calibrators specified in config file")
        return
    scan_list = \
        list_scans(config['datadir'],
                   config['calibrator_directories'])

    scan_list.sort()

    outfile = args.output
    if outfile is None:
        outfile = args.config.replace(".ini", "_cal.hdf5")
    caltable = CalibratorTable()
    caltable.from_scans(scan_list, freqsplat=args.splat, nofilt=args.nofilt)
    caltable.update()

    if args.show:
        caltable.show()

    caltable.write(outfile, overwrite=True)
