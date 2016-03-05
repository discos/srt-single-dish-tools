"""
Produce calibrated light curves.

``SDTlcurve`` is a script that, given a list of cross scans from different
sources, is able to recognize calibrators and use them to convert the observed
counts into a density flux value in Jy.
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from .scan import Scan, list_scans
from .read_config import read_config, sample_config_file, get_config_file
from .fit import fit_baseline_plus_bell
import os
import sys
import glob
import re
import warnings
import traceback

try:
    import pickle
except:
    import cPickle as pickle

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, Column
# For Python 2 and 3 compatibility
try:
    import configparser
except ImportError:
    import ConfigParser as configparser

CALIBRATOR_CONFIG = None


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
    print("Reading calibrator files ")

    configs = {}
    for cfile in calibrator_file_list:
        print(cfile)
        cparser = configparser.ConfigParser()
        cparser.read(cfile)
        cparser["Info"]["Name"]
        if 'CoeffTable' not in list(cparser.sections()):
            configs[cparser["Info"]["Name"]] = {"Kind": "FreqList",
                                                "Frequencies": [],
                                                "Bandwidths": [],
                                                "Fluxes": [],
                                                "Flux Errors": []}

            for section in cparser.sections():
                if not flux_re.match(section):
                    continue
                configs[cparser["Info"]["Name"]]["Frequencies"].append(
                    float(cparser[section]["freq"]))
                configs[cparser["Info"]["Name"]]["Bandwidths"].append(
                    float(cparser[section]["bwidth"]))
                configs[cparser["Info"]["Name"]]["Fluxes"].append(
                    float(cparser[section]["flux"]))
                configs[cparser["Info"]["Name"]]["Flux Errors"].append(
                    float(cparser[section]["eflux"]))
        else:
            configs[cparser["Info"]["Name"]] = \
                {"CoeffTable": cparser["CoeffTable"], "Kind": "CoeffTable"}

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
                 "Width",
                 "Flux Density", "Flux Density Err",
                 "Elevation",
                 "Flux/Counts", "Flux/Counts Err",
                 "RA", "Dec",
                 "Fit RA", "Fit Dec"]

        dtype = ['U200', 'U200', 'U200', 'U200',
                 'U200', np.int, np.longdouble,
                 np.float, np.float,
                 np.float, np.float,
                 np.float,
                 np.float, np.float,
                 np.float,
                 np.float, np.float,
                 np.float, np.float,
                 np.float, np.float]

        for n, d in zip(names, dtype):
            if n not in self.keys():
                self.add_column(Column(name=n, dtype=d))

    def from_scans(self, scan_list=None, verbose=False, freqsplat=None,
                   config_file=None):
        """Load source table from a list of scans."""

        if scan_list is None:
            if config_file is None:
                config_file = get_config_file()
            config = read_config(config_file)
            scan_list = \
                list_scans(config['datadir'], config['list_of_directories'])
            scan_list.sort()
        nscan = len(scan_list)

        for i_s, s in enumerate(scan_list):
            print('{}/{}: Loading {}'.format(i_s + 1, nscan, s))
            scandir, sname = os.path.split(s)
            try:
                # For now, use nosave. HDF5 doesn't store meta, essential for
                # this
                # TODO: experiment with serialize_meta!
                scan = Scan(s, norefilt=True, nosave=True, verbose=verbose,
                            freqsplat=freqsplat)
            except KeyError as e:
                warnings.warn("Error while processing {}: {}".format(s,
                                                                     str(e)))
            except Exception as e:
                traceback.print_exc()
                warnings.warn("Error while processing {}: {}".format(s,
                                                                     str(e)))

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
                    warnings.warn("Fit failed in scan {s}".format(s=s))
                    continue

                bell = model['Bell']
                # pars = model.parameters
                pnames = model.param_names
                counts = model.amplitude_1.value

                if scan_type.startswith("RA"):
                    fit_ra = bell.mean
                    fit_width = bell.stddev * np.cos(np.radians(pnt_dec))
                    fit_dec = None
                elif scan_type.startswith("Dec"):
                    fit_ra = None
                    fit_dec = bell.mean
                    fit_width = bell.stddev

                index = pnames.index("amplitude_1")

                counts_err = uncert[index]

                self.add_row([scandir, sname, scan_type, source, channel, feed,
                              time, frequency, bandwidth, counts, counts_err,
                              fit_width,
                              flux_density, flux_density_err, el,
                              flux_over_counts, flux_over_counts_err,
                              pnt_ra, pnt_dec, fit_ra, fit_dec])


class CalibratorTable(SourceTable):
    """Class containing all information and functions about calibrators."""

    def __init__(self, *args, **kwargs):
        """Initialize the object."""
        SourceTable.__init__(self, *args, **kwargs)

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

    def get_fluxes(self):
        """Get the tabulated flux of the calibrator."""
        if not self.check_not_empty():
            return

        for it, t in enumerate(self['Time']):
            source = self['Source'][it]
            frequency = self['Frequency'][it] / 1000
            bandwidth = self['Bandwidth'][it] / 1000
            flux, eflux = \
                _get_calibrator_flux(source, frequency, bandwidth, time=t)

            self['Flux Density'][it] = flux
            self['Flux Density Err'][it] = eflux
            print(flux, eflux)

    def calibrate(self):
        """Calculate the calibration constants."""
        if not self.check_not_empty():
            return

        flux = self['Flux Density']
        eflux = self['Flux Density Err']
        counts = self['Counts']
        ecounts = self['Counts Err']

        flux_over_counts = flux / counts
        flux_over_counts_err = \
            (ecounts / counts + eflux / flux) * flux_over_counts

        self['Flux/Counts'][:] = flux_over_counts
        self['Flux/Counts Err'][:] = flux_over_counts_err

    def Jy_over_counts(self, channel=None):
        """Get the conversion from counts to Jy."""
        self.check_up_to_date()

        good_chans = np.ones(len(self["Time"]), dtype=bool)
        if channel is not None:
            good_chans = self["Chan"] == channel

        f_c_ratio = self["Flux/Counts"][good_chans]
        f_c_ratio_err = self["Flux/Counts Err"][good_chans]

        good_fc = (f_c_ratio == f_c_ratio) & (f_c_ratio > 0)
        good_fce = (f_c_ratio_err == f_c_ratio_err) & (f_c_ratio_err >= 0)

        good = good_fc & good_fce

        fc = np.mean(f_c_ratio[good])

        fce = np.sqrt(np.sum(f_c_ratio_err[good] ** 2))\
            / len(f_c_ratio_err[good])

        return fc, fce

    def counts_over_Jy(self, channel=None):
        """Get the conversion from Jy to counts."""
        self.check_up_to_date()

        fc, fce = self.Jy_over_counts(channel=channel)
        cf = 1 / fc
        return cf, fce / fc * cf


def decide_symbol(values):
    """Decide symbols for plotting.

    Assigns different symbols to RA scans, Dec scans, backward and forward.
    """
    raplus = values == "RA>"
    ramin = values == "RA<"
    decplus = values == "Dec>"
    decmin = values == "Dec<"
    symbols = np.array(['a' for i in values])
    symbols[raplus] = u"+"
    symbols[ramin] = u"s"
    symbols[decplus] = u"^"
    symbols[decmin] = u"v"
    return symbols


def flux_function(start_frequency, bandwidth, coeffs, ecoeffs):
    """Flux function from Perley & Butler ApJS 204, 19 (2013)."""
    a0, a1, a2, a3 = coeffs
    print(coeffs,ecoeffs)
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

    # Error is not random, should add linearly
    return np.sum(S) * df, np.sum(eS) * df


def _calc_flux_from_coeffs(conf, frequency, bandwidth=1, time=0):
    """Return the flux of a calibrator at a given frequency.

    Uses Perley & Butler ApJS 204, 19 (2013).
    """
    import io
    coefftable = conf["CoeffTable"]["coeffs"]
    fobj = io.BytesIO(coefftable.encode())
    table = Table.read(fobj, format='ascii.csv')

    idx = np.argmin(np.abs(table["time"] - time))

    a0, a0e = table['a0', 'a0e'][idx]
    a1, a1e = table['a1', 'a1e'][idx]
    a2, a2e = table['a2', 'a2e'][idx]
    a3, a3e = table['a3', 'a3e'][idx]
    coeffs = np.array([a0, a1, a2, a3])

    ecoeffs = np.array([a0e, a1e, a2e, a3e])

    return flux_function(frequency, bandwidth, coeffs, ecoeffs)


calist = ['3C147', '3C48', '3C123', '3C295', '3C286', 'NGC7027']
colors = ['k', 'b', 'r', 'g', 'c', 'm']
colors = dict(zip(calist, colors))


def get_fluxes(basedir, scandir, channel='Ch0', feed=0, plotall=False,
               verbose=True, freqsplat=None):
    """Get fluxes from all scans in path."""
    # dname = os.path.basename(scandir)
    scan_list = \
        list_scans(basedir, [scandir])

    scan_list.sort()
    output_table = Table(names=["Dir", "File", "Scan Type", "Source", "Time",
                                "Frequency", "Bandwidth",
                                "Counts", "Counts Err",
                                "Width",
                                "Flux Density", "Flux Density Err",
                                "Kind",
                                "Elevation",
                                "Flux/Counts", "Flux/Counts Err",
                                "RA", "Dec",
                                "Fit RA", "Fit Dec"],
                         dtype=['U200', 'U200', 'U200', 'U200', np.longdouble,
                                np.float, np.float,
                                np.float, np.float,
                                np.float,
                                np.float, np.float,
                                "U200",
                                np.float,
                                np.float, np.float,
                                np.float, np.float,
                                np.float, np.float])

    if plotall:
        figures = []
        plotted_kinds = []

    nscan = len(scan_list)

    for i_s, s in enumerate(scan_list):
        print('{}/{}: Loading {}'.format(i_s + 1, nscan, s))
        sname = os.path.basename(s)
        try:
            # For now, use nosave. HDF5 doesn't store meta, essential for this
            scan = Scan(s, norefilt=True, nosave=True, verbose=verbose,
                        freqsplat=freqsplat)
        except:
            warnings.warn('{} is an invalid file'.format(s))
            continue
        ras = np.degrees(scan['ra'][:, feed])
        decs = np.degrees(scan['dec'][:, feed])
        time = np.mean(scan['time'][:])
        el = np.degrees(np.mean(scan['el'][:, feed]))
        source = scan.meta['SOURCE']
        backend = scan.meta['backend']
        pnt_ra = np.degrees(scan.meta['RA'])
        pnt_dec = np.degrees(scan.meta['Dec'])

        frequency = scan[channel].meta['frequency']

        bandwidth = scan[channel].meta['bandwidth']

        # Note: Formulas are in GHz here.
        flux, eflux = \
            _get_calibrator_flux(source, frequency / 1000,
                                 bandwidth / 1000, time=time)

        if flux is None:
            flux_density = 1 / bandwidth
            flux_density_err = 0
            flux = 1.
            eflux = 0
            kind = "Source"
        else:
            flux *= 1000
            eflux *= 1000
            # Config gives flux density (... Hz^-1). Normalize by bandwidth
            flux_density = flux / bandwidth
            flux_density_err = eflux / bandwidth
            kind = "Calibrator"

        y = scan[channel]

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
        scan_type = xvariab + scan_direction

        model, fit_info = fit_baseline_plus_bell(x, y, kind='gauss')

        try:
            uncert = fit_info['param_cov'].diagonal() ** 0.5
        except:
            warnings.warn("Fit failed in scan {s}".format(s=s))
            continue

        baseline = model['Baseline']
        bell = model['Bell']
        # pars = model.parameters
        pnames = model.param_names
        counts = model.amplitude_1.value
        if xvariab == "RA":
            fit_ra = bell.mean
            fit_width = bell.stddev * np.cos(np.radians(pnt_dec))
            fit_dec = None
            to_plot = pnt_ra
        elif xvariab == "Dec":
            fit_ra = None
            fit_dec = bell.mean
            to_plot = pnt_dec
            fit_width = bell.stddev

        index = pnames.index("amplitude_1")

        counts_err = uncert[index]

        if plotall:
            figure_name = '{}_{}_{}'.format(source, xvariab, backend)
            first = False
            if figure_name not in figures:
                first = True
                figures.append(figure_name)
                plotted_kinds.append(kind)

            plt.figure(figure_name)

            if first:
                plt.axvline(to_plot)
            data = plt.plot(x, np.array(y - baseline(x)),
                            label='{:.2f}'.format(el))
            plt.plot(x, bell(x), color=plt.getp(data[0], 'color'))
            plt.title('{} (baseline-subtracted)'.format(source))
            plt.xlabel(xvariab)

        flux_over_counts = flux / counts
        flux_over_counts_err = \
            (counts_err / counts + eflux / flux) * flux_over_counts

        output_table.add_row([scandir, sname, scan_type, source, time,
                              frequency, bandwidth, counts, counts_err,
                              fit_width,
                              flux_density, flux_density_err, kind, el,
                              flux_over_counts, flux_over_counts_err,
                              pnt_ra, pnt_dec, fit_ra, fit_dec])

    if plotall:
        for i_f, f in enumerate(figures):
            fig = plt.figure(f)
            if plotted_kinds[i_f] == 'Calibrator':
                plt.legend()
            fig.savefig(f + ".png")

    return output_table


def get_full_table(config_file, channel='Ch0', feed=0, plotall=False,
                   picklefile=None, verbose=True, freqsplat=None):
    """Get all fluxes in the directories specified by the config file."""
    config = read_config(config_file)

    dir_list = config['list_of_directories']
    tables = {}

    ndir = len(dir_list)
    for i_d, d in enumerate(dir_list):
        print('\n-----------------\n')
        print('{}/{}: Loading data in {}'.format(i_d + 1, ndir, d))
        print('\n-----------------\n')

        output_table = get_fluxes(config['datadir'], d, channel=channel,
                                  feed=feed, plotall=plotall, verbose=verbose,
                                  freqsplat=freqsplat)
        tables[d] = output_table

    full_table = Table(vstack(list(tables.values())))
    if picklefile is not None:
        with open(picklefile, 'wb') as f:
            pickle.dump(full_table, f)
    return full_table


def show_calibration(full_table, feed=0, plotall=False):
    """Show the results of calibration."""
    import matplotlib as mpl

    dir_list = list(set(full_table["Dir"]))

    calibrator_table = full_table[full_table["Kind"] == "Calibrator"]

    source_table = full_table[full_table["Kind"] == "Source"]

    for d in dir_list:
        subtable = calibrator_table[calibrator_table["Dir"] == d]
        if len(subtable) == 0:
            continue

        symbols = decide_symbol(subtable["Scan Type"])
        for k in colors.keys():
            if k in subtable["Source"][0]:
                source_color = colors[k]
                break
        else:
            source_color = 'grey'

        # ----------------------- Pointing vs. ELEVATION -------------------
        plt.figure("Pointing Error vs Elevation")
        good_ra = subtable["Fit RA"] == subtable["Fit RA"]
        good_dec = subtable["Fit Dec"] == subtable["Fit Dec"]

        ra_pnt = subtable["RA"]
        dec_pnt = subtable["Dec"]
        ra_fit = subtable["Fit RA"]
        dec_fit = subtable["Fit Dec"]

        el = subtable["Elevation"]
        ra_err = (ra_fit - ra_pnt) / np.cos(np.radians(dec_pnt))
        dec_err = dec_fit - dec_pnt
        pointing_err = np.sqrt(np.mean(ra_err[good_ra])**2 +
                               np.mean(dec_err[good_dec])**2)

        plt.scatter(np.mean(el), pointing_err * 60, color='k', marker='o')
        for _e, _r, _d,  _s in zip(el, ra_err, dec_err, symbols):
            plt.scatter(_e, _r * 60, color='r', marker=_s)
            plt.scatter(_e, _d * 60, color='b', marker=_s)

        fc = np.mean(subtable["Flux/Counts"]) / subtable["Bandwidth"][0]
        fce = np.sqrt(
            np.sum(subtable["Flux/Counts Err"] ** 2)) / \
            len(subtable) / subtable["Bandwidth"][0]
        fce = np.max([fce, np.std(fc)])

        # ----------------------- Calibration vs. ELEVATION -------------------
        plt.figure("Vs Elevation")
        plt.errorbar(np.mean(subtable["Elevation"]), fc, yerr=fce,
                     ecolor=source_color,
                     elinewidth=3)

        # ----------------------- Width vs. ELEVATION -------------------
        ras = np.char.rstrip(subtable["Scan Type"], "><") == "RA"
        decs = np.char.rstrip(subtable["Scan Type"], "><") == "Dec"

        plt.figure("Width Vs Elevation")
        plt.scatter(subtable["Elevation"][ras], subtable["Width"][ras],
                    color=source_color, marker='o')
        plt.scatter(subtable["Elevation"][decs], subtable["Width"][decs],
                    color=source_color, marker='^')

    plt.figure("Vs Elevation")
    plt.ylabel("Jansky / Counts")
    plt.xlabel("Elevation")

    plt.figure("Width Vs Elevation")
    plt.xlabel("Elevation")
    plt.ylabel("Gaussian Width (deg)")

    plt.figure("Pointing Error vs Elevation")
    plt.title("Pointing Error vs Elevation (black: total; red: RA; blue: Dec)")
    plt.xlabel('Elevation')
    plt.ylabel('Pointing error (arcmin)')

    rap_symb = mpl.lines.Line2D([0], [0], linestyle="none", c='r', marker='+')
    dep_symb = mpl.lines.Line2D([0], [0], linestyle="none", c='b', marker='^')
    ram_symb = mpl.lines.Line2D([0], [0], linestyle="none", c='r', marker='s')
    dem_symb = mpl.lines.Line2D([0], [0], linestyle="none", c='b', marker='v')
    tot_symb = mpl.lines.Line2D([0], [0], linestyle="none", c='k', marker='o')

    plt.legend([rap_symb, dep_symb, ram_symb, dem_symb, tot_symb],
               ['RA>', 'Dec>', 'RA<', 'Dec<', 'Tot'], numpoints=1)

    f_c_ratio = calibrator_table["Flux/Counts"]

    good = (f_c_ratio == f_c_ratio) & (f_c_ratio > 0)
    fc = np.mean(f_c_ratio[good])
    f_c_ratio_err = calibrator_table["Flux/Counts Err"]
    good = (f_c_ratio_err == f_c_ratio_err) & (f_c_ratio_err > 0)
    fce = np.sqrt(np.sum(f_c_ratio_err[good] ** 2))\
        / len(calibrator_table)

    source_table["Flux Density"] = \
        source_table["Counts"] * fc / source_table["Bandwidth"]
    source_table["Flux Density Err"] = \
        (source_table["Counts Err"] / source_table["Counts"]) * \
        source_table["Flux Density"]
    source_table["Flux Density Systematic"] = \
        fce / fc * source_table["Flux Density"]

    sources = list(set(source_table["Source"]))
    source_colors = ['k', 'b', 'r', 'g', 'c', 'm']
    for i_s, s in enumerate(sources):
        c = source_colors[i_s % len(sources)]
        filtered = source_table[source_table["Source"] == s]

        if len(filtered) > 20:
            plt.figure(s)
            from astropy.visualization import hist
            hist(filtered["Flux Density"], bins='knuth',
                 histtype='stepfilled')

            plt.xlabel("Flux values")
        plt.figure(s + '_callc')
        plt.errorbar(filtered['Time'], filtered['Flux Density'],
                     yerr=filtered['Flux Density Err'], label=s, fmt=None,
                     ecolor=c, color=c)
        plt.fill_between(
            filtered['Time'],
            filtered['Flux Density'] - filtered['Flux Density Systematic'],
            filtered['Flux Density'] + filtered['Flux Density Systematic'],
            color=c, alpha=0.1,
            label=s + '-systematic')

        plt.legend()
        plt.xlabel('Time (MJD)')
        plt.ylabel('Flux (Jy)')
        plt.savefig(s + '_callc.png')
        np.savetxt(s + '_callc_data.txt',
                   np.array([filtered['Time'],
                             filtered['Flux Density'],
                             filtered['Flux Density Err']]).T)

        plt.figure(s + '_lc')
        plt.errorbar(filtered['Time'], filtered['Counts'],
                     yerr=filtered['Counts Err'], label=s, fmt=None)

        plt.xlabel('Time (MJD)')
        plt.ylabel('Counts')
        plt.savefig(s + '_lc.png')
        np.savetxt(s + '_lc_data.txt',
                   np.array([filtered['Time'],
                             filtered['Counts'],
                             filtered['Counts Err']]).T)

        plt.legend()
    plt.show()


def test_calibration_tp():
    """Test that the calibration executes completely."""
    import pickle
    curdir = os.path.abspath(os.path.dirname(__file__))
    config_file = \
        os.path.abspath(os.path.join(curdir, '..', '..',
                                     'TEST_DATASET',
                                     'test_calib.ini'))
    full_table = get_full_table(config_file, plotall=True,
                                picklefile='data_tp.pickle')

    with open('data_tp.pickle', 'rb') as f:
        full_table = pickle.load(f)
    show_calibration(full_table)


def test_calibration_roach():
    """Test that the calibration executes completely, ROACH version."""
    curdir = os.path.abspath(os.path.dirname(__file__))
    config_file = \
        os.path.abspath(os.path.join(curdir, '..', '..',
                                     'TEST_DATASET',
                                     'test_calib_roach.ini'))
    full_table = get_full_table(config_file, plotall=True,
                                picklefile='data_r2.pickle')

    with open('data_r2.pickle', 'rb') as f:
        full_table = pickle.load(f)
    show_calibration(full_table)


def main_lc_calibrator(args=None):
    """Main function."""
    import argparse
    import os

    description = ('Load a series of scans from a config file '
                   'and produce a map.')
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("--sample-config", action='store_true', default=False,
                        help='Produce sample config file')

    parser.add_argument("-c", "--config", type=str, default=None,
                        help='Config file')

    parser.add_argument("--pickle-file", type=str, default='db.pickle',
                        help='Name for the intermediate pickle file')

    parser.add_argument("--splat", type=str, default=None,
                        help=("Spectral scans will be scrunched into a single "
                              "channel containing data in the given frequency "
                              "range, starting from the frequency of the first"
                              " bin. E.g. '0:1000' indicates 'from the first "
                              "bin of the spectrum up to 1000 MHz above'. ':' "
                              "or 'all' for all the channels."))

    parser.add_argument("--refilt", default=False,
                        action='store_true',
                        help='Re-run the scan filtering')

    args = parser.parse_args(args)

    if args.sample_config:
        sample_config_file()
        sys.exit()

    assert args.config is not None, "Please specify the config file!"

    if not os.path.exists(args.pickle_file):
        full_table = get_full_table(args.config, plotall=True,
                                    picklefile=args.pickle_file,
                                    freqsplat=args.splat)

    with open(args.pickle_file, 'rb') as f:
        full_table = pickle.load(f)
        full_table.sort('Time')
    show_calibration(full_table)
