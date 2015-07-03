"""Functions for the calibration of scans"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from .scan import Scan, list_scans
from .read_config import read_config
from .fit import fit_baseline_plus_bell
import os
import glob
import re

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
# For Python 2 and 3 compatibility
try:
    import configparser
except ImportError:
    import ConfigParser as configparser

CALIBRATOR_CONFIG = None

def read_calibrator_config():
    """Read the configuration of calibrators in data/calibrators"""
    flux_re = re.compile(r'^Flux')
    curdir = os.path.dirname(__file__)
    calibdir = os.path.join(curdir, '..', '..', 'data', 'calibrators')
    calibrator_file_list = glob.glob(os.path.join(calibdir, '*.ini'))
    configs = {}
    for cfile in calibrator_file_list:
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


def flux_function(frequency, bandwidth, coeffs, ecoeffs):
    a0, a1, a2, a3 = coeffs
    a0e, a1e, a2e, a3e = ecoeffs
    f0 = frequency - bandwidth / 2
    f1 = frequency + bandwidth / 2

    fs = np.linspace(f0, f1, 21)
    df = np.diff(fs)[0]

    logf = np.log10(fs)
    logS = a0 + a1 * logf + a2 * logf**2 + a3 * logf**3

    S = 10 ** logS

    return (np.sum(S) * df)


def calc_flux_from_coeffs(conf, frequency, bandwidth=1, time=0):
    import io
    coefftable = conf["CoeffTable"]["coeffs"]
    fobj = io.BytesIO(coefftable.encode())
    table = Table.read(fobj, format='ascii.csv')
    # print(table)
    idx = np.argmin(np.abs(table["time"] - time))
    # print(idx, table[idx])
    a0, a0e = table['a0', 'a0e'][idx]
    a1, a1e = table['a1', 'a1e'][idx]
    a2, a2e = table['a2', 'a2e'][idx]
    a3, a3e = table['a3', 'a3e'][idx]
    coeffs = np.array([a0, a1, a2, a3])

    ecoeffs = np.array([a0e, a1e, a2e, a3e])

    return flux_function(frequency, bandwidth, coeffs, ecoeffs), 0


def get_calibrator_flux(calibrator, frequency, bandwidth=1, time=0):
    global CALIBRATOR_CONFIG

    if CALIBRATOR_CONFIG is None:
        CALIBRATOR_CONFIG = read_calibrator_config()

    calibrators = CALIBRATOR_CONFIG.keys()
    if calibrator not in calibrators:
        return None, None

    conf = CALIBRATOR_CONFIG[calibrator]
    # find closest value among frequencies
    if conf["Kind"] == "FreqList":
        idx = (np.abs(np.array(conf["Frequencies"]) - frequency)).argmin()
        return conf["Fluxes"][idx] * bandwidth, \
            conf["Flux Errors"][idx] * bandwidth
    elif conf["Kind"] == "CoeffTable":
        return calc_flux_from_coeffs(conf, frequency, bandwidth, time)


calist =            ['3C147', '3C48', '3C123', '3C295', '3C286', 'NGC7027']
fluxlist = np.array([5.1885, 3.6722, 11.0837, 4.1610, 5.7194, 5.6107])
efluxlist = fluxlist * 0.05
colors = ['k', 'b', 'r', 'g', 'c', 'm']
colors = dict(zip(calist, colors))


def get_fluxes(basedir, scandir, channel='Ch0', feed=0, plotall=False):
    """Get fluxes from all scans in path."""
    dname = os.path.basename(scandir)
    scan_list = \
        list_scans(basedir, [scandir])

    scan_list.sort()
    output_table = Table(names=["Dir", "File", "Source", "Time",
                                "Frequency", "Bandwidth",
                                "Counts", "Counts Err",
                                "Flux Density", "Flux Density Err",
                                "Kind",
                                "Elevation",
                                "Flux/Counts", "Flux/Counts Err"],
                         dtype=['U200', 'U200', 'U200', np.longdouble,
                                np.float, np.float,
                                np.float, np.float,
                                np.float, np.float,
                                "U200",
                                np.float,
                                np.float, np.float])

    if plotall:
        figures = []
    for s in scan_list:
        sname = os.path.basename(s)
        try:
            # For now, use nosave. HDF5 doesn't store meta, essential for this
            scan = Scan(s, norefilt=True, nosave=True)
        except:
            print('{} is an invalid file'.format(s))
            continue
        ras = np.degrees(scan['ra'][:, feed])
        decs = np.degrees(scan['dec'][:, feed])
        time = np.mean(scan['time'][:])
        el = np.degrees(np.mean(scan['el'][:, feed]))
        source = scan.meta['SOURCE']

        frequency = scan[channel].meta['frequency']

        bandwidth = scan[channel].meta['bandwidth']

        # Note: Formulas are in GHz here.
        flux, eflux = \
            get_calibrator_flux(source, frequency / 1000,
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
            xvariab = 'DEC'
        model, fit_info = fit_baseline_plus_bell(x, y, kind='gauss')

        try:
            uncert = fit_info['param_cov'].diagonal() ** 0.5
        except:
            print("fit failed")
            continue

        pars = model.parameters
        pnames = model.param_names
        counts = model.amplitude_1.value

        index = pnames.index("amplitude_1")

        counts_err = uncert[index]

        if plotall:
            baseline = model['Baseline']
            bell = model['Bell']
            figure_name = source + '_' + xvariab
            if figure_name not in figures:
                figures.append(figure_name)
            plt.figure(figure_name)
            plt.plot(x, y - baseline(x), label='{:.2f}'.format(el))
            plt.plot(x, bell(x))
            plt.title('{} (baseline-subtracted)'.format(source))
            plt.xlabel(xvariab)

        flux_over_counts = flux / counts
        flux_over_counts_err = \
            (counts_err / counts + eflux / flux) * flux_over_counts

        output_table.add_row([scandir, sname, source, time,
                              frequency, bandwidth, counts, counts_err,
                              flux_density, flux_density_err, kind, el,
                              flux_over_counts, flux_over_counts_err])

    if plotall:
        for f in figures:
            fig = plt.figure(f)
            plt.legend()
            fig.savefig(f + ".png")
        # plt.close(fig)
    return output_table


def get_full_table(config_file, channel='Ch0', feed=0, plotall=False):
    """Get all fluxes in the directories specified by the config file"""
    config = read_config(config_file)

    dir_list = config['list_of_directories']
    tables = {}
    for d in dir_list:
        output_table = get_fluxes(config['datadir'], d, channel=channel,
                                  feed=feed, plotall=plotall)
        tables[d] = output_table

    return Table(vstack(list(tables.values())))


def show_calibration(full_table, feed=0, plotall=False):
    """Show the results of calibration"""

    dir_list = list(set(full_table["Dir"]))

    calibrator_table = full_table[full_table["Kind"] == "Calibrator"]

    source_table = full_table[full_table["Kind"] == "Source"]

    for d in dir_list:
        subtable = calibrator_table[calibrator_table["Dir"] == d]
        if len(subtable) == 0:
            continue

        fig = plt.figure("Time evolution")

        fc = np.mean(subtable["Flux/Counts"]) / subtable["Bandwidth"][0]
        fce = np.sqrt(np.sum(subtable["Flux/Counts Err"] ** 2)) / len(subtable) / subtable["Bandwidth"][0]

        plt.errorbar(np.mean(subtable["Time"]), fc, yerr=fce,
                     ecolor=colors[subtable["Source"][0]],
                     elinewidth=3)

        plt.figure("Vs Elevation")
        plt.errorbar(np.mean(subtable["Elevation"]), fc, yerr=fce,
                     ecolor=colors[subtable["Source"][0]],
                     elinewidth=3)

        counts = np.mean(subtable["Counts"])
        counts_err = \
            np.sqrt(np.sum(subtable["Counts Err"] ** 2)) / len(subtable)
        plt.figure("Vs Flux")
        plt.errorbar(counts, fc, yerr=fce, xerr=counts_err,
                     ecolor=colors[subtable["Source"][0]],
                     elinewidth=3)

    plt.figure("Time evolution")
    plt.xlabel("Time")
    plt.ylabel("Jansky / Counts")
    plt.figure("Vs Elevation")
    plt.ylabel("Jansky / Counts")
    plt.xlabel("Elevation")
    plt.figure("Vs Flux")
    plt.ylabel("Jansky / Counts")
    plt.xlabel("Flux Density")

    fc = np.mean(calibrator_table["Flux/Counts"])
    fce = np.sqrt(np.sum(calibrator_table["Flux/Counts Err"] ** 2))\
        / len(calibrator_table)

    source_table["Flux Density"] = \
        source_table["Counts"] * fc / source_table["Bandwidth"]
    source_table["Flux Density Err"] = \
        (source_table["Counts Err"] / source_table["Counts"] + fce / fc) * \
        source_table["Flux Density"] / source_table["Bandwidth"]

    sources = list(set(source_table["Source"]))
    for s in sources:
        filtered = source_table[source_table["Source"] == s]
        print(filtered[("Source", "Flux Density", "Flux Density Err",
                        "Counts", "Counts Err")])
        if len(filtered) > 20:
            plt.figure(s)
            plt.hist(filtered["Flux Density"], bins=len(filtered))
            plt.xlabel("Flux values")

    plt.show()


def test_calibration_tp():
    curdir = os.path.abspath(os.path.dirname(__file__))
    config_file = \
        os.path.abspath(os.path.join(curdir, '..', '..',
                                     'TEST_DATASET',
                                     'test_calib.ini'))
    full_table = get_full_table(config_file, plotall=True)
    show_calibration(full_table)


def test_calibration_roach():
    curdir = os.path.abspath(os.path.dirname(__file__))
    config_file = \
        os.path.abspath(os.path.join(curdir, '..', '..',
                                     'TEST_DATASET',
                                     'test_calib_roach.ini'))
    full_table = get_full_table(config_file)
    show_calibration(full_table)
