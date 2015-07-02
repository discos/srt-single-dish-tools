"""Functions for the calibration of scans"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from .scan import Scan, list_scans
from .read_config import read_config
from .fit import fit_baseline_plus_bell
import os

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack

calist = ['3C147', '3C48', '3C123', '3C295', '3C286', 'NGC7027']
fluxlist = np.array([5.1885, 3.6722, 11.0837, 4.1610, 5.7194, 5.6107])
efluxlist = fluxlist * 0.05
colors = ['k', 'b', 'r', 'g', 'c', 'm']
colors = dict(zip(calist, colors))


def get_fluxes(basedir, scandir, channel='Ch0', feed=0, plot=False):
    """Get fluxes from all scans in path."""
    dname = os.path.basename(scandir)
    scan_list = \
        list_scans(basedir, [scandir])

    scan_list.sort()
    output_table = Table(names=["Dir", "File", "Source", "Time", "Counts",
                                "Counts Err", "Flux", "Flux Err", "Kind",
                                "Elevation", "Flux/Counts",
                                "Flux/Counts Err"],
                         dtype=['U200', 'U200', 'U200', np.longdouble, np.float,
                                np.float, np.float, np.float,
                                "U200", np.float, np.float, np.float])
    if plot:
        fig = plt.figure(dname)

    for s in scan_list:
        sname = os.path.basename(s)
        try:
            scan = Scan(s, norefilt=True)
        except:
            print('{} is an invalid file'.format(s))
            continue
        ras = scan['ra'][:, feed]
        decs = scan['dec'][:, feed]
        time = np.mean(scan['time'][:])
        el = np.degrees(np.mean(scan['el'][:, feed]))
        source = scan.meta['SOURCE']
        try:
            index = calist.index(source)
            flux = fluxlist[index]
            eflux = efluxlist[index]
            kind = "Calibrator"
        except:
            flux = 1.
            eflux = 0
            kind = "Source"

        y = scan[channel]

        ravar = np.max(ras) - np.min(ras)
        decvar = np.max(decs) - np.min(decs)
        if ravar > decvar:
            x = ras
        else:
            x = decs
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

        if plot:
            plt.plot(x - np.min(x), y)
            plt.plot(x - np.min(x), model(x))
            plt.title(source)

        flux_over_counts = flux / counts
        flux_over_counts_err = \
            (counts_err / counts + eflux / flux) * flux_over_counts

        output_table.add_row([scandir, sname, source, time, counts, counts_err,
                              flux, eflux, kind, el, flux_over_counts,
                              flux_over_counts_err])

    if plot:
        fig.savefig(dname + ".png")
        plt.close(fig)
    return output_table


def get_full_table(config_file, channel='Ch0', feed=0, plot=False):
    """Get all fluxes in the directories specified by the config file"""
    config = read_config(config_file)

    dir_list = config['list_of_directories']
    tables = {}
    for d in dir_list:
        output_table = get_fluxes(config['datadir'], d, channel=channel,
                                  feed=feed, plot=plot)
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

        fc = np.mean(subtable["Flux/Counts"])
        fce = np.sqrt(np.sum(subtable["Flux/Counts Err"] ** 2)) / len(subtable)
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
    plt.xlabel("Flux")

    fc = np.mean(calibrator_table["Flux/Counts"])
    fce = np.sqrt(np.sum(calibrator_table["Flux/Counts Err"] ** 2))\
        / len(calibrator_table)

    source_table["Flux"] = source_table["Counts"] * fc
    source_table["Flux Err"] = \
        (source_table["Counts Err"] / source_table["Counts"] + fce / fc) * \
        source_table["Flux"]

    sources = list(set(source_table["Source"]))
    for s in sources:
        filtered = source_table[source_table["Source"] == s]
        print(filtered[("Source", "Flux", "Flux Err", "Counts", "Counts Err")])
        if len(filtered) > 20:
            plt.figure(s)
            plt.hist(filtered["Flux"], bins=len(filtered))
            plt.xlabel("Flux values")

    plt.show()


def test_calibration_tp():
    curdir = os.path.abspath(os.path.dirname(__file__))
    config_file = \
        os.path.abspath(os.path.join(curdir, '..', '..',
                                     'TEST_DATASET',
                                     'test_calib.ini'))
    full_table = get_full_table(config_file)
    show_calibration(full_table)


def test_calibration_roach():
    curdir = os.path.abspath(os.path.dirname(__file__))
    config_file = \
        os.path.abspath(os.path.join(curdir, '..', '..',
                                     'TEST_DATASET',
                                     'test_calib_roach.ini'))
    full_table = get_full_table(config_file)
    show_calibration(full_table)
