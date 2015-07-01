"""Functions for the calibration of scans"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from .scan import Scan, list_scans
from .read_config import read_config
from .fit import fit_baseline_plus_bell
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column

calist=['3C147','3C48','3C123','3C295','3C286','NGC7027']
fluxlist=[5.1885,3.6722,11.0837,4.1610,5.7194,5.6107]
colors=['k', 'b', 'r', 'g', 'c', 'm']
colors = dict(zip(calist,colors))

def calibration(config_file, channel='Ch0', feed=0):
    config = read_config(config_file)

    dir_list = config['list_of_directories']
    tables = {}
    for d in dir_list:
        scan_list = \
            list_scans(config['datadir'], [d])

        scan_list.sort()

        output_table = Table(names=["Source", "Time", "Amplitude", "Kind",
                                    "Elevation", "Counts/flux",
                                    "Counts/flux Err"],
                             dtype=['S10', np.longdouble, np.float, "S10",
                                    np.float, np.float, np.float])

        for s in scan_list:
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
                kind = "Calibrator"
            except:
                flux = 1.
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
            pars = model.parameters
            pnames = model.param_names
            counts = model.amplitude_1.value

            index = pnames.index("amplitude_1")

            counts_err = uncert[index]
            # plt.plot(x, y)
            # plt.plot(x, model(x))
            # plt.title(source)
            # plt.show()

            count_over_flux = flux / counts
            count_over_flux_err = counts_err / counts * count_over_flux

            output_table.add_row([source, time, counts, kind,
                                  el, count_over_flux, count_over_flux_err])
        tables[d] = output_table

    for d in dir_list:
        print(d)
        print(tables[d])
        if tables[d]["Kind"][0] == b"Source":
            continue

        plt.figure("Time evolution")

        cf = np.mean(tables[d]["Counts/flux"])
        cfe = np.sqrt(np.sum(tables[d]["Counts/flux Err"] ** 2)) / len(tables[d])
        plt.errorbar(np.mean(tables[d]["Time"]), cf, yerr=cfe,
                     ecolor=colors[tables[d]["Source"][0].decode()])

        plt.figure("Vs Elevation")
        plt.errorbar(np.mean(tables[d]["Elevation"]), cf, yerr=cfe, ecolor=colors[tables[d]["Source"][0].decode()])

    plt.show()

def test_calibration():
    curdir = os.path.abspath(os.path.dirname(__file__))
    config_file = \
        os.path.abspath(os.path.join(curdir, '..', '..',
                                     'TEST_DATASET',
                                     'test_calib.ini'))
    calibration(config_file)
