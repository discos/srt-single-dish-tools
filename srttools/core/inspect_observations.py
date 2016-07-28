"""From a given list of directories, read the relevant information and link observations to calibrators."""

from __future__ import (absolute_import, division,
                        print_function)
import os
import glob
import warnings
import numpy as np
from astropy.table import Table, Column
from .io import read_data
from .calibration import read_calibrator_config
import sys

def standard_string(s):
    if sys.version_info >= (3, 0, 0):
        # for Python 3
        if isinstance(s, bytes):
            s = s.decode('ascii')  # or  s = str(s)[2:-1]
    else:
        # for Python 2
        if isinstance(s, unicode):
            s = str(s)
    return s

def inspect_directories(directories):
    info = Table()
    names = ["Dir", "Sample File", "Source", "Receiver", "Backend",
             "Time", "Frequency", "Bandwidth"]

    dtype = ['S200', 'S200', 'S200', 'S200', 'S200',
             np.double, np.float, np.float]

    for n, d in zip(names, dtype):
        if n not in info.keys():
            info.add_column(Column(name=n, dtype=d))

    for d in directories:
        fits_files = glob.glob(os.path.join(d, '*.fits'))

        for f in fits_files:
            try:
                print("Reading {}".format(f), end="\r")
            
                data = read_data(f)
                backend = data.meta['backend']
                receiver = data.meta['receiver']
                frequency = data["Ch0"].meta['frequency']
                bandwidth = data["Ch0"].meta['bandwidth']
                source = data.meta['SOURCE']
                time = data[0]['time']

                info.add_row([d, f, source, receiver, backend,
                              time, frequency, bandwidth])
                break
            except:
                continue

    return(info)


def split_observation_table(info, max_calibrator_delay=0.4, max_source_delay=0.2):
    grouped_table = info.group_by(["Receiver", "Backend"])
    indices = grouped_table.groups.indices

    groups = {}
    for i, ind in enumerate(zip(indices[:-1], indices[1:])):
        start_row = grouped_table[ind[0]]
        print("Group {}, Backend = {}, Receiver = {}".format(i, standard_string(start_row["Backend"]),
                                                             standard_string(start_row["Receiver"])))
        s = split_by_source(grouped_table[ind[0]:ind[1]],
                            max_calibrator_delay=max_calibrator_delay,
                            max_source_delay=max_source_delay)

        receiver = start_row["Receiver"]
        backend = start_row["Backend"]

        label = standard_string(receiver) + ',' + standard_string(backend)

        groups[label] = s

    return groups


def split_by_source(info, max_calibrator_delay=0.4, max_source_delay=0.2):
    cal_config = read_calibrator_config()
    calibrators = cal_config.keys()

    sources = list(set(info["Source"]))
    # Find observation blocks of a given source
    retval = {}
    for s in sources:
        if standard_string(s) in calibrators:
            continue
        condition = info["Source"] == s
        s = standard_string(s)
        retval[s] = {}
        filtered_table = info[condition]

        start_idxs = []
        end_idxs = []
        for i, f in enumerate(filtered_table):
            if i == 0:
                start_idxs.append(0)
                continue
            if filtered_table[i]["Time"] - filtered_table[i-1]["Time"] > max_source_delay:
                start_idxs.append(i)
                end_idxs.append(i)
        end_idxs.append(len(filtered_table))

        contiguous = list(zip(start_idxs, end_idxs))

        for i, cont in enumerate(contiguous):
            retval[s]["Obs{}".format(i)] = {}
            print("---------------")
            print("{}, observation {}\n".format(s, i + 1))
            ft = filtered_table[cont[0]:cont[1]]

            observation_start = ft[0]["Time"]
            observation_end = ft[-1]["Time"]

            print("Source observations:")
            retval[s]["Obs{}".format(i)]["Src"] = []
            for c in range(cont[0], cont[1]):
                print(standard_string(filtered_table[c]["Dir"]))
                retval[s]["Obs{}".format(i)]["Src"].append(standard_string(filtered_table[c]["Dir"]))

            print("")
            print("Calibrator observations:")
            retval[s]["Obs{}".format(i)]["Cal"] = []
            condition1 = np.abs(info["Time"] - observation_start) < max_calibrator_delay
            condition2 = np.abs(info["Time"] - observation_end) < max_calibrator_delay
            condition = condition1 & condition2

            for row in info[condition]:
                if standard_string(row["Source"]) in calibrators:
                    print(standard_string(row["Dir"]))
                    retval[s]["Obs{}".format(i)]["Cal"].append(standard_string(row["Dir"]))

            print("")
            print("---------------\n")
    return retval


def main_inspector(args=None):
    import argparse

    description = ('From a given list of directories, read the relevant information'
                   ' and link observations to calibrators. A single file is read for'
                   ' each directory.')
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("directories", nargs='+',
                        help="Directories to inspect",
                        default=None, type=str)
    parser.add_argument("-g", "--group-by", default=None, type=str, nargs="+")

    args = parser.parse_args(args)

    info = inspect_directories(args.directories)
    info.write('table.csv')
    split_observation_table(info)

    if args.group_by is not None:
        rearranged_info = info.group_by(args.group_by)
        rearranged_info.write('rearranged_table.csv')

