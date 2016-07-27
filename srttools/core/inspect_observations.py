"""From a given list of directories, read the relevant information and link observations to calibrators."""

from __future__ import (absolute_import, division,
                        print_function)
import os
import glob
import warnings
import numpy as np
from astropy.table import Table, Column
from .io import read_data


def inspect_directories(directories):
    info = Table()
    names = ["Dir", "File", "Source", "Receiver", "Backend",
             "Time", "Frequency", "Bandwidth"]

    dtype = ['S200', 'S200', 'S200', 'S200', 'S200',
             np.double, np.float, np.float]

    for n, d in zip(names, dtype):
        if n not in info.keys():
            info.add_column(Column(name=n, dtype=d))

    for d in directories:
        fits_files = glob.glob(os.path.join(d, '*.fits'))

        for f in fits_files:
            print("Reading {}".format(f), end="\r")
            try:
                data = read_data(f)
                backend = data.meta['backend']
                receiver = data.meta['receiver']
                frequency = data["Ch0"].meta['frequency']
                bandwidth = data["Ch0"].meta['bandwidth']
                source = data.meta['SOURCE']
                time = data[0]['time']

                info.add_row([d, f, source, receiver, backend,
                              time, frequency, bandwidth])
            except:
                pass
    print(info)


def main_inspector(args=None):
    import argparse

    description = ('From a given list of directories, read the relevant information'
                   ' and link observations to calibrators.')
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("directories", nargs='+',
                        help="Directories to inspect",
                        default=None, type=str)

    args = parser.parse_args(args)

    inspect_directories(args.directories)