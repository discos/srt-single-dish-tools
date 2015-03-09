from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
from .io import read_data
import glob
from .read_config import read_config
import os


class Scan():
    '''Class containing a single scan'''
    def __init__(self, fname):
        self.filename = fname
        self.table = read_data(fname)
        for col in self.table.columns:
            setattr(self, col, self.table.field(col))

    def baseline_subtract(self):
        pass

    def __repr__(self):
        reprstring = '\n\n----Scan from file {} ----\n'.format(self.filename)
        reprstring += repr(self.table)
        return reprstring


class ScanSet():
    def __init__(self, config_file=None):
        self.config = read_config(config_file)
        self.scans = []
        datadir = self.config['datadir']
        dirlist = self.config['list_of_directories']
        for d in dirlist:
            for f in glob.glob(os.path.join(datadir, d, '*')):
                self.scans.append(Scan(f))

    def __repr__(self):
        reprstring = ''
        for s in self.scans:
            reprstring += repr(s)

        return reprstring


def test_scan():
    '''Test that data are read.'''
    import os
    import matplotlib.pyplot as plt
    curdir = os.path.abspath(os.path.dirname(__file__))
    datadir = os.path.join(curdir, '..', '..', 'TEST_DATASET')

    fname = os.path.join(datadir, '20140603-103246-scicom-3C157',
                         '20140603-103246-scicom-3C157_003_003.fits')

    scan = Scan(fname)

    for i in range(2):
        plt.plot(scan.time, scan.data[:, i])
    plt.show()


def test_scanset():
    '''Test that sets of data are read.'''
    import os
    curdir = os.path.abspath(os.path.dirname(__file__))
    config = os.path.join(curdir, '..', '..', 'TEST_DATASET',
                          'test_config.ini')

    scanset = ScanSet(config)

    print(scanset)
