from .io import read_data
import glob
from .read_config import read_config


class Scan():
    '''Class containing a single scan'''
    def __init__(self, fname):
        self.table = read_data(fname)
        for col in self.table.columns:
            setattr(self, col, self.table.field(col))

    def baseline_subtract():
        pass


class ScanSet():
    def __init__(self):
        config = read_config


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


