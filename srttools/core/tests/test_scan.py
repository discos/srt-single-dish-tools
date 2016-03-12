# -*- coding: utf-8 -*-

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ..read_config import read_config

import matplotlib.pyplot as plt
import unittest

from ..scan import Scan, ScanSet


class Test1_Scan(unittest.TestCase):
    @classmethod
    def setup_class(klass):
        import os
        global DEBUG_MODE
        DEBUG_MODE = True

        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, '..', '..', 'TEST_DATASET')

        klass.fname = \
            os.path.abspath(
                os.path.join(klass.datadir, '20140603-103246-scicom-3C157',
                             '20140603-103246-scicom-3C157_003_003.fits'))

        klass.config_file = \
            os.path.abspath(os.path.join(klass.curdir, '..', '..',
                                         'TEST_DATASET',
                                         'test_config.ini'))

        read_config(klass.config_file)

    def step_1_scan(self):
        '''Test that data are read.'''

        scan = Scan(self.fname)

        scan.write('scan.hdf5', overwrite=True)

    def step_2_read_scan(self):
        scan = Scan('scan.hdf5')
        plt.ion()
        for col in scan.chan_columns():
            plt.plot(scan['time'], scan[col])
        plt.ioff()
        plt.show()

        return scan

    def test_all(self):
        self.step_1_scan()
        self.step_2_read_scan()
