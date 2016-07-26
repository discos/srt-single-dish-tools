from ..calibration import CalibratorTable
from ..read_config import read_config
from ..scan import list_scans
import numpy as np
import matplotlib.pyplot as plt
import unittest
from astropy.table import Table
from ..imager import ScanSet
import os
import glob
import shutil


class Test2_Calibration(unittest.TestCase):
    @classmethod
    def setup_class(klass):
        import os
        global DEBUG_MODE
        DEBUG_MODE = True

        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, 'data')

        klass.config_file = \
            os.path.abspath(os.path.join(klass.datadir, "calibrators.ini"))

        klass.config = read_config(klass.config_file)

    def step_1_calibrate(self):
        """Simple calibration from scans."""
        scan_list = \
            list_scans(self.config['datadir'],
                       self.config['calibrator_directories'])

        scan_list.sort()

        caltable = CalibratorTable()
        caltable.from_scans(scan_list)
        caltable.update()
        caltable.write(os.path.join(self.datadir, 'calibrators.hdf5'),
                       path="config", overwrite=True)

    def step_999_cleanup(self):
        """Clean up the mess."""
        os.unlink(os.path.join(self.datadir, 'calibrators.hdf5'))
        for d in self.config['list_of_directories']:
            hfiles = glob.glob(os.path.join(self.config['datadir'], d, '*.hdf5'))
            for h in hfiles:
                os.unlink(h)

            dirs = glob.glob(os.path.join(self.config['datadir'], d, '*_scanfit'))
            for dirname in dirs:
                shutil.rmtree(dirname)

    def test_all(self):
        self.step_1_calibrate()

        self.step_999_cleanup()
