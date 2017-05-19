from __future__ import division, print_function
from ..calibration import CalibratorTable
from ..read_config import read_config
from ..scan import list_scans
import pytest

import os
import glob
import shutil


class TestCalibration(object):
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

        klass.scan_list = \
            list_scans(klass.config['datadir'],
                       klass.config['calibrator_directories'])

        klass.scan_list.sort()

    def test_check_class(self):
        caltable = CalibratorTable()
        caltable.from_scans(self.scan_list)
        with pytest.warns(UserWarning):
            caltable.check_up_to_date()

    def test_calibration(self):
        """Simple calibration from scans."""

        caltable = CalibratorTable()
        caltable.from_scans(self.scan_list)

        caltable.update()
        caltable.Jy_over_counts('Ch0')
        caltable.Jy_over_counts('Ch0', elevation=19)
        caltable.counts_over_Jy('Ch0')

    def test_calibration_write_and_plot(self):
        """Simple calibration from scans."""

        caltable = CalibratorTable()
        caltable.from_scans(self.scan_list)

        caltable.update()
        caltable.write(os.path.join(self.datadir, 'calibrators.hdf5'),
                       path="config", overwrite=True)
        caltable.show()

    @classmethod
    def teardown_class(klass):
        """Clean up the mess."""
        os.unlink('calibration_summary.png')
        os.unlink(os.path.join(klass.datadir, 'calibrators.hdf5'))
        for d in klass.config['list_of_directories']:
            hfiles = \
                glob.glob(os.path.join(klass.config['datadir'], d, '*.hdf5'))
            for h in hfiles:
                os.unlink(h)

            dirs = \
                glob.glob(os.path.join(klass.config['datadir'], d,
                                       '*_scanfit'))
            for dirname in dirs:
                shutil.rmtree(dirname)
