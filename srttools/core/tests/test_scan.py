# -*- coding: utf-8 -*-

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ..read_config import read_config

import unittest

from ..scan import Scan
import os


class Test1_Scan(unittest.TestCase):
    @classmethod
    def setup_class(klass):
        import os
        global DEBUG_MODE
        DEBUG_MODE = True

        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, 'data')

        klass.fname = \
            os.path.abspath(
                os.path.join(klass.datadir, 'gauss_dec',
                             'Dec0.fits'))

        klass.config_file = \
            os.path.abspath(os.path.join(klass.datadir, 'test_config.ini'))
        print(klass.config_file)

        read_config(klass.config_file)

    def step_1_scan(self):
        '''Test that data are read.'''

        scan = Scan(self.fname)

        scan.write('scan.hdf5', overwrite=True)

    def cleanup(self):
        """Cleanup."""
        os.unlink('scan.hdf5')

    def test_all(self):
        self.step_1_scan()
        self.cleanup()
