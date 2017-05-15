# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function)

from ..read_config import read_config

import pytest

from ..scan import Scan
from ..io import print_obs_info_fitszilla
import os


class Test1_Scan(object):
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
        h5file = klass.fname.replace('.fits', '.hdf5')
        if os.path.exists(h5file):
            os.unlink(h5file)

        klass.config_file = \
            os.path.abspath(os.path.join(klass.datadir, 'test_config.ini'))
        print(klass.config_file)

        read_config(klass.config_file)

    def test_print_info(self, capsys):
        print_obs_info_fitszilla(self.fname)
        out, err = capsys.readouterr()
        assert 'bandwidth' in out.lower()

    def test_scan(self):
        '''Test that data are read.'''

        scan = Scan(self.fname)

        scan.write('scan.hdf5', overwrite=True)

    @classmethod
    def teardown_class(klass):
        """Cleanup."""
        os.unlink('scan.hdf5')


class Test2_Scan(object):
    @classmethod
    def setup_class(klass):
        import os
        global DEBUG_MODE
        DEBUG_MODE = True

        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, 'data')

        klass.fname = \
            os.path.abspath(
                os.path.join(klass.datadir, 'spectrum',
                             'roach_template.fits'))

        h5file = klass.fname.replace('.fits', '.hdf5')
        if os.path.exists(h5file):
            os.unlink(h5file)
        klass.config_file = \
            os.path.abspath(os.path.join(klass.datadir, 'spectrum.ini'))
        print(klass.config_file)

        read_config(klass.config_file)

    def test_scan(self):
        '''Test that data are read.'''

        scan = Scan(self.fname)

        scan.write('scan.hdf5', overwrite=True)
        scan.baseline_subtract('rough')

    @classmethod
    def teardown_class(klass):
        """Cleanup."""
        os.unlink('scan.hdf5')

