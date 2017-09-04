# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function)

import numpy as np

from srttools.read_config import read_config
import os
import astropy.units as u


class TestConfig(object):
    @classmethod
    def setup_class(cls):
        cls.curdir = os.path.abspath(os.path.dirname(__file__))
        cls.datadir = os.path.join(cls.curdir, 'data')

    def test_read_config(self):
        """Test that config file are read."""

        fname = os.path.join(self.datadir, 'test_config.ini')

        config = read_config(fname)

        np.testing.assert_almost_equal(config['pixel_size'].to(u.rad).value,
                                       np.radians(0.5 / 60))
        assert config['interpolation'] == 'spline'

    def test_read_incomplete_config(self):
        """Test that config file are read."""
        fname = os.path.join(self.datadir, 'test_config_incomplete.ini')

        config = read_config(fname)

        np.testing.assert_almost_equal(config['pixel_size'].to(u.rad).value,
                                       np.radians(1 / 60))
        assert config['interpolation'] == 'linear'
