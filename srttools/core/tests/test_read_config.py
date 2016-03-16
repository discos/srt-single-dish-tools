# -*- coding: utf-8 -*-

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np

from ..read_config import read_config


def test_read_config():
    """Test that config file are read."""
    import os
    curdir = os.path.abspath(os.path.dirname(__file__))
    datadir = os.path.join(curdir, 'data')

    fname = os.path.join(datadir, 'test_config.ini')

    config = read_config(fname)

    np.testing.assert_almost_equal(config['pixel_size'], np.radians(0.5 / 60))
    assert config['interpolation'] == 'spline'



def test_read_incomplete_config():
    """Test that config file are read."""
    import os
    curdir = os.path.abspath(os.path.dirname(__file__))
    datadir = os.path.join(curdir, 'data')

    fname = os.path.join(datadir, 'test_config_incomplete.ini')

    config = read_config(fname)

    np.testing.assert_almost_equal(config['pixel_size'], np.radians(1 / 60))
    assert config['interpolation'] == 'linear'

