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

    assert np.all(config['npix'] == [64, 64])
    assert config['interpolation'] == 'spline'



def test_read_incomplete_config():
    """Test that config file are read."""
    import os
    curdir = os.path.abspath(os.path.dirname(__file__))
    datadir = os.path.join(curdir, 'data')

    fname = os.path.join(datadir, 'test_config_incomplete.ini')

    config = read_config(fname)

    assert np.all(config['npix'] == [32, 32])
    assert config['interpolation'] == 'linear'

