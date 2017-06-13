from ..simulate import simulate_map, simulate_scan, save_scan
from ..simulate import _default_map_shape
import numpy as np
import shutil
from ..io import mkdir_p
import os
import pytest


class TestSimulate(object):
    @classmethod
    def setup_class(cls):
        cls.outdir = os.path.join('sim')
        cls.emptydir = os.path.join('sim', 'empty')
        for d in [cls.emptydir]:
            mkdir_p(d)

    def test_sim_scan(self):
        """Test the simulation of a single scan."""
        times, position, shape = simulate_scan()
        save_scan(times, position, np.zeros_like(position),
                  {'Ch0': shape, 'Ch1': shape},
                  os.path.join(self.outdir, 'output.fits'))

    def test_sim_map_empty(self):
        """Test the simulation of an empty map."""
        simulate_map(width_ra=2, width_dec=2., outdir=self.emptydir)

    def test_raises_wrong_map_shape(self):
        with pytest.raises(ValueError):
            _default_map_shape(np.zeros((3, 4)), np.ones((3, 6)))

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.outdir)
