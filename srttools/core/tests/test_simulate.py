from ..simulate import simulate_map, simulate_scan, save_scan
import numpy as np
import shutil
from ..io import mkdir_p
import os


def _2d_gauss(x, y):
    return 100 * np.exp(- (x ** 2 + y ** 2) / 8)


class TestSimulate(object):
    @classmethod
    def setup_class(cls):
        cls.outdir = os.path.join('sim')
        cls.emptydir = mkdir_p(os.path.join('sim', 'empty'))
        cls.gaussdir = mkdir_p(os.path.join('sim', 'gauss'))

    def test_sim_scan(self):
        """Test the simulation of a single scan."""
        times, position, shape = simulate_scan()
        save_scan(times, position, np.zeros_like(position), {'Ch0': shape, 'Ch1': shape},
                  os.path.join(self.outdir, 'output.fits'))

    def test_sim_map_empty(self):
        """Test the simulation of an empty map."""
       simulate_map(width_ra=5, width_dec=6., outdir=self.emptydir)

    def test_sim_map_gauss(self):
        """Test the simulation of an map with a large Gaussian PSF."""
        simulate_map(width_ra=5, width_dec=6., count_map=_2d_gauss,
                     outdir=self.gaussdir)

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.outdir)

