import numpy as np
import shutil
import os
import pytest
from srttools.core.io import mkdir_p
from srttools.core.simulate import simulate_map, simulate_scan, save_scan
from srttools.core.simulate import _default_map_shape


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

    def test_sim_scan_other_columns(self):
        """Test the simulation of a single scan."""
        times, position, shape = simulate_scan()
        save_scan(times, position, np.zeros_like(position),
                  {'Ch0': shape, 'Ch1': shape},
                  os.path.join(self.outdir, 'output.fits'),
                  scan_type="Any",
                  other_columns={'Pippo': np.zeros_like(shape)})

    def test_sim_map_empty(self):
        """Test the simulation of an empty map."""
        simulate_map(width_ra=2, width_dec=2., outdir=self.emptydir)

    def test_sim_map_empty_messy(self):
        """Test the simulation of an empty map."""
        simulate_map(width_ra=2, width_dec=2., outdir=self.emptydir,
                     baseline='messy')

    def test_sim_map_empty_slope(self):
        """Test the simulation of an empty map."""
        simulate_map(width_ra=2, width_dec=2., outdir=self.emptydir,
                     baseline='slope')

    def test_raises_wrong_map_shape(self):
        with pytest.raises(ValueError):
            _default_map_shape(np.zeros((3, 4)), np.ones((3, 6)))

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.outdir)
