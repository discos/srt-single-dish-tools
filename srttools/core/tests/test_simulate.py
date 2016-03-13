from ..simulate import simulate_map, simulate_scan, save_scan
import numpy as np
import shutil
from ..io import mkdir_p
import os



def _2d_gauss(x, y):
    return 100 * np.exp(- (x ** 2 + y ** 2) / 8)

def test_sim_scan():
    """Test the simulation of a single scan."""
    outdir = os.path.join('sim')
    mkdir_p(outdir)
    times, position, shape = simulate_scan()
    save_scan(times, position, np.zeros_like(position), {'Ch0': shape, 'Ch1': shape},
              os.path.join(outdir, 'output.fits'))
    shutil.rmtree(outdir)

def test_sim_map_empty():
    """Test the simulation of an empty map."""
    outdir = os.path.join('sim', 'empty')
    simulate_map(width_ra=5, width_dec=6., outdir=outdir)
    shutil.rmtree(outdir)

def test_sim_map_gauss():
    """Test the simulation of an map with a large Gaussian PSF."""
    outdir = os.path.join('sim', 'gauss')
    simulate_map(width_ra=5, width_dec=6., count_map=_2d_gauss, outdir=outdir)
    shutil.rmtree(outdir)

