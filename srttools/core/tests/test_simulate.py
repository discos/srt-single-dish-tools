from ..simulate import simulate_map, simulate_scan, save_scan
from ..scan import Scan
import numpy as np


def _2d_gauss(x, y):
    return 100 * np.exp(- (x ** 2 + y ** 2) / 8)

def test_sim_scan():
    times, position, shape = simulate_scan()
    save_scan(times, position, np.zeros_like(position), {'Ch0': shape, 'Ch1': shape}, 'output.fits')
    s = Scan('output.fits')

def test_sim_map_empty():
    simulate_map(width_ra=5, width_dec=6.)

def test_sim_map_gauss():
    simulate_map(width_ra=5, width_dec=6., count_map=_2d_gauss, outdir='gauss/')

