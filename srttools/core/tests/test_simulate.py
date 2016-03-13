from ..simulate import simulate_map, simulate_scan, save_scan
from ..scan import Scan
import numpy as np


def test_sim_scan():
    times, position, shape = simulate_scan()
    save_scan(times, position, np.zeros_like(position), {'Ch0': shape, 'Ch1': shape}, 'output.fits')
    s = Scan('output.fits')

def test_sim_map():
    simulate_map(width_ra=5, width_dec=6.)
