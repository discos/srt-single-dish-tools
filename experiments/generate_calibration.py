import numpy as np
import os
from srttools.core.simulate import simulate_scan, save_scan


def gauss(x):
    return 100 * np.exp(-x ** 2 /(2*0.1**2))


for i in range(12):
    factor = (i // 2) * 2 + 1
    center = factor * 7
    times, position, shape = simulate_scan(shape=gauss)

    if i % 2 == 1:
        position = position[::-1]

    if i % 4 in [0, 1]:
        ra = position / np.cos(np.radians(center)) + center
        dec = np.zeros_like(position) + center
        direction = "RA"
    else:
        dec = position + center
        ra = np.zeros_like(position) + center
        direction = "Dec"

    print(i, direction, center)

    save_scan(times + 57400 * 86400, ra, dec,
              {'Ch0': shape, 'Ch1': shape},
              os.path.join('calibrator_ra', 'calibrator{}.fits'.format(i)))
