import numpy as np
import os
from srttools.core.simulate import simulate_scan, save_scan


def gauss(x):
    amp = np.random.normal(100, 1)
    print(amp)
    return amp * np.exp(-x ** 2 /(2*0.1**2))


def gauss_outlier(x):
    return 140 * np.exp(-x ** 2 /(2*0.1**2))


for i in range(12):
    factor = (i // 2) * 2 + 1
    center = factor * 7
    noise = np.random.normal(1, 0.1)
    times, position, shape = simulate_scan(shape=gauss, noise_amplitude=noise)

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

    save_scan(times + (i//4) * 3600 + 57400 * 86400, ra, dec,
              {'Ch0': shape, 'Ch1': shape},
              os.path.join('calibrators', 'calibrator{}.fits'.format(i)),
              scan_type=direction)

# produce outlier!
times, position, shape = simulate_scan(shape=gauss_outlier)
center = 21
ra = position / np.cos(np.radians(center)) + center
dec = np.zeros_like(position) + center

save_scan(times + 57400 * 86400, ra, dec,
          {'Ch0': shape, 'Ch1': shape},
          os.path.join('calibrators', 'calibrator{}.fits'.format(i)),
          scan_type=direction)

