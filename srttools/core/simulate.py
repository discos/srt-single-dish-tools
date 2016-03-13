"""Functions to simulate scans and maps."""

import numpy as np
import numpy.random as ra
import os
from astropy.io import fits
from astropy.table import Table, vstack
from .scan import Scan
from .io import mkdir_p

def simulate_scan(dt=0.04, length=120., speed=4., shape=None, noise_amplitude=1.):
    """Simulate a scan.

    Parameters
    ----------
    dt : float
        The integration time in seconds
    length : float
        Length of the scan in arcminutes
    speed : float
        Speed of the scan in arcminutes / second
    shape : function
        Function that describes the shape of the scan. If None, a
        constant scan is assumed. The zero point of the scan is in the
        *center* of it
    noise_amplitude : float
        Noise level in counts
    """
    if shape is None:
        shape = lambda x : 100
    nbins = np.rint(length / speed / dt)

    times = np.arange(nbins) * dt
    position = np.arange(-nbins / 2, nbins / 2) / nbins * length / 60  # In degrees!

    return times, position, shape(position) + ra.normal(0, noise_amplitude, position.shape)


def save_scan(times, ra, dec, channels, filename='out.fits', other_columns=None):
    """Save a simulated scan in fitszilla format.

    Parameters
    ----------
    times : iterable
        times corresponding to each bin center
    ra : iterable
        RA corresponding to each bin center
    dec : iterable
        Dec corresponding to each bin center
    channels : {'Ch0': array([...]), 'Ch1': array([...]), ...}
        Dictionary containing the count array. Keys represent the name of the channel
    filename : str
        Output file name
    """
    template = os.path.abspath(os.path.join(os.getcwd(), '..', 'data', 'scan_template.fits'))
    lchdulist = fits.open(template)
    datahdu = lchdulist['DATA TABLE']
    data_table_data = Table(datahdu.data)

    newtable = Table(names=['time', 'raj2000', 'decj2000'], data=[times, np.radians(ra), np.radians(dec)])
    for ch in channels.keys():
        newtable[ch] = channels[ch]
    if other_columns is None:
        other_columns = {}
    for col in other_columns.keys():
        newtable[col] = other_columns[col]

    data_table_data = vstack([data_table_data, newtable])

    nrows = len(data_table_data)

    hdu = fits.BinTableHDU.from_columns(datahdu.data.columns, nrows=nrows)
    for colname in datahdu.data.columns.names:
        hdu.data[colname][:] = data_table_data[colname]

    datahdu.data = hdu.data
    # print(datahdu)
    # lchdulist['DATA TABLE'].name = 'TMP'
    # lchdulist.append(datahdu)
    lchdulist.writeto(filename, clobber=True)


def simulate_map(dt=0.04, length_ra=120., length_dec=120., speed=4., spacing=0.5, count_map=None, noise_amplitude=1.,
                 width_ra=None, width_dec=None, outdir='sim/'):
    """Simulate a map.

    Parameters
    ----------
    dt : float
        The integration time in seconds
    length : float
        Length of the scan in arcminutes
    speed : float
        Speed of the scan in arcminutes / second
    shape : function
        Function that describes the shape of the scan. If None, a
        constant scan is assumed. The zero point of the scan is in the
        *center* of it
    noise_amplitude : float
        Noise level in counts
    """

    mkdir_p(outdir)
    if count_map is None:
        count_map = lambda x, y : 100

    nbins_ra = np.rint(length_ra / speed / dt)
    nbins_dec = np.rint(length_dec / speed / dt)

    times = np.arange(nbins_ra) * dt
    position_ra = np.arange(-nbins_ra / 2, nbins_ra / 2) / nbins_ra * length_ra / 60  # In degrees!
    position_dec = np.arange(-nbins_dec / 2, nbins_dec / 2) / nbins_dec * length_dec / 60  # In degrees!

    if width_dec is None:
        width_dec = length_ra
    if width_ra is None:
        width_ra = length_dec
    # Dec scans
    print(nbins_ra, nbins_dec)
    for i_d, start_dec in enumerate(np.arange(-width_dec / 2 , width_dec / 2 + spacing, spacing) / 60):
        counts = count_map(position_ra, start_dec) + ra.normal(0, noise_amplitude, position_ra.shape)
        save_scan(times, position_ra, np.zeros_like(position_ra) + start_dec, {'Ch0': counts, 'Ch1': counts},
                  filename=os.path.join(outdir, 'Ra{}.fits'.format(i_d)))

    print(i_d)
    # RA scans
    for i_r, start_ra in enumerate(np.arange(-width_ra / 2 , width_ra / 2 + spacing, spacing) / 60):
        counts = count_map(start_ra, position_dec) + ra.normal(0, noise_amplitude, position_dec.shape)
        save_scan(times, np.zeros_like(position_dec) + start_ra, position_dec, {'Ch0': counts, 'Ch1': counts},
                  filename=os.path.join(outdir, 'Dec{}.fits'.format(i_r)))
    print(i_r)
