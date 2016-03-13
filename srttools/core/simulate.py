"""Functions to simulate scans and maps."""

import numpy as np
import numpy.random as ra
import os
from astropy.io import fits
from astropy.table import Table, vstack
from .scan import Scan

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


def simulate_map(dt=0.04, length_ra=120., length_dec=120., speed=4., spacing=0.5, shape=None, noise_amplitude=1.):
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
    if shape is None:
        shape = lambda x : 100

    ra_scans = {}
    dec_scans = {}
    for start_dec in np.arange(0, length_dec + spacing, spacing):
        ra_scans[start_dec]


def test_save_scan():
    times, position, shape = simulate_scan()
    save_scan(times, position, np.zeros_like(position), {'Ch0': shape, 'Ch1': shape}, 'output.fits')
    s = Scan('output.fits')
    print(s)