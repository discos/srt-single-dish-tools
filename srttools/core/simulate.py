"""Functions to simulate scans and maps."""

from __future__ import (absolute_import, division,
                        print_function)

import numpy as np
import numpy.random as ra
import os
from astropy.io import fits
from astropy.table import Table, vstack
from .io import mkdir_p, locations
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
import six
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(x):
        return x


def _default_flat_shape(x):
    """A flat shape.

    Examples
    --------
    >>> _default_flat_shape(4314)
    100.0
    >>> _default_flat_shape(np.arange(3))
    array([ 100.,  100.,  100.])
    """
    return 100 + np.zeros(np.asarray(x).shape)


def _default_map_shape(x, y):
    """A flat map shape.

    Examples
    --------
    >>> _default_map_shape(4314, 234)
    100.0
    >>> _default_map_shape(np.zeros((3, 4)), np.ones((3, 4)))
    array([[ 100.,  100.,  100.,  100.],
           [ 100.,  100.,  100.,  100.],
           [ 100.,  100.,  100.,  100.]])
    """
    x = np.asarray(x)
    y = np.asarray(y)
    # It will raise a ValueError when x and y are not compatible
    return 100 + np.zeros_like(y) * np.zeros_like(x)


def simulate_scan(dt=0.04, length=120., speed=4., shape=None,
                  noise_amplitude=1., center=0.):
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
    center : float
        Center coordinate in degrees
    """
    if shape is None:
        shape = _default_flat_shape

    nbins = np.rint(length / speed / dt)

    times = np.arange(nbins) * dt
    # In degrees!
    position = np.arange(-nbins / 2, nbins / 2) / nbins * length / 60

    return times, position + center, shape(position) + \
        ra.normal(0, noise_amplitude, position.shape)


def save_scan(times, ra, dec, channels, filename='out.fits',
              other_columns=None, scan_type=None, src_ra=None, src_dec=None,
              srcname='Dummy'):
    """Save a simulated scan in fitszilla format.

    Parameters
    ----------
    times : iterable
        times corresponding to each bin center, in seconds
    ra : iterable
        RA corresponding to each bin center
    dec : iterable
        Dec corresponding to each bin center
    channels : {'Ch0': array([...]), 'Ch1': array([...]), ...}
        Dictionary containing the count array. Keys represent the name of the
        channel
    filename : str
        Output file name
    srcname : str
        Name of the source
    """
    if src_ra is None:
        src_ra = np.mean(ra)
    if src_dec is None:
        src_dec = np.mean(dec)

    curdir = os.path.abspath(os.path.dirname(__file__))
    template = os.path.abspath(os.path.join(curdir, '..', 'data',
                                            'scan_template.fits'))
    lchdulist = fits.open(template)
    datahdu = lchdulist['DATA TABLE']
    lchdulist[0].header['SOURCE'] = "Dummy"
    lchdulist[0].header['ANTENNA'] = "SRT"
    lchdulist[0].header['HIERARCH RIGHTASCENSION'] = np.radians(src_ra)
    lchdulist[0].header['HIERARCH DECLINATION'] = np.radians(src_dec)
    if scan_type is not None:
        lchdulist[0].header['HIERARCH SubScanType'] = scan_type

    data_table_data = Table(datahdu.data)

    obstimes = Time((times / 86400 + 57000) * u.day, format='mjd', scale='utc')

    coords = SkyCoord(ra, dec, unit=u.degree, location=locations['srt'],
                      obstime=obstimes)

    altaz = coords.altaz
    el = altaz.alt.rad
    az = altaz.az.rad
    newtable = Table(names=['time', 'raj2000', 'decj2000', "el", "az"],
                     data=[obstimes.value, np.radians(ra), np.radians(dec),
                           el, az])

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
    lchdulist[0].header['SOURCE'] = srcname
    lchdulist.writeto(filename, clobber=True)


def simulate_map(dt=0.04, length_ra=120., length_dec=120., speed=4.,
                 spacing=0.5, count_map=None, noise_amplitude=1.,
                 width_ra=None, width_dec=None, outdir='sim/',
                 baseline="flat", mean_ra=180, mean_dec=70,
                 srcname='Dummy'):

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
    spacing : float
        Spacing between scans, in arcminutes
    baseline : str
        "flat", "slope" (linearly increasing/decreasing) or "messy"
        (random walk)
    count_map : function
        Flux distribution function, centered on zero
    outdir : str or iterable (str, str)
        If a single string, put all files in that directory; if two strings,
        put RA and DEC scans in the two directories.
    """
    import matplotlib.pyplot as plt

    if isinstance(outdir, six.string_types):
        outdir = (outdir, outdir)
    outdir_ra = outdir[0]
    outdir_dec = outdir[1]

    mkdir_p(outdir_ra)
    mkdir_p(outdir_dec)

    if count_map is None:
        count_map = _default_map_shape

    if baseline == "flat":
        mmin = mmax = 0
        qmin = qmax = 0
        stochastic_amp = 0
    elif baseline == "slope":
        mmin, mmax = -5, 5
        qmin, qmax = 0, 150
        stochastic_amp = 0
    elif baseline == "messy":
        mmin, mmax = 0, 0
        qmin, qmax = 0, 0
        stochastic_amp = 20

    nbins_ra = np.int(np.rint(length_ra / speed / dt))
    nbins_dec = np.int(np.rint(length_dec / speed / dt))

    times_ra = np.arange(nbins_ra) * dt
    times_dec = np.arange(nbins_dec) * dt

    ra_array = np.arange(-nbins_ra / 2,
                         nbins_ra / 2) / nbins_ra * length_ra / 60
    dec_array = np.arange(-nbins_dec / 2,
                          nbins_dec / 2) / nbins_dec * length_dec / 60
    # In degrees!
    if width_dec is None:
        width_dec = length_dec
    if width_ra is None:
        width_ra = length_ra
    # Dec scans
    fig = plt.figure()

    delta_decs = np.arange(-width_dec/2, width_dec/2 + spacing, spacing)/60
    print("Simulating dec scans...")
    for i_d, delta_dec in enumerate(tqdm(delta_decs)):

        start_dec = mean_dec + delta_dec
        m = ra.uniform(mmin, mmax)
        q = ra.uniform(qmin, qmax)
        signs = np.random.choice([-1, 1], nbins_ra)
        stochastic = \
            np.cumsum(signs) * stochastic_amp / np.sqrt(nbins_ra)

        baseline = m * ra_array + q + stochastic
        counts = count_map(ra_array, delta_dec) + \
            ra.normal(0, noise_amplitude, ra_array.shape) + \
            baseline

        actual_ra = mean_ra + ra_array / np.cos(np.radians(start_dec))

        save_scan(times_ra, actual_ra, np.zeros_like(actual_ra) + start_dec,
                  {'Ch0': counts, 'Ch1': counts},
                  filename=os.path.join(outdir_ra, 'Ra{}.fits'.format(i_d)),
                  src_ra=mean_ra, src_dec=mean_dec, srcname=srcname)
        plt.plot(ra_array, counts)
    fig.savefig(os.path.join(outdir_ra, "allscans_ra.png"))
    plt.close(fig)

    fig = plt.figure()
    delta_ras = np.arange(-width_ra / 2, width_ra / 2 + spacing,
                          spacing) / 60
    print("Simulating RA scans...")
    # RA scans
    for i_r, delta_ra in enumerate(tqdm(delta_ras)):
        start_ra = delta_ra / np.cos(np.radians(mean_dec)) + mean_ra
        m = ra.uniform(mmin, mmax)
        q = ra.uniform(qmin, qmax)

        signs = np.random.choice([-1, 1], nbins_dec)
        stochastic = \
            np.cumsum(signs) * stochastic_amp / np.sqrt(nbins_dec)

        baseline = m * dec_array + q + stochastic
        counts = count_map(delta_ra, dec_array) + \
            ra.normal(0, noise_amplitude, dec_array.shape) + \
            baseline

        save_scan(times_dec, np.zeros_like(dec_array) + start_ra,
                  dec_array + mean_dec,
                  {'Ch0': counts, 'Ch1': counts},
                  filename=os.path.join(outdir_dec, 'Dec{}.fits'.format(i_r)),
                  src_ra=mean_ra, src_dec=mean_dec, srcname=srcname)

        plt.plot(dec_array, counts)

    fig.savefig(os.path.join(outdir_dec, "allscans_dec.png"))
    plt.close(fig)
