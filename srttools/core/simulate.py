"""Functions to simulate scans and maps."""

from __future__ import (absolute_import, division,
                        print_function)

import numpy as np
import numpy.random as ra
import os
from astropy.io import fits
from astropy.table import Table, vstack
from .scan import Scan
from .io import mkdir_p, locations
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
from astropy.time import Time
import astropy.units as u


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
        def shape(x): return 100

    nbins = np.rint(length / speed / dt)

    times = np.arange(nbins) * dt
    # In degrees!
    position = np.arange(-nbins / 2, nbins / 2) / nbins * length / 60

    return times, position + center, shape(position) + \
        ra.normal(0, noise_amplitude, position.shape)


def save_scan(times, ra, dec, channels, filename='out.fits',
              other_columns=None, scan_type=None):
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
    """
    curdir = os.path.abspath(os.path.dirname(__file__))
    template = os.path.abspath(os.path.join(curdir, '..', 'data',
                                            'scan_template.fits'))
    lchdulist = fits.open(template)
    datahdu = lchdulist['DATA TABLE']
    lchdulist[0].header['SOURCE'] = "Dummy"
    lchdulist[0].header['ANTENNA'] = "SRT"
    lchdulist[0].header['HIERARCH RIGHTASCENSION'] = np.radians(np.mean(ra))
    lchdulist[0].header['HIERARCH DECLINATION'] = np.radians(np.mean(dec))
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
    lchdulist.writeto(filename, clobber=True)


def simulate_map(dt=0.04, length_ra=120., length_dec=120., speed=4.,
                 spacing=0.5, count_map=None, noise_amplitude=1.,
                 width_ra=None, width_dec=None, outdir='sim/',
                 baseline="flat", mean_ra=180, mean_dec=70):

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
    """
    import matplotlib.pyplot as plt

    mkdir_p(outdir)
    if count_map is None:
        def count_map(x, y): return 100

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

    times = np.arange(nbins_ra) * dt
    # In degrees!
    position_ra = mean_ra + \
        np.arange(-nbins_ra / 2, nbins_ra / 2) / nbins_ra * length_ra / 60
    position_dec = mean_dec + \
        np.arange(-nbins_dec / 2, nbins_dec / 2) / nbins_dec * length_dec / 60

    if width_dec is None:
        width_dec = length_ra
    if width_ra is None:
        width_ra = length_dec
    # Dec scans
    fig = plt.figure()

    for i_d, start_dec in enumerate(mean_dec + np.arange(-width_dec/2,
                                                         width_dec/2 + spacing,
                                                         spacing)/60):
        m = ra.uniform(mmin, mmax)
        q = ra.uniform(qmin, qmax)
        signs = np.random.choice([-1, 1], nbins_ra)
        stochastic = \
            np.cumsum(signs) * stochastic_amp / np.sqrt(nbins_ra)

        baseline = m * position_ra + q + stochastic
        counts = count_map(position_ra, start_dec) + \
            ra.normal(0, noise_amplitude, position_ra.shape) + \
            baseline

        save_scan(times, position_ra, np.zeros_like(position_ra) + start_dec,
                  {'Ch0': counts, 'Ch1': counts},
                  filename=os.path.join(outdir, 'Ra{}.fits'.format(i_d)))
        plt.plot(position_ra, counts)

    # RA scans
    for i_r, start_ra in enumerate(mean_ra + np.arange(-width_ra / 2,
                                                       width_ra / 2 + spacing,
                                                       spacing) / 60):
        m = ra.uniform(mmin, mmax)
        q = ra.uniform(qmin, qmax)

        signs = np.random.choice([-1, 1], nbins_dec)
        stochastic = \
            np.cumsum(signs) * stochastic_amp / np.sqrt(nbins_dec)

        baseline = m * position_dec + q + stochastic
        counts = count_map(start_ra, position_dec) + \
            ra.normal(0, noise_amplitude, position_dec.shape) + \
            baseline

        save_scan(times, np.zeros_like(position_dec) + start_ra, position_dec,
                  {'Ch0': counts, 'Ch1': counts},
                  filename=os.path.join(outdir, 'Dec{}.fits'.format(i_r)))

        plt.plot(position_dec, counts)

    fig.savefig(os.path.join(outdir, "allscans.png"))
    plt.close(fig)
