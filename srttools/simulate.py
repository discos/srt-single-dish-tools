"""Functions to simulate scans and maps."""

from __future__ import (absolute_import, division,
                        print_function)

import numpy as np
import numpy.random as ra
import os
from astropy.io import fits
from astropy.table import Table, vstack
import astropy.units as u
import six
import collections

from .io import mkdir_p, locations
from .utils import tqdm
from astropy.coordinates import SkyCoord
from astropy.time import Time
try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

__all__ = ["simulate_scan", "save_scan", "simulate_map"]


def _is_number(x):
    """"Test if a string or other is a number

    Examples
    --------
    >>> _is_number('3')
    True
    >>> _is_number(3.)
    True
    >>> _is_number('a')
    False
    """
    try:
        float(x)
        return True
    except ValueError:
        return False


def _default_flat_shape(x):
    """A flat shape.

    Examples
    --------
    >>> _default_flat_shape(4314)
    100.0
    >>> np.allclose(_default_flat_shape(np.arange(3)),
    ...             np.array([100., 100., 100.]))
    True
    """
    return 100 + np.zeros(np.asarray(x).shape)


def _2d_gauss(x, y, sigma=2.5 / 60.):
    """A Gaussian beam"""
    return np.exp(-(x ** 2 + y ** 2) / (2 * sigma**2))


def calibrator_scan_func(x):
    return 100 * _2d_gauss(x, 0, sigma=2.5 / 60)


def sim_crossscans(ncross, caldir, scan_func=calibrator_scan_func,
                   srcname='DummyCal', channel_ratio=0.8, baseline="flat"):
    src_ra = 185
    src_dec = 75
    speed = 2.  # arcmin/s
    dt = 0.04
    dtheta = speed * dt
    length = 4 / dtheta

    timedelta = 0

    for i in tqdm(range(ncross)):
        times, ras, scan0 = \
            simulate_scan(dt=dt, length=length, speed=speed, shape=scan_func,
                          noise_amplitude=0.2, center=0,
                          baseline=baseline)
        _, _, scan1 = \
            simulate_scan(dt=dt, length=length, speed=speed, shape=scan_func,
                          noise_amplitude=0.2, center=0,
                          baseline=baseline)

        ras = ras / np.cos(np.radians(src_dec)) + src_ra
        if i % 2 != 0:
            ras = ras[::-1]

        decs = np.zeros_like(ras) + src_dec

        save_scan(times + timedelta, ras, decs,
                  {'Ch0': scan0, 'Ch1': scan1 * channel_ratio},
                  filename=os.path.join(caldir, '{}_Ra.fits'.format(i)),
                  src_ra=src_ra, src_dec=src_dec, srcname=srcname,
                  counts_to_K=(0.03, 0.03 / channel_ratio))
        timedelta += times[-1] + 1

        times, decs, scan0 = \
            simulate_scan(dt=dt, length=length, speed=speed, shape=scan_func,
                          noise_amplitude=0.2, center=src_dec,
                          baseline=baseline)
        _, _, scan1 = \
            simulate_scan(dt=dt, length=length, speed=speed, shape=scan_func,
                          noise_amplitude=0.2, center=src_dec,
                          baseline=baseline)

        if i % 2 != 0:
            decs = decs[::-1]

        ras = np.zeros_like(decs) + src_ra

        save_scan(times + timedelta, ras, decs,
                  {'Ch0': scan0, 'Ch1': scan1 * channel_ratio},
                  filename=os.path.join(caldir, '{}_Dec.fits'.format(i)),
                  src_ra=src_ra, src_dec=src_dec, srcname=srcname,
                  counts_to_K=(0.03, 0.03 / channel_ratio))
        timedelta += times[-1] + 1


def _default_map_shape(x, y):
    """A flat map shape.

    Examples
    --------
    >>> _default_map_shape(4314, 234)
    100
    >>> res = np.array([[ 100.,  100.,  100.,  100.],
    ...                 [ 100.,  100.,  100.,  100.],
    ...                 [ 100.,  100.,  100.,  100.]])
    >>> np.allclose(_default_map_shape(np.zeros((3, 4)), np.ones((3, 4))), res)
    True
    """
    x = np.asarray(x)
    y = np.asarray(y)
    # It will raise a ValueError when x and y are not compatible
    return 100 + np.zeros_like(y) * np.zeros_like(x)


def simulate_scan(dt=0.04, length=120., speed=4., shape=None,
                  noise_amplitude=1., center=0., baseline="flat"):
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
    baseline : str
        "flat", "slope" (linearly increasing/decreasing), "messy"
        (random walk) or a number (which gives an amplitude to the random-walk
        baseline, that is 20 for "messy").
    """
    if shape is None:
        shape = _default_flat_shape

    nbins = np.rint(length / speed / dt)

    times = np.arange(nbins) * dt
    # In degrees!
    position = np.arange(-nbins / 2, nbins / 2) / nbins * length / 60

    scan_baseline = _create_baseline(position, baseline)

    return times, position + center, shape(position) + \
        ra.normal(0, noise_amplitude, position.shape) + scan_baseline


def save_scan(times, ra, dec, channels, filename='out.fits',
              other_columns=None, scan_type=None, src_ra=None, src_dec=None,
              srcname='Dummy', counts_to_K=0.03):
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
    counts_to_K : float, array or dict
        Conversion factor between counts and K. If array, it has to be the same
        length as channels.keys()
    """
    if src_ra is None:
        src_ra = np.mean(ra)
    if src_dec is None:
        src_dec = np.mean(dec)
    # If it's a single value, make it into a list
    if not isinstance(counts_to_K, collections.Iterable):
        counts_to_K = counts_to_K * np.ones(len(list(channels.keys())))
    # If it's a list, make it into a dict
    if not hasattr(counts_to_K, 'keys'):
        counts_to_K = dict([(ch, counts_to_K[i])
                            for i, ch in enumerate(channels.keys())])

    curdir = os.path.abspath(os.path.dirname(__file__))
    template = os.path.abspath(os.path.join(curdir, 'data',
                                            'scan_template.fits'))
    lchdulist = fits.open(template)
    datahdu = lchdulist['DATA TABLE']
    temphdu = lchdulist['ANTENNA TEMP TABLE']
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

    temptable = Table()
    for ch in channels.keys():
        temptable[ch] = newtable[ch] * counts_to_K[ch]

    thdu = fits.BinTableHDU.from_columns(temphdu.data.columns, nrows=nrows)
    for colname in temphdu.data.columns.names:
        thdu.data[colname][:] = temptable[colname]

    temphdu.data = thdu.data

    lchdulist[0].header['SOURCE'] = srcname
    lchdulist.writeto(filename, overwrite=True)
    lchdulist.close()


def _create_baseline(x, baseline_kind="flat"):

    if baseline_kind == "flat":
        mmin = mmax = 0
        qmin = qmax = 0
        stochastic_amp = 0
    elif baseline_kind == "slope":
        mmin, mmax = -5, 5
        qmin, qmax = 0, 150
        stochastic_amp = 0
    elif baseline_kind == "messy":
        mmin, mmax = 0, 0
        qmin, qmax = 0, 0
        stochastic_amp = 20
    elif _is_number(baseline_kind):
        mmin, mmax = 0, 0
        qmin, qmax = 0, 0
        stochastic_amp = float(baseline_kind)
    else:
        raise ValueError("baseline has to be 'flat', 'slope', 'messy' or a "
                         "number")

    n = len(x)
    m = ra.uniform(mmin, mmax)
    q = ra.uniform(qmin, qmax)
    signs = np.random.choice([-1, 1], n)

    stochastic = \
        np.cumsum(signs) * stochastic_amp / np.sqrt(n)

    baseline = m * x + q

    return baseline + stochastic


def simulate_map(dt=0.04, length_ra=120., length_dec=120., speed=4.,
                 spacing=0.5, count_map=None, noise_amplitude=1.,
                 width_ra=None, width_dec=None, outdir='sim/',
                 baseline="flat", mean_ra=180, mean_dec=70,
                 srcname='Dummy', channel_ratio=1):

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
        "flat", "slope" (linearly increasing/decreasing), "messy"
        (random walk) or a number (which gives an amplitude to the random-walk
        baseline, that is 20 for "messy").
    count_map : function
        Flux distribution function, centered on zero
    outdir : str or iterable (str, str)
        If a single string, put all files in that directory; if two strings,
        put RA and DEC scans in the two directories.
    channel_ratio : float
        Ratio between the counts in the two channels
    """

    if isinstance(outdir, six.string_types):
        outdir = (outdir, outdir)
    outdir_ra = outdir[0]
    outdir_dec = outdir[1]

    mkdir_p(outdir_ra)
    mkdir_p(outdir_dec)

    if count_map is None:
        count_map = _default_map_shape

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
    if HAS_MPL:
        fig = plt.figure()

    delta_decs = np.arange(-width_dec/2, width_dec/2 + spacing, spacing)/60
    print("Simulating dec scans...")
    for i_d, delta_dec in enumerate(tqdm(delta_decs)):

        start_dec = mean_dec + delta_dec

        counts_clean = count_map(ra_array, delta_dec)

        baseline0 = _create_baseline(ra_array, baseline)
        baseline1 = _create_baseline(ra_array, baseline)

        counts0 = counts_clean + \
            ra.normal(0, noise_amplitude, ra_array.shape) + baseline0
        counts1 = counts_clean + \
            ra.normal(0, noise_amplitude, ra_array.shape) + baseline1

        actual_ra = mean_ra + ra_array / np.cos(np.radians(start_dec))

        if i_d % 2 != 0:
            actual_ra = actual_ra[::-1]
        save_scan(times_ra, actual_ra, np.zeros_like(actual_ra) + start_dec,
                  {'Ch0': counts0, 'Ch1': counts1 * channel_ratio},
                  filename=os.path.join(outdir_ra, 'Ra{}.fits'.format(i_d)),
                  src_ra=mean_ra, src_dec=mean_dec, srcname=srcname,
                  counts_to_K=(0.03, 0.03 / channel_ratio))
        if HAS_MPL:
            plt.plot(ra_array, counts0)
            plt.plot(ra_array, counts1)

    if HAS_MPL:
        fig.savefig(os.path.join(outdir_ra, "allscans_ra.png"))
        plt.close(fig)

        fig = plt.figure()
    delta_ras = np.arange(-width_ra / 2, width_ra / 2 + spacing,
                          spacing) / 60
    print("Simulating RA scans...")
    # RA scans
    for i_r, delta_ra in enumerate(tqdm(delta_ras)):
        start_ra = delta_ra / np.cos(np.radians(mean_dec)) + mean_ra

        counts_clean = count_map(delta_ra, dec_array)

        baseline0 = _create_baseline(dec_array, baseline)
        baseline1 = _create_baseline(dec_array, baseline)

        counts0 = counts_clean + \
            ra.normal(0, noise_amplitude, ra_array.shape) + baseline0
        counts1 = counts_clean + \
            ra.normal(0, noise_amplitude, ra_array.shape) + baseline1

        if i_r % 2 != 0:
            dec_array = dec_array[::-1]
        save_scan(times_dec, np.zeros_like(dec_array) + start_ra,
                  dec_array + mean_dec,
                  {'Ch0': counts0, 'Ch1': counts1 * channel_ratio},
                  filename=os.path.join(outdir_dec, 'Dec{}.fits'.format(i_r)),
                  src_ra=mean_ra, src_dec=mean_dec, srcname=srcname)

        if HAS_MPL:
            plt.plot(dec_array, counts0)
            plt.plot(dec_array, counts1)

    if HAS_MPL:
        fig.savefig(os.path.join(outdir_dec, "allscans_dec.png"))
        plt.close(fig)


def main_simulate(args=None):
    """Preprocess the data."""
    import argparse

    description = ('Simulate a single scan or a map with a point source.')
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-s', "--source-flux", type=float, default=1,
                        help='Source flux in Jy')

    parser.add_argument('-n', "--noise-amplitude", type=float, default=1,
                        help='White noise amplitude')

    parser.add_argument('-b', "--baseline", type=str, default='flat',
                        help='Baseline kind: "flat", "slope" (linearly '
                             'increasing/decreasing), "messy" '
                             '(random walk) or a number (which gives an '
                             'amplitude to the random-walk baseline, that '
                             'would be 20 for "messy")')

    parser.add_argument('-g', '--geometry', nargs=4, type=float,
                        default=[120, 120, 120, 120],
                        help='Geometry specification: length_ra, length_dec, '
                             'width_ra, width_dec, in arcmins. A square map of'
                             ' 2 degrees would be specified as 120 120 120 '
                             '120. A cross-like map, 2x2 degrees wide but only'
                             ' along 1-degree stripes, is specified as 120 120'
                             ' 60 60')

    parser.add_argument('--beam-width', type=float, default=2.5,
                        help='Gaussian beam width in arcminutes')

    parser.add_argument('--spacing', type=float, default=0.5,
                        help='Spacing between scans in arcminutes '
                             '(default 0.5)')

    parser.add_argument('-o', "--outdir-root", type=str, default='sim',
                        help='Output directory root. Here, source and '
                             'calibrator scans/maps will be saved in '
                             'outdir/gauss_ra, outdir/gauss_dec, '
                             'outdir/calibrator1, outdir/calibrator2, where '
                             'outdir is the outdir root')

    parser.add_argument("--scan-speed", type=float, default=4.,
                        help='Scan speed in arcminutes/second')

    parser.add_argument("--integration-time", type=float, default=0.04,
                        help='Integration time in seconds')

    parser.add_argument("--no-cal", action='store_true', default=False,
                        help="Don't simulate calibrators")

    parser.add_argument("--debug", action='store_true', default=False,
                        help='Plot stuff and be verbose')

    args = parser.parse_args(args)

    def local_gauss_src_func(x, y):
        return args.source_flux * 100 * _2d_gauss(x, y,
                                                  sigma=args.beam_width/60)

    def calibrator_scan_func(x):
        return 100 * _2d_gauss(x, 0, sigma=args.beam_width/60)

    if not args.no_cal:
        cal1 = os.path.join(args.outdir_root, 'calibrator1')
        mkdir_p(cal1)
        sim_crossscans(5, cal1, scan_func=calibrator_scan_func,
                       channel_ratio=0.9, baseline=args.baseline)
        cal2 = os.path.join(args.outdir_root, 'calibrator2')
        mkdir_p(cal2)
        sim_crossscans(5, cal2, scan_func=calibrator_scan_func,
                       srcname='DummyCal2', channel_ratio=0.9,
                       baseline=args.baseline)

    simulate_map(dt=args.integration_time, length_ra=args.geometry[0],
                 length_dec=args.geometry[1], speed=args.scan_speed,
                 spacing=args.spacing, noise_amplitude=args.noise_amplitude,
                 width_ra=args.geometry[2], width_dec=args.geometry[3],
                 outdir=(os.path.join(args.outdir_root, 'gauss_ra'),
                         os.path.join(args.outdir_root, 'gauss_dec')),
                 baseline=args.baseline, mean_ra=180, mean_dec=70,
                 srcname='Dummy', channel_ratio=0.9,
                 count_map=local_gauss_src_func)
