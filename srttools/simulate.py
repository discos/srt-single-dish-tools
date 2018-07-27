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
from .utils import tqdm, jit
from astropy.coordinates import SkyCoord
from astropy.time import Time
try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

__all__ = ["simulate_scan", "save_scan", "simulate_map"]


DEFAULT_CAL_OFFSET = 211


summary_header = """
SIMPLE  =                    T / file does conform to FITS standard
BITPIX  =                    8 / number of bits per data pixel
NAXIS   =                    0 / number of data axes
EXTEND  =                    T / FITS dataset may contain extensions
COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy
COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H
HIERARCH BackendName = 'NULL    ' / Backend name
CREATOR = 'NULL    '           / Software (incl. version)
DATE-OBS= '2016-10-03T14:59:08.753' / Observation time
HIERARCH Declination = 0.253290907695677 / Target declination (radians)
EQUINOX =                   0. / Equinox of RA, Dec
EXPTIME =                   0. / Total integration time (seconds)
FITSVER = 'V.1.11  '           / FITS version
LST     =                    0 / Local sidereal time
HIERARCH LogFileName = 'NULL    ' / Name of the log file
HIERARCH NUSEBANDS =         0 / Number of sections
OBJECT  = 'W51     '           / Target source name
OBSID   = 'NULL    '           / Observer or operator initials
PROJID  = 'NULL    '           / ProjectID
HIERARCH RESTFREQ1 =  22235.08 / Rest frequency (MHz)
HIERARCH RESTFREQ2 =  22235.08 / Rest frequency (MHz)
HIERARCH RESTFREQ3 =  22235.08 / Rest frequency (MHz)
HIERARCH RESTFREQ4 =  22235.08 / Rest frequency (MHz)
HIERARCH ReceiverCode = 'CCB     ' / Receiver name
HIERARCH RightAscension = 5.07757730974885 / Target right ascension (radians)
SCANGEOM= 'NULL    '           / Scan geometry
SCANMODE= 'NULL    '           / Mapping mode
SCANTYPE= 'NULL    '           / Scan astronomical type
SCANXVEL=                   0. / Tracking rate (optional,OTF)
SWTCHMOD= 'NULL    '           / Switch mode
HIERARCH ScheduleName = 'NULL    ' / Name of the schedule
TELESCOP= 'SRT     '           / Telescope name
VDEF    = 'OP      '           / Radial velocity definition
VFRAME  = 'LSRK    '           / Radial velocity reference frame
VRAD    =                    0 / Radial velocity
WOBUSED =                    0 / Wobbler used?

"""


def _apply_spectrum_to_data(spec_func, counts, nbin, bw=1000):
    if nbin == 1:
        return counts
    single = False
    if not isinstance(counts, collections.Iterable):
        counts = [counts]
        single = True
    counts = np.asarray(counts)
    df = bw / nbin
    freqs = np.arange(0, bw, df)
    single_spec = spec_func(freqs)
    spec = np.zeros((len(counts), len(freqs)))
    for i, c in enumerate(counts):
        spec[i, :] += c * single_spec

    if single:
        return spec[0]
    return spec


def _standard_source_spectrum(counts, nbin, bw=1000, sigma=1):
    def spec_func(f):
        f = f - bw / 2
        return np.exp(-(f ** 2) / (2 * sigma ** 2))
    return _apply_spectrum_to_data(spec_func, counts, nbin, bw)


def _standard_bkg_spectrum(counts, nbin, bw=1000):
    def spec_func(f):
        sp = 1 + 0.1 * np.sin(2 * np.pi * 5 / bw * f) * (1 - f / bw)
        sp -= 0.5 * f / bw
        return sp
    return _apply_spectrum_to_data(spec_func, counts, nbin, bw)


def create_summary(filename, key_dict=None):
    if key_dict is None:
        key_dict = {}
    header = fits.Header.fromstring(summary_header, sep='\n')
    for key, value in key_dict.items():
        header[key] = value

    primary_hdu = fits.PrimaryHDU(header=header)
    hdul = fits.HDUList([primary_hdu])
    hdul.writeto(filename, overwrite=True)
    return filename


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
    except (ValueError, TypeError):
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


@jit(nopython=True)
def _2d_gauss(x, y, sigma=2.5 / 60.):
    """A Gaussian beam"""
    return np.exp(-(x ** 2 + y ** 2) / (2 * sigma**2))


@jit(nopython=True)
def calibrator_scan_func(x):
    return 100 * _2d_gauss(x, 0, sigma=2.5 / 60)


def sim_crossscans(ncross, caldir, scan_func=calibrator_scan_func,
                   srcname='DummyCal', channel_ratio=0.8, baseline="flat",
                   nbin=1):
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
                          baseline=baseline, nbin=nbin)
        _, _, scan1 = \
            simulate_scan(dt=dt, length=length, speed=speed, shape=scan_func,
                          noise_amplitude=0.2, center=0,
                          baseline=baseline, nbin=nbin)

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
                          baseline=baseline, nbin=nbin)
        _, _, scan1 = \
            simulate_scan(dt=dt, length=length, speed=speed, shape=scan_func,
                          noise_amplitude=0.2, center=src_dec,
                          baseline=baseline, nbin=nbin)

        if i % 2 != 0:
            decs = decs[::-1]

        ras = np.zeros_like(decs) + src_ra

        save_scan(times + timedelta, ras, decs,
                  {'Ch0': scan0, 'Ch1': scan1 * channel_ratio},
                  filename=os.path.join(caldir, '{}_Dec.fits'.format(i)),
                  src_ra=src_ra, src_dec=src_dec, srcname=srcname,
                  counts_to_K=(0.03, 0.03 / channel_ratio))
        timedelta += times[-1] + 1

    create_summary(os.path.join(caldir, 'summary.fits'),
                   {'RightAscension': np.radians(src_ra),
                    'Declination': np.radians(src_dec),
                    'Object': srcname})


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


def sim_position_switching(caldir, srcname='Dummy', nbin=1,
                           offset=np.radians(3), strategy=None):
    dt = 0.04
    src_ra = 185
    src_dec = 75
    if strategy is None:
        strategy = [1, 1, 1]

    last_time = 0
    for n_on in range(strategy[0]):
        times, ras, on = \
            simulate_scan(dt=dt, length=0, speed=1, baseline=(0, 10, 0),
                          noise_amplitude=0.2, center=src_dec, nbin=nbin)

        times += last_time
        last_time = times[0]
        decs = np.zeros_like(ras) + src_dec
        save_scan(times, ras, decs,
                  {'Ch0': on, 'Ch1': on},
                  filename=os.path.join(caldir, 'ON_{}.fits'.format(n_on)),
                  src_ra=src_ra, src_dec=src_dec, srcname=srcname,
                  counts_to_K=(0.03, 0.03),
                  other_keywords={'SIGNAL': 'SIGNAL'})


    for n_off in range(strategy[1]):
        times, _, off = \
            simulate_scan(dt=dt, length=0, speed=1, baseline=(0, 10, 0),
                          shape=lambda x: 0,
                          noise_amplitude=0.2, center=src_dec, nbin=nbin)

        times += last_time
        last_time = times[0]
        save_scan(times, ras + offset, decs, {'Ch0': off, 'Ch1': off},
                  filename=os.path.join(caldir, 'OFF_{}.fits'.format(n_off)),
                  src_ra=src_ra, src_dec=src_dec, srcname=srcname,
                  counts_to_K=(0.03, 0.03),
                  other_keywords={'RightAscension Offset': offset,
                                  'SIGNAL': 'REFERENCE'})

    for n_cal in range(strategy[2]):
        times, _, cal = \
            simulate_scan(dt=dt, length=0, speed=1, baseline=(0, 10, 0),
                          shape=lambda x: 0,
                          noise_amplitude=0.2, center=src_dec, nbin=nbin,
                          calon=True)

        times += last_time
        last_time = times[0]
        save_scan(times, ras + offset, decs, {'Ch0': cal, 'Ch1': cal},
                  filename=os.path.join(caldir, 'CAL_{}.fits'.format(n_cal)),
                  src_ra=src_ra, src_dec=src_dec, srcname=srcname,
                  counts_to_K=(0.03, 0.03),
                  other_keywords={'RightAscension Offset': offset,
                                  'SIGNAL': 'REFERENCE'},
                  other_columns={'flag_cal': 1})

    create_summary(os.path.join(caldir, 'summary.fits'),
                   {'RightAscension': np.radians(src_ra),
                    'Declination': np.radians(src_ra),
                    'Object': srcname})
    return caldir


def simulate_scan(dt=0.04, length=120., speed=4., shape=None,
                  noise_amplitude=1., center=0., baseline="flat",
                  nbin=1, calon=False, nsamples=None):
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
    baseline : str, number or tuple
        "flat", "slope" (linearly increasing/decreasing), "messy"
        (random walk), a number (which gives an amplitude to the random-walk
        baseline, that is 20 for "messy"), or a tuple (m, q, messy_amp) giving
        the maximum and minimum absolute-value slope and intercept, and the
        random-walk amplitude.
    """
    if shape is None:
        shape = _default_flat_shape

    if nsamples is None and length == 0:
        nsamples = 100
    elif nsamples is None:
        nsamples = np.rint(length / speed / dt)

    times = np.arange(nsamples) * dt
    # In degrees!
    position = np.arange(-nsamples / 2, nsamples / 2) / nsamples * length / 60

    scan_baseline = _create_baseline(position, baseline)

    signal = _standard_source_spectrum(shape(position), nbin)
    bkg = _standard_bkg_spectrum(scan_baseline, nbin)
    scan_shape = signal + bkg
    scan_shape += ra.normal(0, noise_amplitude, scan_shape.shape)

    if calon:
        scan_shape += DEFAULT_CAL_OFFSET

    return times, position + center, scan_shape


def save_scan(times, ra, dec, channels, filename='out.fits',
              other_columns=None, other_keywords=None, scan_type=None,
              src_ra=None, src_dec=None, srcname='Dummy', counts_to_K=0.03):
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
    if other_columns is None:
        other_columns = {}
    if other_keywords is None:
        other_keywords = {}
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
    secthdu = lchdulist['SECTION TABLE']

    lchdulist[0].header['SOURCE'] = "Dummy"
    lchdulist[0].header['ANTENNA'] = "SRT"
    lchdulist[0].header['HIERARCH RIGHTASCENSION'] = np.radians(src_ra)
    lchdulist[0].header['HIERARCH DECLINATION'] = np.radians(src_dec)
    if scan_type is not None:
        lchdulist[0].header['HIERARCH SubScanType'] = scan_type

    for key in other_keywords.keys():
        lchdulist[0].header[key] = other_keywords[key]

    data_table_data = Table(datahdu.data)
    data_table_data.remove_column('Ch0')
    data_table_data.remove_column('Ch1')

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

    for col in other_columns.keys():
        newtable[col] = other_columns[col]

    data_table_data = vstack([data_table_data, newtable])

    hdu = fits.BinTableHDU(data_table_data, header=datahdu.header)
    nrows = len(data_table_data)

    datahdu.data = hdu.data

    temptable = Table()
    for ch in channels.keys():
        dummy_data = newtable[ch]
        if len(dummy_data.shape) == 2:
            dummy_data = np.sum(dummy_data, axis=1)
        temptable[ch] = dummy_data * counts_to_K[ch]

    thdu = fits.BinTableHDU.from_columns(temphdu.data.columns, nrows=nrows)
    for colname in temphdu.data.columns.names:
        thdu.data[colname][:] = temptable[colname]

    temphdu.data = thdu.data

    shape = channels['Ch0'].shape
    if len(shape) == 2:
        secthdu.data['bins'] = shape[1]

    lchdulist[0].header['SOURCE'] = srcname
    lchdulist.writeto(filename, overwrite=True)
    lchdulist.close()


def _single_value_as_tuple(value, nvals=2):
    """If a value is single, return as a tuple.

    Examples
    --------
    >>> np.all(_single_value_as_tuple(1) == (1, 1))
    True
    >>> np.all(_single_value_as_tuple((1, 1, 1)) == (1, 1, 1))
    True
    >>> np.all(_single_value_as_tuple(1, nvals=3) == (1, 1, 1))
    True
    """
    if isinstance(value, collections.Iterable):
        return value
    return tuple([value] * nvals)


def _create_baseline(x, baseline_kind="flat"):
    """
    Parameters
    ----------
    x : float, array-like
        The x values for the baseline
    baseline : str, number of tuple
        "flat", "slope" (linearly increasing/decreasing), "messy"
        (random walk), a number (which gives an amplitude to the random-walk
        baseline, that is 20 for "messy"), or a tuple (m, q, messy_amp) giving
        the maximum and minimum absolute-value slope and intercept, and the
        random-walk amplitude.
    """
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
    elif isinstance(baseline_kind, collections.Iterable) and not \
            isinstance(baseline_kind, six.string_types):
        m = _single_value_as_tuple(baseline_kind[0], nvals=2)
        q = _single_value_as_tuple(baseline_kind[1], nvals=2)
        mmin, mmax = m[0], m[1]
        qmin, qmax = q[0], q[1]
        stochastic_amp = float(baseline_kind[2])
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
                 srcname='Dummy', channel_ratio=1, nbin=1):

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

        counts_clean = \
            _standard_source_spectrum(count_map(ra_array, delta_dec),
                                      nbin=nbin)

        baseline0 = \
            _standard_bkg_spectrum(_create_baseline(ra_array, baseline),
                                   nbin=nbin)
        baseline1 = \
            _standard_bkg_spectrum(_create_baseline(ra_array, baseline),
                                   nbin=nbin)

        counts0 = counts_clean + \
            ra.normal(0, noise_amplitude, counts_clean.shape) + baseline0
        counts1 = counts_clean + \
            ra.normal(0, noise_amplitude, counts_clean.shape) + baseline1

        actual_ra = mean_ra + ra_array / np.cos(np.radians(start_dec))

        if i_d % 2 != 0:
            actual_ra = actual_ra[::-1]
        fname = os.path.join(outdir_ra, 'Ra{}.fits'.format(i_d))
        other_keywords = {'Declination Offset': delta_dec}
        save_scan(times_ra, actual_ra, np.zeros_like(actual_ra) + start_dec,
                  {'Ch0': counts0, 'Ch1': counts1 * channel_ratio},
                  filename=fname, other_keywords=other_keywords,
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

        counts_clean = \
            _standard_source_spectrum(count_map(delta_ra, dec_array),
                                      nbin=nbin)

        baseline0 = \
            _standard_bkg_spectrum(_create_baseline(dec_array, baseline),
                                   nbin=nbin)
        baseline1 = \
            _standard_bkg_spectrum(_create_baseline(dec_array, baseline),
                                   nbin=nbin)

        counts0 = counts_clean + \
            ra.normal(0, noise_amplitude, counts_clean.shape) + baseline0
        counts1 = counts_clean + \
            ra.normal(0, noise_amplitude, counts_clean.shape) + baseline1

        if i_r % 2 != 0:
            dec_array = dec_array[::-1]
        other_keywords = {'RightAscension Offset': delta_ra}
        save_scan(times_dec, np.zeros_like(dec_array) + start_ra,
                  dec_array + mean_dec,
                  {'Ch0': counts0, 'Ch1': counts1 * channel_ratio},
                  other_keywords=other_keywords,
                  filename=os.path.join(outdir_dec, 'Dec{}.fits'.format(i_r)),
                  src_ra=mean_ra, src_dec=mean_dec, srcname=srcname)

        if HAS_MPL:
            plt.plot(dec_array, counts0)
            plt.plot(dec_array, counts1)

    if HAS_MPL:
        fig.savefig(os.path.join(outdir_dec, "allscans_dec.png"))
        plt.close(fig)

    print("Creating summary...")
    create_summary(os.path.join(outdir_ra, 'summary.fits'),
                   {'RightAscension': np.radians(mean_ra),
                    'Declination': np.radians(mean_dec),
                    'Object': srcname})
    if outdir_ra == outdir_dec:
        return outdir_ra, outdir_ra
    create_summary(os.path.join(outdir_dec, 'summary.fits'),
                   {'RightAscension': np.radians(mean_ra),
                    'Declination': np.radians(mean_dec),
                    'Object': srcname})
    return outdir_ra, outdir_dec


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

    parser.add_argument("--spectral-bins", type=int, default=1,
                        help='Simulate a spectrum with this number of bins')

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
                       channel_ratio=0.9, baseline=args.baseline,
                       nbin=args.spectral_bins)
        cal2 = os.path.join(args.outdir_root, 'calibrator2')
        mkdir_p(cal2)
        sim_crossscans(5, cal2, scan_func=calibrator_scan_func,
                       srcname='DummyCal2', channel_ratio=0.9,
                       baseline=args.baseline,
                       nbin=args.spectral_bins)

    simulate_map(dt=args.integration_time, length_ra=args.geometry[0],
                 length_dec=args.geometry[1], speed=args.scan_speed,
                 spacing=args.spacing, noise_amplitude=args.noise_amplitude,
                 width_ra=args.geometry[2], width_dec=args.geometry[3],
                 outdir=(os.path.join(args.outdir_root, 'gauss_ra'),
                         os.path.join(args.outdir_root, 'gauss_dec')),
                 baseline=args.baseline, mean_ra=180, mean_dec=70,
                 srcname='Dummy', channel_ratio=0.9,
                 count_map=local_gauss_src_func,
                 nbin=args.spectral_bins)
