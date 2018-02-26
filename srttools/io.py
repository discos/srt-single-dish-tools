"""Input/output functions."""
from __future__ import (absolute_import, division,
                        print_function)
import astropy.io.fits as fits
from astropy.table import Table
import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation, AltAz, Angle, ICRS
import os
from astropy.time import Time
import logging
import warnings


__all__ = ["mkdir_p", "detect_data_kind", "correct_offsets", "observing_angle",
           "get_rest_angle", "print_obs_info_fitszilla", "read_data_fitszilla",
           "read_data", "root_name"]


locations = {'srt': EarthLocation(4865182.7660, 791922.6890, 4035137.1740,
                                  unit=u.m),
             'medicina': EarthLocation(Angle("11:38:49", u.deg),
                                       Angle("44:31:15", u.deg),
                                       25 * u.meter),
             'greenwich': EarthLocation(lat=51.477*u.deg, lon=0*u.deg)}


def mkdir_p(path):
    """Safe mkdir function.

    Parameters
    ----------
    path : str
        Name of the directory/ies to create

    Notes
    -----
    Found at
    http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
    """
    import errno
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def _check_derotator(derot_angle):
    # Check that derotator angle is outside any plausible value
    if np.any(np.abs(derot_angle) > 2*360):
        return False
    return True


def detect_data_kind(fname):
    """Placeholder for function that recognizes data format."""
    if fname.endswith('.hdf5'):
        return 'hdf5'
    else:
        return 'fitszilla'


def correct_offsets(obs_angle, xoffset, yoffset):
    """Correct feed offsets for derotation angle.

    All angles are in radians.

    Examples
    --------
    >>> x = 2 ** 0.5
    >>> y = 2 ** 0.5
    >>> angle = np.pi / 4
    >>> xoff, yoff = correct_offsets(angle, x, y)
    >>> np.allclose([xoff, yoff], 2 ** 0.5)
    True

    """
    sep = np.sqrt(xoffset**2. + yoffset**2.)

    new_xoff = sep * np.cos(obs_angle)
    new_yoff = sep * np.sin(obs_angle)

    return new_xoff, new_yoff


def observing_angle(rest_angle, derot_angle):
    """Calculate the observing angle of the multifeed.

    If values have no units, they are assumed in radians

    Parameters
    ----------
    rest_angle : float or Astropy quantity, angle
        rest angle of the feeds
    derot_angle : float or Astropy quantity, angle
        derotator angle

    Examples
    --------
    >>> observing_angle(0 * u.rad, 2 * np.pi * u.rad).to(u.rad).value
    0.0
    >>> observing_angle(0, 2 * np.pi).to(u.rad).value
    0.0
    """
    if not hasattr(rest_angle, 'unit'):
        rest_angle *= u.rad
    if not hasattr(derot_angle, 'unit'):
        derot_angle *= u.rad
    return rest_angle + (2 * np.pi * u.rad - derot_angle)


def _rest_angle_default(n_lat_feeds):
    """Default rest angles for a multifeed, in units of a circle

    Assumes uniform coverage.

    Examples
    --------
    >>> np.allclose(_rest_angle_default(5),
    ...             np.array([1., 0.8, 0.6, 0.4, 0.2]))
    True
    >>> np.allclose(_rest_angle_default(6) * 360,
    ...             np.array([360., 300., 240., 180., 120., 60.]))
    True
    """
    return np.arange(1, 0, -1 / n_lat_feeds)


def get_rest_angle(xoffsets, yoffsets):
    """Calculate the rest angle for multifeed.

    The first feed is assumed to be at position 0, for it the return value is 0

    Examples
    --------
    >>> xoffsets = [0.0, -0.0382222, -0.0191226, 0.0191226, 0.0382222,
    ...             0.0191226, -0.0191226]
    >>> yoffsets = [0.0, 0.0, 0.0331014, 0.0331014, 0.0, -0.0331014,
    ...             -0.0331014]
    >>> np.allclose(get_rest_angle(xoffsets, yoffsets).to(u.deg).value,
    ...             np.array([0., 180., 120., 60., 360., 300., 240.]))
    True
    """
    if len(xoffsets) <= 2:
        return np.array([0]*len(xoffsets))
    xoffsets = np.asarray(xoffsets)
    yoffsets = np.asarray(yoffsets)
    n_lat_feeds = len(xoffsets) - 1
    rest_angle_default = _rest_angle_default(n_lat_feeds) * 2 * np.pi * u.rad
    w_0 = np.where((xoffsets[1:] > 0) & (yoffsets[1:] == 0.))[0][0]
    return np.concatenate(([0],
                           np.roll(rest_angle_default.to(u.rad).value,
                                   w_0))) * u.rad


def print_obs_info_fitszilla(fname):
    """Placeholder for function that prints out oberving information."""
    lchdulist = fits.open(fname)
    section_table_data = lchdulist['SECTION TABLE'].data
    sample_rates = section_table_data['sampleRate']

    print('Sample rates:', sample_rates)

    rf_input_data = lchdulist['RF INPUTS'].data
    print('Feeds          :', rf_input_data['feed'])
    print('IFs            :', rf_input_data['ifChain'])
    print('Polarizations  :', rf_input_data['polarization'])
    print('Frequencies    :', rf_input_data['frequency'])
    print('Bandwidths     :', rf_input_data['bandWidth'])

    lchdulist.close()


def read_data_fitszilla(fname):
    """Open a fitszilla FITS file and read all relevant information."""

    # Open FITS file
    lchdulist = fits.open(fname)

    # ----------- Extract generic observation information ------------------
    source = lchdulist[0].header['SOURCE']
    site = lchdulist[0].header['ANTENNA'].lower()
    receiver = lchdulist[0].header['HIERARCH RECEIVER CODE']
    ra = lchdulist[0].header['HIERARCH RIGHTASCENSION'] * u.rad
    dec = lchdulist[0].header['HIERARCH DECLINATION'] * u.rad
    ra_offset = dec_offset = az_offset = el_offset = 0 * u.rad
    if 'HIERARCH RightAscension Offset' in lchdulist[0].header:
        ra_offset = \
            lchdulist[0].header['HIERARCH RightAscension Offset'] * u.rad
    if 'HIERARCH Declination Offset' in lchdulist[0].header:
        dec_offset = lchdulist[0].header['HIERARCH Declination Offset'] * u.rad
    if 'HIERARCH Azimuth Offset' in lchdulist[0].header:
        az_offset = lchdulist[0].header['HIERARCH Azimuth Offset'] * u.rad
    if 'HIERARCH Elevation Offset' in lchdulist[0].header:
        el_offset = lchdulist[0].header['HIERARCH Elevation Offset'] * u.rad

    # Check. If backend is not specified, use Total Power
    try:
        backend = lchdulist[0].header['HIERARCH BACKEND NAME']
    except Exception:
        backend = 'TP'

    # ----------- Read the list of channel ids ------------------
    section_table_data = lchdulist['SECTION TABLE'].data
    chan_ids = section_table_data['id']
    nbin_per_chan = section_table_data['bins']
    if len(list(set(nbin_per_chan))) > 1:
        lchdulist.close()

        raise ValueError('Only datasets with the same nbin per channel are '
                         'supported at the moment')
    nbin_per_chan = list(set(nbin_per_chan))[0]
    types = section_table_data['type']
    if 'stokes' in types:
        is_polarized = True
    else:
        is_polarized = False

    # ----------- Read the list of RF inputs, feeds, polarization, etc. --
    rf_input_data = lchdulist['RF INPUTS'].data
    feeds = rf_input_data['feed']
    IFs = rf_input_data['ifChain']
    polarizations = rf_input_data['polarization']
    chan_names = ['Feed{}_{}'.format(f, p)
                  for f, p in zip(feeds, polarizations)]
    frequencies = rf_input_data['frequency']
    bandwidths = rf_input_data['bandWidth']

    # ----- Read the offsets of different feeds (nonzero only if multifeed)--
    feed_input_data = lchdulist['FEED TABLE'].data
    # Add management of historical offsets
    xoffsets = feed_input_data['xOffset'] * u.rad
    yoffsets = feed_input_data['yOffset'] * u.rad

    relpowers = feed_input_data['relativePower']

    # -------------- Read data!-----------------------------------------
    datahdu = lchdulist['DATA TABLE']
    data_table_data = Table(datahdu.data)
    for col in data_table_data.colnames:
        if col == col.lower():
            continue
        data_table_data.rename_column(col, col.lower())

    is_spectrum = 'SPECTRUM' in list(datahdu.header.values())
    if is_spectrum:
        nchan = len(chan_ids)

        _, nbins = data_table_data['spectrum'].shape

        if nbin_per_chan * nchan * 2 == nbins and not is_polarized:
            warnings.warn('Data appear to contain polarization information '
                          'but are classified as simple, not stokes, in the '
                          'Section table.')
            is_polarized = True
        if nbin_per_chan * nchan != nbins and \
                nbin_per_chan * nchan * 2 != nbins:
            lchdulist.close()

            raise ValueError('Something wrong with channel subdivision: '
                             '{} bins/channel, {} channels, '
                             '{} total bins'.format(nbin_per_chan, nchan,
                                                    nbins))

        for ic, ch in enumerate(chan_names):
            data_table_data[ch] = \
                data_table_data['spectrum'][:, ic * nbin_per_chan:
                                            (ic + 1) * nbin_per_chan]
        if is_polarized:
            if len(list(set(feeds))) > 1:
                raise ValueError('Polarized data are only supported for single'
                                 ' feed observations')
            feed = feeds[0]
            data_table_data['Feed{}_Q'.format(feed)] = \
                data_table_data['spectrum'][:, 2 * nbin_per_chan:
                                               3 * nbin_per_chan]
            data_table_data['Feed{}_U'.format(feed)] = \
                data_table_data['spectrum'][:, 3 * nbin_per_chan:
                                               4 * nbin_per_chan]
            chan_names += ['Feed{}_Q'.format(feed),
                           'Feed{}_U'.format(feed)]
    else:
        for ic, ch in enumerate(chan_names):
            data_table_data[ch] = \
                data_table_data['ch{}'.format(chan_ids[ic])]

    # ----------- Read temperature data, if possible ----------------
    tempdata = lchdulist['ANTENNA TEMP TABLE'].data
    try:
        for ic, ch in enumerate(chan_names):
            td = tempdata['ch{}'.format(chan_ids[ic])]
            data_table_data[ch + '-Temp'] = td
    except Exception as e:
        logging.warning("Could not read temperature information from file."
                        "Exception: {}".format(str(e)))
        for ic, ch in enumerate(chan_names):
            data_table_data[ch + '-Temp'] = 0.

    info_to_retrieve = \
        ['time', 'derot_angle'] + [ch + '-Temp' for ch in chan_names]

    new_table = Table()

    new_table.meta['SOURCE'] = source
    new_table.meta['site'] = site
    new_table.meta['backend'] = backend
    new_table.meta['receiver'] = receiver
    new_table.meta['RA'] = ra
    new_table.meta['Dec'] = dec
    for i, off in zip("ra,dec,el,az".split(','),
                      [ra_offset, dec_offset, el_offset, az_offset]):
        new_table.meta[i + "_offset"] = off

    for info in info_to_retrieve:
        new_table[info] = data_table_data[info]

    if not _check_derotator(new_table['derot_angle']):
        logging.warning('Derotator angle looks weird. Setting to 0')
        new_table['derot_angle'][:] = 0

    # Duplicate raj and decj columns (in order to be corrected later)
    new_table['ra'] = \
        np.tile(data_table_data['raj2000'],
                (np.max(feeds) + 1, 1)).transpose()
    new_table['dec'] = \
        np.tile(data_table_data['decj2000'],
                (np.max(feeds) + 1, 1)).transpose()
    new_table['el'] = \
        np.tile(data_table_data['el'],
                (np.max(feeds) + 1, 1)).transpose()
    new_table['az'] = \
        np.tile(data_table_data['az'],
                (np.max(feeds) + 1, 1)).transpose()

    for info in ['ra', 'dec', 'az', 'el', 'derot_angle']:
        new_table[info].unit = u.radian

    rest_angles = get_rest_angle(xoffsets, yoffsets)

    for i in range(0, new_table['el'].shape[1]):
        # offsets < 0.001 arcseconds: don't correct (usually feed 0)
        if np.abs(xoffsets[i]) < np.radians(0.001 / 60.) * u.rad and \
           np.abs(yoffsets[i]) < np.radians(0.001 / 60.) * u.rad:
            continue

        # Calculate observing angle
        obs_angle = observing_angle(rest_angles[i], new_table['derot_angle'])

        xoffs, yoffs = correct_offsets(obs_angle, xoffsets[i], yoffsets[i])

        new_table['el'][:, i] += yoffs.to(u.rad).value
        new_table['az'][:, i] += \
            xoffs.to(u.rad).value / np.cos(new_table['el'][:, i])

        obstimes = Time(new_table['time'] * u.day, format='mjd', scale='utc')

        coords = AltAz(az=Angle(new_table['az'][:, i]),
                       alt=Angle(new_table['el'][:, i]),
                       location=locations[site],
                       obstime=obstimes)

        # According to line_profiler, coords.icrs is *by far* the longest
        # operation in this function, taking between 80 and 90% of the
        # execution time. Need to study a way to avoid this.
        coords_deg = coords.transform_to(ICRS)
        new_table['ra'][:, i] = np.radians(coords_deg.ra)
        new_table['dec'][:, i] = np.radians(coords_deg.dec)

    for ic, ch in enumerate(chan_ids):
        chan_name = chan_names[ic]
        if bandwidths[ic] < 0:
            frequencies[ic] -= bandwidths[ic]
            bandwidths[ic] *= -1
            for i in range(
                    data_table_data[chan_name].shape[0]):
                data_table_data[chan_name][i, :] = \
                    data_table_data[chan_name][i, ::-1]

        new_table[chan_name] = \
            data_table_data[chan_name] * relpowers[feeds[ic]]

        newmeta = \
            {'polarization': polarizations[ic],
             'feed': int(feeds[ic]),
             'IF': int(IFs[ic]),
             'frequency': float(frequencies[ic]),
             'bandwidth': float(bandwidths[ic]),
             'xoffset': float(xoffsets[feeds[ic]].to(u.rad).value) * u.rad,
             'yoffset': float(yoffsets[feeds[ic]].to(u.rad).value) * u.rad,
             'relpower': float(relpowers[feeds[ic]])
             }
        new_table[chan_name].meta.update(newmeta)

        new_table[chan_name + '-filt'] = \
            np.ones(len(data_table_data[chan_name]), dtype=bool)

    if is_polarized:
        for feed in list(set(feeds)):
            for stokes_par in 'QU':
                chan_name = 'Feed{}_{}'.format(feed, stokes_par)
                new_table[chan_name] = \
                    data_table_data[chan_name]

                newmeta = \
                    {'polarization': stokes_par,
                     'feed': int(feed),
                     'IF': -1,
                     'frequency': float(frequencies[feed * 0]),
                     'bandwidth': float(bandwidths[feed * 0]),
                     'xoffset': float(
                         xoffsets[feed].to(u.rad).value) * u.rad,
                     'yoffset': float(
                         yoffsets[feed].to(u.rad).value) * u.rad,
                     'relpower': 1.
                     }
                new_table[chan_name].meta.update(newmeta)

                new_table[chan_name + '-filt'] = \
                    np.ones(len(data_table_data[chan_name]), dtype=bool)

    lchdulist.close()
    return new_table


def read_data(fname):
    """Read the data, whatever the format, and return them."""
    kind = detect_data_kind(fname)
    if kind == 'fitszilla':
        return read_data_fitszilla(fname)
    elif kind == 'hdf5':
        return Table.read(fname, path='scan')


def root_name(fname):
    """Return the file name without extension."""
    return os.path.splitext(fname)[0]
