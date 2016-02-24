"""Input/output functions."""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import astropy.io.fits as fits
from astropy.table import Table
import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation, AltAz
import os
from astropy.time import Time

DEBUG_MODE = False

locations = {'SRT': EarthLocation(4865182.7660, 791922.6890, 4035137.1740,
                                  unit=u.m),
             'Greenwich': EarthLocation(lat=51.477*u.deg, lon=0*u.deg)}


def _standard_offsets():
    # Case with fake multibeam data in files. Damn it
    radius = 0.0006671036
    # 0 for feed 0, radius for the other six
    radii = np.array([0] + [radius]*6)
    # Feeds 1--6 are at angles -60, -120, etc. Here I use angle 0 for
    # convenience for feed 0, but it has no effect since radii[0] is 0
    feed_angles = -np.arange(0, 7, 1) * np.pi * 2/6

    xoffsets = radii * np.cos(feed_angles)
    yoffsets = radii * np.sin(feed_angles)
    return xoffsets, yoffsets


def detect_data_kind(fname):
    """Placeholder for function that recognizes data format."""
    if fname.endswith('.hdf5'):
        return 'hdf5'
    else:
        return 'fitszilla'


def correct_offsets(derot_angle, xoffset, yoffset):
    """Correct feed offsets for derotation angle."""
    # Clockwise rotation of angle derot_angle
    new_xoff = xoffset * np.cos(derot_angle) - yoffset * np.sin(derot_angle)
    new_yoff = xoffset * np.sin(derot_angle) + yoffset * np.cos(derot_angle)

    return new_xoff, new_yoff


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
    global DEBUG_MODE

    # Open FITS file
    lchdulist = fits.open(fname)

    # ----------- Extract generic observation information ------------------
    source = lchdulist[0].header['SOURCE']
    receiver = lchdulist[0].header['HIERARCH RECEIVER CODE']
    ra = lchdulist[0].header['HIERARCH RIGHTASCENSION']
    dec = lchdulist[0].header['HIERARCH DECLINATION']
    # Check. If backend is not specified, use Total Power
    try:
        backend = lchdulist[0].header['HIERARCH BACKEND NAME']
    except:
        backend = 'TP'

    # ----------- Read the list of channel ids ------------------
    section_table_data = lchdulist['SECTION TABLE'].data
    chan_ids = section_table_data['id']

    # ----------- Read the list of RF inputs, feeds, polarization, etc. --
    rf_input_data = lchdulist['RF INPUTS'].data
    feeds = rf_input_data['feed']
    IFs = rf_input_data['ifChain']
    polarizations = rf_input_data['polarization']
    frequencies = rf_input_data['frequency']
    bandwidths = rf_input_data['bandWidth']

    # ----- Read the offsets of different feeds (nonzero only if multifeed)--
    feed_input_data = lchdulist['FEED TABLE'].data
    xoffsets = feed_input_data['xOffset']
    yoffsets = feed_input_data['yOffset']

    if DEBUG_MODE and len(xoffsets) > 1:
        xoffsets, yoffsets = _standard_offsets()

    relpowers = feed_input_data['relativePower']

    # -------------- Read data!-----------------------------------------
    datahdu = lchdulist['DATA TABLE']
    data_table_data = Table(datahdu.data)
    is_spectrum = 'SPECTRUM' in list(datahdu.header.values())
    if is_spectrum:
        nchan = len(chan_ids)

        nrows, nbins = data_table_data['SPECTRUM'].shape
        nbin_per_chan = nbins // nchan
        assert nbin_per_chan * nchan == nbins, \
            'Something wrong with channel subdivision'
        for ic, ch in enumerate(chan_ids):
            data_table_data['Ch{}'.format(ch)] = \
                data_table_data['SPECTRUM'][:, ic * nbin_per_chan:
                                            (ic + 1) * nbin_per_chan]

    info_to_retrieve = ['time', 'derot_angle']

    new_table = Table()

    new_table.meta['SOURCE'] = source
    new_table.meta['backend'] = backend
    new_table.meta['receiver'] = receiver
    new_table.meta['RA'] = ra
    new_table.meta['Dec'] = dec

    for info in info_to_retrieve:
        new_table[info] = data_table_data[info]

    if DEBUG_MODE:
        # Case with fake multibeam data in files. Still damn it
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

    # Coordinate correction. Will it work?
    for i in range(0, new_table['el'].shape[1]):
        # offsets < 0.001 arcseconds: don't correct (usually feed 0)
        if xoffsets[i] < np.radians(0.001 / 60.) and \
           yoffsets[i] < np.radians(0.001 / 60.):
            continue
        xoffs, yoffs = correct_offsets(new_table['derot_angle'],
                                       xoffsets[i],
                                       yoffsets[i])

        new_table['el'][:, i] += yoffs
        # TODO: Not sure about this cosine factor
        new_table['az'][:, i] += xoffs / np.cos(new_table['el'][:, i])

        obstimes = Time(new_table['time'] * u.day, format='mjd', scale='utc')
        coords = AltAz(az=new_table['az'][:, i],
                       alt=new_table['el'][:, i], unit=u.radian,
                       location=locations['SRT'],
                       obstime=obstimes)

        # According to line_profiler, coords.icrs is *by far* the longest
        # operation in this function, taking between 80 and 90% of the
        # execution time. Need to study a way to avoid this.
        coords_deg = coords.icrs
        new_table['ra'][:, i] = np.radians(coords_deg.ra)
        new_table['dec'][:, i] = np.radians(coords_deg.dec)

    for ic, ch in enumerate(chan_ids):
        if bandwidths[ic] < 0:
            frequencies[ic] -= bandwidths[ic]
            bandwidths[ic] *= -1
            for i in range(data_table_data['Ch{}'.format(ch)].shape[0]):
                data_table_data['Ch{}'.format(ch)][i, :] = \
                    data_table_data['Ch{}'.format(ch)][i, ::-1]
        new_table['Ch{}'.format(ch)] = \
            data_table_data['Ch{}'.format(ch)] * relpowers[feeds[ic]]

        new_table['Ch{}'.format(ch)].meta = {'polarization': polarizations[ic],
                                             'feed': feeds[ic],
                                             'IF': IFs[ic],
                                             'frequency': frequencies[ic],
                                             'bandwidth': bandwidths[ic],
                                             'xoffset': xoffsets[feeds[ic]],
                                             'yoffset': yoffsets[feeds[ic]],
                                             'relpower': relpowers[feeds[ic]]
                                             }
        new_table['Ch{}_feed'.format(ch)] = \
            np.zeros(len(data_table_data), dtype=np.uint8) + feeds[ic]

        new_table['Ch{}-filt'.format(ch)] = \
            np.ones(len(data_table_data['Ch{}'.format(ch)]), dtype=bool)
    lchdulist.close()
    return new_table


def read_data(fname):
    """Read the data, whatever the format, and return them."""
    kind = detect_data_kind(fname)
    if kind == 'fitszilla':
        return read_data_fitszilla(fname)
    elif kind == 'hdf5':
        return Table.read(fname)


def root_name(fname):
    """Return the file name without extension."""
    import os
    return os.path.splitext(fname)[0]


def profile_coords():
    """Same test above, with profiling."""
    global DEBUG_MODE
    DEBUG_MODE = True
    curdir = os.path.abspath(os.path.dirname(__file__))
    datadir = os.path.join(curdir, '..', '..', 'TEST_DATASET')

    fname = os.path.join(datadir, '20150410-001307-scicom-W44',
                         '20150410-001307-scicom-W44_002_003.fits')
    new_table = read_data(fname)
    return new_table
