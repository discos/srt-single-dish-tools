from __future__ import print_function, division
from astropy.io import fits
import astropy.units as u
from astropy.time import Time
from astropy.table import Table
import numpy as np
import copy
import warnings

from .io import get_coords_from_altaz_offset, correct_offsets
from .io import get_rest_angle, observing_angle, locations


def convert_to_complete_fitszilla(fname, outname):
    assert outname != fname, 'Files cannot have the same name'

    lchdulist = fits.open(fname)

    feed_input_data = lchdulist['FEED TABLE'].data
    xoffsets = feed_input_data['xOffset'] * u.rad
    yoffsets = feed_input_data['yOffset'] * u.rad
    # ----------- Extract generic observation information ------------------
    site = lchdulist[0].header['ANTENNA'].lower()
    location = locations[site]

    rest_angles = get_rest_angle(xoffsets, yoffsets)

    datahdu = lchdulist['DATA TABLE']
    data_table_data = Table(datahdu.data)

    new_table = Table()
    info_to_retrieve = \
        ['time', 'derot_angle', 'el', 'az', 'raj2000', 'decj2000']
    for info in info_to_retrieve:
        new_table[info.replace('j2000', '')] = data_table_data[info]

    el_save = new_table['el']
    az_save = new_table['az']
    derot_angle = new_table['derot_angle']
    el_save.unit = u.rad
    az_save.unit = u.rad
    derot_angle.unit = u.rad
    times = new_table['time']

    for i in range(len(xoffsets)):
        obs_angle = observing_angle(rest_angles[i], derot_angle)

        # offsets < 0.001 arcseconds: don't correct (usually feed 0)
        if np.abs(xoffsets[i]) < np.radians(0.001 / 60.) * u.rad and \
                np.abs(yoffsets[i]) < np.radians(0.001 / 60.) * u.rad:
            continue
        el = copy.deepcopy(el_save)
        az = copy.deepcopy(az_save)
        xoffs, yoffs = correct_offsets(obs_angle, xoffsets[i], yoffsets[i])
        obstimes = Time(times * u.day, format='mjd', scale='utc')

        # el and az are also changed inside this function (inplace is True)
        ra, dec = \
            get_coords_from_altaz_offset(obstimes, el, az, xoffs, yoffs,
                                         location=location, inplace=True)
        ra = fits.Column(array=ra, name='raj2000', format='1D')
        dec = fits.Column(array=dec, name='decj2000', format='1D')
        el = fits.Column(array=el, name='el', format='1D')
        az = fits.Column(array=az, name='az', format='1D')
        new_data_extension = \
            fits.BinTableHDU.from_columns([ra, dec, el, az])
        new_data_extension.name = 'Coord{}'.format(i)
        lchdulist.append(new_data_extension)

    lchdulist.writeto(outname + '.fits', overwrite=True)


def converter_main(args=None):
    import argparse

    description = ('Load a series of scans and convert them to various'
                   'formats')
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("files", nargs='*',
                        help="Single files to process",
                        default=None, type=str)

    parser.add_argument("-f", "--format", type=str, default='fitsmod',
                        help='Format of output files (default '
                             'fitszilla_mod, indicating a fitszilla with '
                             'converted coordinates for feed number *n* in '
                             'a separate COORDn extensions)')

    args = parser.parse_args(args)

    for fname in args.files:
        outroot = fname.replace('.fits', '_' + args.format)
        if args.format == 'fitsmod':
            convert_to_complete_fitszilla(fname, outroot)
        else:
            warnings.warn('Unknown output format')
