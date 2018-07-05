from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
import os
import numpy as np
from srttools.io import mkdir_p, locations, read_data_fitszilla, \
    get_chan_columns, classify_chan_columns
from srttools.utils import scantype, force_move_file, minmax, median_diff
from srttools.fit import detrend_spectroscopic_data
import warnings


def default_scan_info_table():
    return Table(names=['scan_id', 'start', 'stop',
                        'ra_min', 'ra_max', 'ra_d',
                        'dec_min', 'dec_max', 'dec_d',
                        'az_min', 'az_max', 'az_d',
                        'el_min', 'el_max', 'el_d',
                        'glon_min', 'glon_max', 'glon_d',
                        'glat_min', 'glat_max', 'glat_d',
                        'is_skydip', 'kind', 'direction'],

                 dtype=[int, float, float,
                        float, float, float, float, float, float,
                        float, float, float, float, float, float,
                        float, float, float, float, float, float,
                        bool, 'S10', 'S5'])


def get_subscan_info(subscan):
    info = default_scan_info_table()
    scan_id = subscan.meta['SubScanID']
    start, stop = minmax(subscan['time'])
    ramin, ramax = minmax(subscan['ra'])
    decmin, decmax = minmax(subscan['dec'])
    azmin, azmax = minmax(subscan['az'])
    elmin, elmax = minmax(subscan['el'])
    is_skydip = subscan.meta['is_skydip']

    d_ra = median_diff(subscan['ra'])
    d_dec = median_diff(subscan['dec'])
    d_az = median_diff(subscan['az'])
    d_el = median_diff(subscan['el'])

    ravar = (ramax - ramin) * np.cos(np.mean((decmin, decmax)))
    decvar = decmax - decmin
    azvar = (azmax - azmin) * np.cos(np.mean((elmin, elmax)))
    elvar = elmax - elmin

    tot_eq = np.sqrt(ravar ** 2 + decvar ** 2)
    tot_hor = np.sqrt(elvar ** 2 + azvar ** 2)
    ravar /= tot_eq
    decvar /= tot_hor

    directions = np.array(['ra', 'dec', 'az', 'el'])
    allvars = np.array([ravar, decvar, azvar, elvar])

    if tot_eq > 2 * tot_hor:
        kind = 'point'
        direction = ''
    else:
        kind = 'line'
        direction = directions[np.argmax(allvars)]

    info.add_row([scan_id, start, stop,
                  ramin, ramax, d_ra, decmin, decmax, d_dec,
                  azmin, azmax, d_az, elmin, elmax, d_el,
                  0, 0, 0, 0, 0, 0, is_skydip, kind, direction])

    return info


def format_direction(direction):
    """
    Examples
    --------
    >>> format_direction('ra')
    'ra'
    >>> format_direction('el')
    'alat'
    >>> format_direction('az')
    'alon'
    """
    lowerdir = direction.lower()
    if lowerdir == 'el':
        return 'alat'
    elif lowerdir == 'az':
        return 'alon'
    return direction


def get_observing_strategy_from_subscan_info(info):
    """Get observing strategy from subscan information."""
    kinds = info['kind']
    skydips = info['is_skydip']
    lines = info[kinds == 'line']
    points = info[kinds == 'point']
    ctype = 'RA/DEC'
    durations = (info['stop'] - info['start']) * 86400

    xspc, yspc = (0, 0)

    zigzag = False

    stype = 'MAP'
    direction = 'Unkn'
    length = 0

    if np.all(skydips):
        stype = 'SKYDIP'
        mode = 'OTF'
        geom = 'LINE'
        direction = 'ALAT'
    elif len(lines) > len(points):
        mode = 'OTF'
        ra_lines = lines[lines['direction'] == 'ra']
        dec_lines = lines[lines['direction'] == 'dec']
        az_lines = lines[lines['direction'] == 'az']
        el_lines = lines[lines['direction'] == 'el']

        directions = np.array(['ra', 'dec', 'az', 'el'])
        nsub = np.array([len(lines[lines['direction'] == d])
                         for d in directions])

        direction = directions[np.argmax(nsub)]
        if direction in ['ra', 'dec']:
            lon_lines, dlon, lat_lines, dlat = ra_lines, 'ra', dec_lines, 'dec'
        elif direction in ['az', 'el']:
            lon_lines, dlon, lat_lines, dlat = az_lines, 'az', el_lines, 'el'
        else:
            raise ValueError('Unknown scan direction')

        ctype = format_direction(dlon) + '/' + format_direction(dlat)
        sample_dist_lon = lon_lines[dlon + '_d']
        sample_dist_lat = lon_lines[dlat + '_d']

        if len(lon_lines) == len(lat_lines):
            geom = 'CROSS'
            zigzag = True
            length = \
                np.median(lon_lines[dlon + '_max'] - lon_lines[dlon + '_min'])
        elif len(lon_lines) > len(lat_lines):
            geom = 'LINE'
            # if we see an inversion of direction, set zigzag to True
            zigzag = np.any(sample_dist_lon[:-1] * sample_dist_lon[1:] < 0)
            length = \
                np.median(lon_lines[dlon + '_max'] - lon_lines[dlon + '_min'])
            direction = format_direction(dlon)
            xspc = 0
            yspc = median_diff(info[dlat + '_min'], sorting=True)
        else:
            geom = 'LINE'
            zigzag = np.any(sample_dist_lat[:-1] * sample_dist_lat[1:] < 0)
            length = \
                np.median(lat_lines[dlat + '_max'] - lat_lines[dlat + '_min'])
            direction = format_direction(dlat)
            yspc = 0
            xspc = median_diff(info[dlon + '_min'], sorting=True)

    else:
        mode = 'RASTER'
        geom = 'SINGLE'

    results = type('results', (), {})()
    results.mode = mode
    results.geom = geom
    results.sep = (xspc, yspc)
    results.zigzag = zigzag
    results.length = length
    results.type = stype
    results.ctype = ctype
    results.stype = stype
    results.scanvel = length / np.median(durations)
    results.direction = direction
    results.nobs = len(info['scan_id'])
    results.scantime = np.median(durations)
    return results


def _copy_hdu_and_adapt_length(hdu, length):
    data = hdu.data
    columns = []
    for col in data.columns:
        newvals = [data[col.name][0]] * length
        newcol = fits.Column(name=col.name, array=newvals,
                             format=col.format)
        columns.append(newcol)
    newhdu = fits.BinTableHDU.from_columns(columns)
    newhdu.header = hdu.header
    return newhdu


keywords_to_reset = [
    '11CD2F', '11CD2I', '11CD2J', '11CD2R', '11CD2S',
    '1CRPX2F', '1CRPX2I', '1CRPX2J', '1CRPX2R', '1CRPX2S', '1CRVL2F',
    '1CRVL2I', '1CRVL2J', '1CRVL2R', '1CRVL2S', '1CTYP2F', '1CTYP2I',
    '1CTYP2J', '1CTYP2R', '1CTYP2S', '1CUNI2F', '1CUNI2I', '1CUNI2J',
    '1CUNI2R', '1CUNI2S', '1SOBS2F', '1SOBS2I', '1SOBS2J', '1SOBS2R',
    '1SOBS2S', '1SPEC2F', '1SPEC2I', '1SPEC2J', '1SPEC2R', '1SPEC2S',
    '1VSOU2R', 'AN', 'ANRX', 'AW', 'AWRX', 'BANDWID', 'BLATOBJ', 'BLONGOBJ',
    'CA', 'CARX', 'DEWCABIN', 'DEWRTMOD', 'DEWUSER', 'DEWZERO', 'DISTANCE',
    'ECCENTR', 'FDELTACA', 'FDELTAIA', 'FDELTAIE', 'FDELTAX', 'FDELTAXT',
    'FDELTAY', 'FDELTAYT', 'FDELTAZ', 'FDELTAZT', 'FDTYPCOD', 'FEBEBAND',
    'FEBEFEED', 'FEGAIN', 'FREQRES', 'FRTHRWHI', 'FRTHRWLO', 'GRPID1',
    'GRPLC1', 'HACA', 'HACA2', 'HACA2RX', 'HACA3', 'HACA3RX', 'HACARX',
    'HASA', 'HASA2', 'HASA2RX', 'HASARX', 'HECA2', 'HECA2RX', 'HECA3',
    'HECA3RX', 'HECE', 'HECE2', 'HECE2RX', 'HECE6', 'HECE6RX', 'HECERX',
    'HESA', 'HESA2', 'HESA2RX', 'HESA3', 'HESA3RX', 'HESA4', 'HESA4RX',
    'HESA5', 'HESA5RX', 'HESARX', 'HESE', 'HESERX', 'HSCA', 'HSCA2',
    'HSCA2RX', 'HSCA5', 'HSCA5RX', 'HSCARX', 'HSSA3', 'HSSA3RX', 'IA', 'IARX',
    'IE', 'IERX', 'INCLINAT', 'LATOBJ', 'LONGASC', 'LONGOBJ', 'LONGSTRN',
    'NFEBE', 'NOPTREFL', 'NPAE', 'NPAERX', 'NPHASES', 'NRX', 'NRXRX', 'NRY',
    'NRYRX', 'NUSEBAND', 'OMEGA', 'OPTPATH', 'ORBEPOCH', 'ORBEQNOX', 'PATLAT',
    'PATLONG', 'PDELTACA', 'PDELTAIA', 'PDELTAIE', 'PERIDATE', 'PERIDIST',
    'REFOFFX', 'REFOFFY', 'REF_ONLN', 'REF_POL', 'RESTFREQ', 'SBSEP',
    'SCANLEN', 'SCANLINE', 'SCANNUM', 'SCANPAR1', 'SCANPAR2', 'SCANROT',
    'SCANRPTS', 'SCANSKEW', 'SCANTIME', 'SCANXSPC', 'SCANXVEL', 'SCANYSPC',
    'SIDEBAND', 'SIG_ONLN', 'SIG_POL', 'SKYFREQ', 'SWTCHMOD', 'TBLANK',
    'TRANSITI', 'TSYNC', 'WCSNM2F', 'WCSNM2I', 'WCSNM2J', 'WCSNM2R',
    'WCSNM2S', 'WOBTHROW', 'WOBUSED']


def pack_data(scan, polar_dict, detrend=False):
    """Pack data into MBFITS-ready format

    Examples
    --------
    >>> scan = {'Feed0_LCP': np.arange(4), 'Feed0_RCP': np.arange(4, 8)}
    >>> polar = {'LCP': 'Feed0_LCP', 'RCP': 'Feed0_RCP'}
    >>> res = pack_data(scan, polar)
    >>> np.allclose(res, [[0, 4], [1, 5], [2, 6], [3, 7]])
    True
    >>> scan = {'Feed0_LCP': np.arange(2), 'Feed0_RCP': np.arange(2, 4),
    ...         'Feed0_Q': np.arange(4, 6), 'Feed0_U': np.arange(6, 8)}
    >>> polar = {'LCP': 'Feed0_LCP', 'RCP': 'Feed0_RCP', 'Q': 'Feed0_Q',
    ...          'U': 'Feed0_U'}
    >>> res = pack_data(scan, polar)
    >>> np.allclose(res, [[0, 2, 4, 6], [1, 3, 5, 7]])
    True
    >>> scan = {'Feed0_LCP': np.ones((2, 4)), 'Feed0_RCP': np.zeros((2, 4))}
    >>> polar = {'LCP': 'Feed0_LCP', 'RCP': 'Feed0_RCP'}
    >>> res = pack_data(scan, polar)
    >>> np.allclose(res, [[[ 1.,  1.,  1.,  1.], [ 0.,  0.,  0.,  0.]],
    ...                   [[ 1.,  1.,  1.,  1.], [ 0.,  0.,  0.,  0.]]])
    True
    """

    polar_list = list(polar_dict.keys())
    if 'LCP' in polar_list:
        data = [scan[polar_dict['LCP']], scan[polar_dict['RCP']]]
        try:
            data.append(scan[polar_dict['Q']])
            data.append(scan[polar_dict['U']])
        except KeyError:
            pass
    else:  # pragma: no cover
        raise ValueError('Polarization kind not implemented yet')

    if detrend:
        new_data = []
        for d in data:
            detr, _ = detrend_spectroscopic_data(0, d, 'als')
            new_data.append(detr)
        data = new_data
    return np.stack(data, axis=1)


def reset_all_keywords(header):
    """Set a specific list of keywords to zero or empty string.

    Examples
    --------
    >>> from astropy.io.fits import Header
    >>> h = Header({'SCANNUM': 5, 'OPTPATH': 'dafafa', 'a': 'blabla'})
    >>> h2 = reset_all_keywords(h)
    >>> h2['SCANNUM']
    0
    >>> h2['OPTPATH']
    ''
    >>> # This is not in the list of keywords to eliminate
    >>> h2['a']
    'blabla'
    """
    import six
    for key in keywords_to_reset:
        if key in header:
            if isinstance(header[key], six.string_types):
                header[key] = ''
            else:
                header[key] = type(header[key])(0)
    return header


class MBFITS_creator():
    def __init__(self, dirname, test=False):
        self.dirname = dirname
        self.test = test
        mkdir_p(dirname)
        curdir = os.path.dirname(__file__)
        datadir = os.path.join(curdir, '..', 'data')
        self.template_dir = os.path.join(datadir, 'mbfits_template')

        self.FEBE = {}

        self.GROUPING = 'GROUPING.fits'
        with fits.open(os.path.join(self.template_dir,
                                    'GROUPING.fits'),
                       memmap=False) as grouping_template:
            grouping_template[1].data = grouping_template[1].data[:1]

            grouping_template.writeto(
                os.path.join(self.dirname, self.GROUPING), overwrite=True)

        self.SCAN = 'SCAN.fits'
        with fits.open(os.path.join(self.template_dir,
                                    'SCAN.fits'),
                       memmap=False) as scan_template:
            scan_template[1].data['FEBE'][0] = 'EMPTY'

            scan_template.writeto(os.path.join(self.dirname, self.SCAN),
                                  overwrite=True)
        self.date_obs = Time.now()
        self.scan_info = default_scan_info_table()
        self.nfeeds = None
        self.ra = 0
        self.dec = 0
        self.site = None
        self.lst = 1e32

    def fill_in_summary(self, summaryfile):
        print('Loading {}'.format(summaryfile))
        with fits.open(summaryfile, memmap=False) as hdul:
            header = hdul[0].header
            hdudict = dict(header.items())

        self.ra = np.degrees(hdudict['RightAscension'])
        self.dec = np.degrees(hdudict['Declination'])
        self.restfreq = None
        if 'RESTFREQ1' in hdudict:
            self.resfreq = hdudict['RESTFREQ1']

        try:
            self.date_obs = Time(hdudict['DATE-OBS'])
        except KeyError:
            self.date_obs = Time(hdudict['DATE'])
        try:
            self.obsid = int(hdudict['OBSID'])
        except (KeyError, ValueError):
            self.obsid = 9999

        with fits.open(os.path.join(self.dirname, self.GROUPING),
                       memmap=False) as grouphdul:
            groupheader = grouphdul[0].header
            groupdict = dict(groupheader.items())
            for key in hdudict.keys():
                if key in groupdict:
                    groupheader[key] = hdudict[key]
            groupheader['RA'] = self.ra
            groupheader['DEC'] = self.dec
            groupheader['DATE-OBS'] = self.date_obs.value
            groupheader['MJD-OBS'] = self.date_obs.mjd
            groupheader['SCANNUM'] = self.obsid
            grouphdul.writeto('tmp.fits', overwrite=True)

        force_move_file('tmp.fits', os.path.join(self.dirname, self.GROUPING))

        with fits.open(os.path.join(self.dirname, self.SCAN),
                       memmap=False) as scanhdul:
            scanheader = reset_all_keywords(scanhdul[1].header)
            scandict = dict(scanheader.items())
            for key in hdudict.keys():
                if key[:5] in ['NAXIS', 'PGCOU', 'GCOUN']:
                    continue
                if key in scandict:
                    scanheader[key] = hdudict[key]
            # Todo: update with correct keywords
            scanheader['DATE-OBS'] = self.date_obs.value
            scanheader['MJD'] = self.date_obs.mjd
            scanheader['SCANNUM'] = self.obsid
            scanhdul.writeto('tmp.fits', overwrite=True)

        force_move_file('tmp.fits', os.path.join(self.dirname, self.SCAN))

    def add_subscan(self, scanfile, detrend=False):
        print('Loading {}'.format(scanfile))

        subscan = read_data_fitszilla(scanfile)
        subscan_info = get_subscan_info(subscan)

        self.scan_info.add_row(subscan_info[0])

        time = Time(subscan['time'] * u.day, scale='utc', format='mjd')
        if self.date_obs.mjd > time[0].mjd:
            self.date_obs = time[0]
        if self.site is None:
            self.site = subscan.meta['site']

        chans = get_chan_columns(subscan)

        combinations = classify_chan_columns(chans)

        if self.nfeeds is None:
            self.nfeeds = len(combinations.keys())

        for feed in combinations:
            felabel = subscan.meta['receiver'] + '{}'.format(feed)
            febe = felabel + '-' + subscan.meta['backend']

            datapar = os.path.join(self.template_dir, '1',
                                   'FLASH460L-XFFTS-DATAPAR.fits')
            with fits.open(datapar, memmap=False) as subs_par_template:
                n = len(subscan)
                # ------------- Update DATAPAR --------------
                subs_par_template[1] = \
                    _copy_hdu_and_adapt_length(subs_par_template[1], n)

                newtable = Table(subs_par_template[1].data)
                newtable['MJD'] = subscan['time']
                newtable['LST'][:] = \
                    time.sidereal_time('apparent',
                                       locations[subscan.meta['site']].lon
                                       ).value
                if newtable['LST'][0] < self.lst:
                    self.lst = newtable['LST'][0]
                newtable['INTEGTIM'][:] = \
                    subscan['Feed0_LCP'].meta['sample_rate']
                newtable['RA'] = subscan['ra'].to(u.deg)
                newtable['DEC'] = subscan['dec'].to(u.deg)
                newtable['AZIMUTH'] = subscan['az'].to(u.deg)
                newtable['ELEVATIO'] = subscan['el'].to(u.deg)
                _, direction = scantype(subscan['ra'], subscan['dec'],
                                        el=subscan['el'], az=subscan['az'])

                direction_cut = \
                    direction.replace('<', '').replace('>', '').lower()
                if direction_cut in ['ra', 'dec']:
                    baslon = subscan['ra'].to(u.deg)
                    baslat = subscan['dec'].to(u.deg)

                    yoff = baslat.value - self.dec
                    # GLS projection
                    xoff = \
                        (baslon.value - self.ra)
                    newtable['LONGOFF'] = xoff * np.cos(np.radians(self.dec))
                    newtable['LATOFF'] = yoff
                elif direction_cut in ['el', 'az']:
                    warnings.warn('AltAz projection not implemented properly')
                    baslon, baslat = \
                        subscan['az'].to(u.deg), subscan['el'].to(u.deg)
                    newtable['LONGOFF'] = 0 * u.deg
                    newtable['LATOFF'] = 0 * u.deg
                else:
                    raise ValueError('Unknown coordinates')

                newtable['CBASLONG'] = baslon
                newtable['CBASLAT'] = baslat
                newtable['BASLONG'] = baslon
                newtable['BASLAT'] = baslat

                newhdu = fits.table_to_hdu(newtable)
                subs_par_template[1].data = newhdu.data
                subs_par_template[1].header['DATE-OBS'] = \
                    time[0].fits.replace('(UTC)', '')
                subs_par_template[1].header['LST'] = newtable['LST'][0]
                subs_par_template[1].header['FEBE'] = febe
                subs_par_template[1].header['SCANDIR'] = \
                    format_direction(direction_cut).upper()
                subs_par_template[1].header['SCANNUM'] = self.obsid

                outdir = str(subscan.meta['SubScanID'])
                mkdir_p(os.path.join(self.dirname, outdir))
                new_datapar = os.path.join(outdir,
                                           febe + '-DATAPAR.fits')
                subs_par_template.writeto('tmp.fits', overwrite=True)

            force_move_file('tmp.fits',
                            os.path.join(self.dirname, new_datapar))

            arraydata = os.path.join(self.template_dir, '1',
                                     'FLASH460L-XFFTS-ARRAYDATA-1.fits')

            new_arraydata_rows = []
            bands = list(combinations[feed].keys())
            for baseband in combinations[feed]:
                nbands = np.max(bands)
                ch = list(combinations[feed][baseband].values())[0]

                packed_data = pack_data(subscan, combinations[feed][baseband],
                                        detrend=detrend)
                # ------------- Update ARRAYDATA -------------
                with fits.open(arraydata, memmap=False) as subs_template:
                    subs_template[1] = \
                        _copy_hdu_and_adapt_length(subs_template[1], n)

                    new_header = \
                        reset_all_keywords(subs_template[1].header)

                    new_header['SCANNUM'] = self.obsid
                    new_header['SUBSNUM'] = subscan.meta['SubScanID']
                    new_header['DATE-OBS'] = self.date_obs.fits
                    new_header['FEBE'] = febe
                    new_header['BASEBAND'] = baseband
                    new_header['NUSEBAND'] = nbands
                    new_header['CHANNELS'] = subscan.meta['channels']
                    new_header['SKYFREQ'] = \
                        subscan[ch].meta['frequency'].to('Hz').value
                    if self.restfreq is not None:
                        new_header['RESTFREQ'] = self.restfreq
                    else:
                        new_header['RESTFREQ'] = new_header['SKYFREQ']
                    bandwidth = subscan[ch].meta['bandwidth'].to('Hz').value
                    new_header['BANDWID'] = bandwidth

                    new_header['FREQRES'] = bandwidth / new_header['CHANNELS']

                    # Todo: check sideband
                    new_header['SIDEBAND'] = 'USB'
                    # Todo: check all these strange keywords. These are
                    # probably NOT the rest frequencies!
                    new_header['1CRVL2F'] = new_header['RESTFREQ']
                    new_header['1CRVL2S'] = new_header['RESTFREQ']
                    for i in ['1CRPX2S', '1CRPX2R', '1CRPX2F', '1CRPX2J']:
                        new_header[i] = (new_header['CHANNELS'] + 1) // 2

                    subs_template[1].header = new_header
                    newtable = Table(subs_template[1].data)
                    newtable['MJD'] = subscan['time']
                    newtable['DATA'] = packed_data
                    newhdu = fits.table_to_hdu(newtable)
                    subs_template[1].data = newhdu.data

                    subname = febe + '-ARRAYDATA-{}.fits'.format(baseband)
                    new_sub = \
                        os.path.join(outdir, subname)
                    subs_template.writeto('tmp.fits', overwrite=True)

                    new_arraydata_rows.append([2, new_sub, 'URL',
                                               'ARRAYDATA-MBFITS',
                                               subscan.meta['SubScanID'], febe,
                                               baseband])

                force_move_file('tmp.fits',
                                os.path.join(self.dirname, new_sub))

            # Finally, update GROUPING file
            with fits.open(os.path.join(self.dirname,
                                        self.GROUPING),
                           memmap=False) as grouping:
                newtable = Table(grouping[1].data)
                if febe not in self.FEBE:

                    nfebe = len(list(self.FEBE.keys()))
                    new_febe = self.add_febe(febe, combinations, feed,
                                             subscan[ch].meta,
                                             bands=bands)

                    grouping[0].header['FEBE{}'.format(nfebe)] = febe
                    grouping[0].header['FREQ{}'.format(nfebe)] = \
                        subscan[ch].meta['frequency'].to('Hz').value
                    grouping[0].header['BWID{}'.format(nfebe)] = \
                        subscan[ch].meta['bandwidth'].to('Hz').value
                    grouping[0].header['LINE{}'.format(nfebe)] = ''

                    newtable.add_row([2, new_febe, 'URL', 'FEBEPAR-MBFITS',
                                      -999, febe, -999])
                    self.FEBE[febe] = new_febe

                newtable.add_row([2, new_datapar, 'URL', 'DATAPAR-MBFITS',
                                  -999, febe, -999])

                for row in new_arraydata_rows:
                    newtable.add_row(row)
                new_hdu = fits.table_to_hdu(newtable)
                grouping[1].data = new_hdu.data
                grouping[0].header['INSTRUME'] = subscan[ch].meta['backend']
                grouping[0].header['TELESCOP'] = self.site

                grouping.writeto('tmp.fits', overwrite=True)

            force_move_file('tmp.fits',
                            os.path.join(self.dirname, self.GROUPING))

            if self.test:
                break

    def add_febe(self, febe, feed_info, feed, meta, bands=None):
        if bands is None:
            bands = [1]
        polar = 'N'
        polar_code = polar[0]

        febe_name = febe + '-FEBEPAR.fits'

        with fits.open(
                os.path.join(self.template_dir,
                             'FLASH460L-XFFTS-FEBEPAR.fits'),
                memmap=False) as febe_template:

            febe_template[1].header = \
                reset_all_keywords(febe_template[1].header)

            febedata = Table(febe_template[1].data)
            # FEBEFEED stores the total number of feeds for the receiver in
            # use.  A receiver outputting two polarisations counts as two
            # feeds.  For an array, count the total no.  of pixels, even if
            # not all in use.
            febedata['USEBAND'] = np.array([bands])
            febedata['NUSEFEED'] = np.array([[2]])
            febedata['USEFEED'] = \
                np.array([[feed * 2 + 1, feed * 2 + 2,
                           feed * 2 + 1, feed * 2 + 2]])
            febedata['BESECTS'] = np.array([[0]])
            febedata['FEEDTYPE'] = np.array([[1, 2, 3, 4]])
            febedata['POLTY'][:] = np.array([polar_code])
            febedata['POLA'][:] = np.array([[0., 0.]])
            new_hdu = fits.table_to_hdu(febedata)

            febe_template[1].data = new_hdu.data
            # TODO: fill in the information given in the subscan[ch]

            new_febe = os.path.join(self.dirname, febe_name)

            febe_template[1].header['DATE-OBS'] = self.date_obs.fits
            febe_template[1].header['FEBE'] = febe
            febe_template[1].header['FEBEFEED'] = self.nfeeds * 2
            febe_template[1].header['NUSEBAND'] = max(bands)
            febe_template[1].header['NPHASES'] = 1
            febe_template[1].header['SWTCHMOD'] = 'NONE'
            febe_template[1].header['SCANNUM'] = self.obsid

            if 'Q' in feed_info[feed][bands[0]].keys():
                febe_template[1].header['FDTYPCOD'] = '1:L, 2:R, 3:Q, 4:U'
            else:
                febe_template[1].header['FDTYPCOD'] = '1:L, 2:R'
            febe_template.writeto('tmp.fits', overwrite=True)
        force_move_file('tmp.fits', new_febe)

        with fits.open(os.path.join(self.dirname, self.SCAN),
                       memmap=False) as scan:
            newtable = Table(scan[1].data)

            if newtable['FEBE'][0].strip() == 'EMPTY':
                newtable['FEBE'][0] = febe
            else:
                newtable.add_row([febe])

            new_hdu = fits.table_to_hdu(newtable)
            scan[1].data = new_hdu.data
            scanheader = scan[1].header
            scanheader['SITELONG'] = np.degrees(meta['SiteLongitude'])
            scanheader['SITELAT'] = np.degrees(meta['SiteLatitude'])
            scanheader['SITEELEV'] = meta['SiteHeight']
            diameter = 64. if meta['site'].lower().strip() == 'srt' else 32.
            scanheader['DIAMETER'] = diameter
            scanheader['PROJID'] = meta['Project_Name']

            scan.writeto('tmp.fits', overwrite=True)
        force_move_file('tmp.fits', os.path.join(self.dirname, self.SCAN))

        return febe_name

    def update_scan_info(self):
        info = \
            get_observing_strategy_from_subscan_info(self.scan_info)

        with fits.open(os.path.join(self.dirname, self.SCAN),
                       memmap=False) as scanhdul:
            scanheader = scanhdul[1].header
            # Todo: update with correct keywords
            scanheader['CTYPE'] = info.ctype
            scanheader['CTYPE1'] = 'RA---GLS'
            scanheader['CTYPE2'] = 'DEC--GLS'
            scanheader['CRVAL1'] = self.ra
            scanheader['CRVAL2'] = self.dec
            scanheader['BLONGOBJ'] = self.ra
            scanheader['BLATOBJ'] = self.dec
            scanheader['LONGOBJ'] = self.ra if not info.ctype[0] == 'A' else 0
            scanheader['LATOBJ'] = self.dec if not info.ctype[0] == 'A' else 0
            scanheader['EQUINOX'] = 2000.
            scanheader['GRPLC1'] = 'GROUPING.fits'
            scanheader['LST'] = self.lst
            scanheader['LATPOLE'] = 90.
            scanheader['LONPOLE'] = 0.
            scanheader['PATLONG'] = 0
            scanheader['MOVEFRAM'] = False
            if info.ctype == 'ALON/ALAT':
                scanheader['WCSNAME'] = 'Absolute horizontal'
            scanheader['SCANTYPE'] = info.stype.upper()
            scanheader['SCANDIR'] = info.direction.upper()
            scanheader['SCANXVEL'] = info.scanvel
            scanheader['SCANTIME'] = info.scantime
            scanheader['SCANMODE'] = info.mode.upper()
            scanheader['SCANGEOM'] = info.geom.upper()
            scanheader['SCANLINE'] = 1
            scanheader['SCANLEN'] = np.degrees(info.length)

            scanheader['SCANYSPC'] = np.degrees(info.sep[1])
            scanheader['SCANXSPC'] = np.degrees(info.sep[0])
            scanheader['SCANPAR1'] = -999
            scanheader['SCANPAR2'] = -999
            scanheader['ZIGZAG'] = info.zigzag
            scanheader['PHASE1'] = 'sig'
            scanheader['PHASE2'] = 'sig'
            scanheader['NOBS'] = info.nobs
            scanheader['NSUBS'] = info.nobs
            scanheader['WOBCYCLE'] = 0.
            scanheader['WOBDIR'] = 'NONE'
            scanheader['WOBMODE'] = 'NONE'
            scanheader['WOBPATT'] = 'NONE'
            scanhdul.writeto('tmp.fits', overwrite=True)
        force_move_file('tmp.fits', os.path.join(self.dirname, self.SCAN))

    def wrap_up_file(self):
        import copy
        prihdu = fits.PrimaryHDU()
        with fits.open(os.path.join(self.dirname, self.GROUPING),
                       memmap=False) as grouhdl:
            prihdu.header = copy.deepcopy(grouhdl[0].header)
            file_list = list(zip(grouhdl[1].data['MEMBER_LOCATION'],
                                 grouhdl[1].data['EXTNAME'],
                                 grouhdl[1].data['FEBE']))

        hdulists = {}
        for febe in self.FEBE.keys():
            hdulists[febe] = fits.HDUList([prihdu])

            with fits.open(os.path.join(self.dirname, self.SCAN),
                           memmap=False) as scanhdul:
                scanhdul[1].data['FEBE'] = [febe]
                newhdu = type(scanhdul[1])()
                newhdu.data = scanhdul[1].data
                newhdu.header = scanhdul[1].header
                hdulists[febe].append(newhdu)

        for fname, ext, febe in file_list:
            if febe == '':
                continue
            with fits.open(os.path.join(self.dirname, fname),
                           memmap=False) as hl:
                newhdu = type(hl[ext])()
                newhdu.data = hl[ext].data
                newhdu.header = hl[ext].header
                hdulists[febe].append(newhdu)

        fnames = {}
        for febe, hdulist in hdulists.items():
            fname = self.dirname + '.' + febe + '.fits'
            hdulist.writeto(fname, overwrite=True)
            hdulist.close()
            fnames[febe] = fname

        return fnames
