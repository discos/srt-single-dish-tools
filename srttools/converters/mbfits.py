from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
import os
import numpy as np
from srttools.io import mkdir_p, locations, read_data_fitszilla, \
    get_chan_columns, get_channel_feed
from srttools.utils import scantype, force_move_file


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


keywords_to_reset = ['11CD2F', '11CD2I', '11CD2J', '11CD2R', '11CD2S',
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
                                    'GROUPING.fits')) as grouping_template:
            grouping_template[1].data = grouping_template[1].data[:1]

            grouping_template.writeto(os.path.join(self.dirname, self.GROUPING),
                                      overwrite=True)

        self.SCAN = 'SCAN.fits'
        with fits.open(os.path.join(self.template_dir,
                                    'SCAN.fits')) as scan_template:
            scan_template[1].data['FEBE'][0] = 'EMPTY'

            scan_template.writeto(os.path.join(self.dirname, self.SCAN),
                                  overwrite=True)
        self.scan_count = 0
        self.date_obs = None

    def fill_in_summary(self, summaryfile):
        with fits.open(summaryfile) as hdul:
            header = hdul[0].header
            hdudict = dict(header.items())

        try:
            self.date_obs = Time(hdudict['DATE-OBS'])
        except KeyError:
            self.date_obs = Time(hdudict['DATE'])

        with fits.open(os.path.join(self.dirname, self.GROUPING)) as grouphdul:
            groupheader = grouphdul[0].header
            groupdict = dict(groupheader.items())
            for key in hdudict.keys():
                if key in groupdict:
                    groupheader[key] = hdudict[key]
            groupheader['RA'] = np.degrees(hdudict['RightAscension'])
            groupheader['DEC'] = np.degrees(hdudict['Declination'])
            groupheader['DATE-OBS'] = self.date_obs.value
            groupheader['MJD-OBS'] = self.date_obs.mjd
            grouphdul.writeto('tmp.fits', overwrite=True)

        force_move_file('tmp.fits', os.path.join(self.dirname, self.GROUPING))

        with fits.open(os.path.join(self.dirname, self.SCAN)) as scanhdul:
            scanheader = scanhdul[1].header
            scandict = dict(scanheader.items())
            for key in hdudict.keys():
                if key[:5] in ['NAXIS', 'PGCOU', 'GCOUN']:
                    continue
                if key in scandict:
                    scanheader[key] = hdudict[key]
            scanheader = reset_all_keywords(scanheader)
            scanheader['DATE-OBS'] = self.date_obs.value
            scanheader['MJD'] = self.date_obs.mjd
            scanhdul.writeto('tmp.fits', overwrite=True)

        force_move_file('tmp.fits', os.path.join(self.dirname, self.SCAN))

    def add_subscan(self, scanfile):
        subscan = read_data_fitszilla(scanfile)
        self.scan_count += 1

        chans = get_chan_columns(subscan)
        for ch in chans:
            feed = get_channel_feed(ch)
            polar = subscan[ch].meta['polarization']
            felabel = subscan.meta['receiver'] + '{}{}'.format(feed, polar)
            febe = felabel + '-' + subscan.meta['backend']

            datapar = os.path.join(self.template_dir, '1',
                                   'FLASH460L-XFFTS-DATAPAR.fits')
            with fits.open(datapar) as subs_par_template:
                n = len(subscan)

                ############ Update DATAPAR ############
                subs_par_template[1] = \
                    _copy_hdu_and_adapt_length(subs_par_template[1], n)

                newtable = Table(subs_par_template[1].data)
                time = Time(subscan['time'] * u.day, scale='utc', format='mjd')
                newtable['MJD'] = subscan['time']
                newtable['LST'][:] = \
                    time.sidereal_time('apparent',
                                       locations[subscan.meta['site']].lon).value
                newtable['INTEGTIM'][:] = \
                    subscan['Feed0_LCP'].meta['sample_rate']
                newtable['RA'] = subscan['ra']
                newtable['DEC'] = subscan['dec']
                newtable['LONGOFF'] = 0.
                newtable['LATOFF'] = 0.
                newtable['AZIMUTH'] = subscan['az']
                newtable['ELEVATIO'] = subscan['el']
                _, direction = scantype(subscan['ra'], subscan['dec'],
                                        el=subscan['el'], az=subscan['az'])
                direction_cut = \
                    direction.replace('<', '').replace('>', '').lower()
                if direction_cut in ['ra', 'dec']:
                    baslon, baslat = subscan['ra'], subscan['dec']
                elif direction_cut in ['el', 'az']:
                    baslon, baslat = subscan['az'], subscan['el']
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

                outdir = str(subscan.meta['SubScanID'])
                mkdir_p(os.path.join(self.dirname, outdir))
                new_datapar = os.path.join(outdir,
                                           febe + '-DATAPAR.fits')
                subs_par_template.writeto('tmp.fits', overwrite=True)

            force_move_file('tmp.fits',
                            os.path.join(self.dirname, new_datapar))

            arraydata = os.path.join(self.template_dir, '1',
                                     'FLASH460L-XFFTS-ARRAYDATA-1.fits')

            ############ Update ARRAYDATA ############
            with fits.open(arraydata) as subs_template:
                subs_template[1] = \
                    _copy_hdu_and_adapt_length(subs_template[1], n)

                subs_template[1].header = \
                    reset_all_keywords(subs_template[1].header)
                subs_template[1].header['SUBSNUM'] = subscan.meta['SubScanID']
                subs_template[1].header['DATE-OBS'] = self.date_obs.value
                subs_template[1].header['FEBE'] = febe
                subs_template[1].header['BASEBAND'] = 1
                subs_template[1].header['CHANNELS'] = subscan.meta['channels']

                newtable = Table(subs_template[1].data)
                newtable['MJD'] = subscan['time']
                newtable['DATA'] = subscan['Feed0_LCP']
                newhdu = fits.table_to_hdu(newtable)
                subs_template[1].data = newhdu.data

                new_sub = \
                    os.path.join(outdir, febe + '-ARRAYDATA-1.fits')
                subs_template.writeto('tmp.fits', overwrite=True)

            force_move_file('tmp.fits', os.path.join(self.dirname, new_sub))

            # Finally, update GROUPING file
            with fits.open(os.path.join(self.dirname,
                                        self.GROUPING)) as grouping:
                newtable = Table(grouping[1].data)
                if febe not in self.FEBE:

                    nfebe = len(list(self.FEBE.keys()))
                    new_febe = self.add_febe(subscan, ch, febe)

                    grouping[0].header['FEBE{}'.format(nfebe)] = febe
                    grouping[0].header['FREQ{}'.format(nfebe)] = \
                        subscan[ch].meta['frequency'] * 1e6
                    grouping[0].header['BWID{}'.format(nfebe)] = \
                        subscan[ch].meta['bandwidth'] * 1e6
                    grouping[0].header['LINE{}'.format(nfebe)] = ''

                    newtable.add_row([2, new_febe, 'URL', 'FEBEPAR-MBFITS',
                                      -999, febe, -999])
                    self.FEBE[febe] = new_febe

                newtable.add_row([2, new_datapar, 'URL', 'DATAPAR-MBFITS',
                                  -999, febe, -999])
                newtable.add_row([2, new_sub, 'URL', 'ARRAYDATA-MBFITS',
                                  self.scan_count, febe, 1])
                new_hdu = fits.table_to_hdu(newtable)
                grouping[1].data = new_hdu.data
                grouping[0].header['INSTRUME'] = subscan[ch].meta['backend']
                grouping[0].header['TELESCOP'] = subscan.meta['site']

                grouping.writeto('tmp.fits', overwrite=True)

            force_move_file('tmp.fits',
                            os.path.join(self.dirname, self.GROUPING))

            if self.test:
                break

    def add_febe(self, subscan, channel, febe):
        feed = get_channel_feed(channel)
        meta = subscan[channel].meta
        polar = meta['polarization']
        polar_code = polar[0]
        if polar_code == 'H':
            polar_code = 'X'
        elif polar_code == 'V':
            polar_code = 'Y'

        febe_name = febe + '-FEBEPAR.fits'

        with fits.open(
                os.path.join(self.template_dir,
                             'FLASH460L-XFFTS-FEBEPAR.fits')) as febe_template:

            febedata = Table(febe_template[1].data)
            febedata['USEBAND'] = [[1]]
            febedata['USEFEED'] = [[feed]]
            febedata['BESECTS'] = [[0]]
            febedata['FEEDTYPE'] = [[1]]
            febedata['POLTY'][:] = [polar_code + polar_code]
            febedata['POLA'][:] = [[0., 0.]]
            new_hdu = fits.table_to_hdu(febedata)
            febe_template.data = new_hdu.data
            # TODO: fill in the information given in the subscan[ch]

            new_febe = os.path.join(self.dirname, febe_name)

            febe_template[1].header['DATE-OBS'] = self.date_obs.value
            febe_template[1].header['FEBE'] = febe
            febe_template[1].header = \
                reset_all_keywords(febe_template[1].header)

            febe_template.writeto('tmp.fits', overwrite=True)
        force_move_file('tmp.fits', new_febe)

        with fits.open(os.path.join(self.dirname, self.SCAN)) as scan:
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

        return new_febe
