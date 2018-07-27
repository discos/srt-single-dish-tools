from astropy.io import fits
from astropy.time import Time
import astropy.units as u
import astropy.constants as c
import os
import numpy as np
from srttools.io import mkdir_p, locations, read_data_fitszilla, \
    get_chan_columns, classify_chan_columns, interpret_chan_name
import glob


model_primary_header = """
SIMPLE  =                    T
BITPIX  =                    8
NAXIS   =                    0
EXTEND  =                    T
BLOCKED =                    T
ORIGIN  = 'SRT'
CREATOR = '      '
END
"""

model_header = """
XTENSION= 'BINTABLE'           / binary table extension
BITPIX  =                    8 / 8-bit bytes
NAXIS   =                    2 / 2-dimensional binary table
NAXIS1  =                  280 / width of table in bytes
NAXIS2  =                   32 / number of rows in table
PCOUNT  =              2621440 / size of special data area
GCOUNT  =                    1 / one data group (required keyword)
TFIELDS =                   39 / number of fields in each row
EXTNAME = 'SINGLE DISH'        / name of this binary table extension
NMATRIX =                    1 / Number of DATA arrays
OBSERVER= ''                   / Observer name(s)
PROJID  = ''                   / Project name
TELESCOP= 'SRT     '           / Telescope name
OBSGEO-X=      4865182.7660    / [m] Antenna ITRF X-coordinate
OBSGEO-Y=      791922.6890     / [m] Antenna ITRF Y-coordinate
OBSGEO-Z=      4035137.1740    / [m] Antenna ITRF Z-coordinate
TTYPE1  = 'SCAN    '           / label for field
TFORM1  = '1I      '           / format of field
TTYPE2  = 'CYCLE   '           / label for field
TFORM2  = '1J      '           / format of field
TTYPE3  = 'DATE-OBS'           / label for field
TFORM3  = '10A     '           / format of field
TTYPE4  = 'TIME    '           / label for field
TFORM4  = '1D      '           / format of field
TUNIT4  = 's       '           / units of field
TTYPE5  = 'EXPOSURE'           / label for field
TFORM5  = '1E      '           / format of field
TUNIT5  = 's       '           / units of field
TTYPE6  = 'OBJECT  '           / label for field
TFORM6  = '16A     '           / format of field
TTYPE7  = 'OBJ-RA  '           / label for field
TFORM7  = '1D      '           / format of field
TUNIT7  = 'deg     '           / units of field
TTYPE8  = 'OBJ-DEC '           / label for field
TFORM8  = '1D      '           / format of field
TUNIT8  = 'deg     '           / units of field
TTYPE9  = 'RESTFRQ '           / label for field
TFORM9  = '1D      '           / format of field
TUNIT9  = 'Hz      '           / units of field
TTYPE10 = 'OBSMODE '           / label for field
TFORM10 = '16A     '           / format of field
TTYPE11 = 'BEAM    '           / label for field
TFORM11 = '1I      '           / format of field
TTYPE12 = 'IF      '           / label for field
TFORM12 = '1I      '           / format of field
TTYPE13 = 'FREQRES '           / label for field
TFORM13 = '1D      '           / format of field
TUNIT13 = 'Hz      '           / units of field
TTYPE14 = 'BANDWID '           / label for field
TFORM14 = '1D      '           / format of field
TUNIT14 = 'Hz      '           / units of field
CTYPE1  = 'FREQ    '           / DATA array axis 1: frequency in Hz.
TTYPE15 = 'CRPIX1  '           / label for field
TFORM15 = '1E      '           / format of field
TTYPE16 = 'CRVAL1  '           / label for field
TFORM16 = '1D      '           / format of field
TUNIT16 = 'Hz      '           / units of field
TTYPE17 = 'CDELT1  '           / label for field
TFORM17 = '1D      '           / format of field
TUNIT17 = 'Hz      '           / units of field
CTYPE2  = 'STOKES  '           / DATA array axis 2: polarization code
CRPIX2  =              1.0E+00 / Polarization code reference pixel
CRVAL2  =             -5.0E+00 / Polarization code at reference pixel (XX)
CDELT2  =             -1.0E+00 / Polarization code axis increment
CTYPE3  = 'RA      '           / DATA array axis 3 (degenerate): RA (mid-int)
CRPIX3  =              1.0E+00 / RA reference pixel
TTYPE18 = 'CRVAL3  '           / label for field
TFORM18 = '1D      '           / format of field
TUNIT18 = 'deg     '           / units of field
CDELT3  =             -1.0E+00 / RA axis increment
CTYPE4  = 'DEC     '           / DATA array axis 4 (degenerate): Dec (mid-int)
CRPIX4  =              1.0E+00 / Dec reference pixel
TTYPE19 = 'CRVAL4  '           / label for field
TFORM19 = '1D      '           / format of field
TUNIT19 = 'deg     '           / units of field
CDELT4  =              1.0E+00 / Dec axis increment
TTYPE20 = 'SCANRATE'           / label for field
TFORM20 = '2E      '           / format of field
TUNIT20 = 'deg/s   '           / units of field
SPECSYS = 'LSRK    '           / Doppler reference frame (transformed)
SSYSOBS = 'TOPOCENT'           / Doppler reference frame of observation
EQUINOX =              2.0E+03 / Equinox of equatorial coordinates
RADESYS = 'FK5     '           / Equatorial coordinate frame
TTYPE21 = 'TSYS    '           / label for field
TFORM21 = '2E      '           / format of field
TUNIT21 = 'K       '           / units of field
TTYPE22 = 'CALFCTR '           / label for field
TFORM22 = '2E      '           / format of field
TTYPE23 = 'DATA    '           / label for field
TFORM23 = '1PE(16384)'         / format of field
TTYPE24 = 'TDIM23  '           / label for field
TFORM24 = '16A     '           / format of field
TUNIT24 = 'K       '           / units of field
TTYPE25 = 'FLAGGED '           / label for field
TFORM25 = '1PB(16384)'         / format of field
TTYPE26 = 'TDIM25  '           / label for field
TFORM26 = '16A     '           / format of field
TTYPE27 = 'TCAL    '           / label for field
TFORM27 = '2E      '           / format of field
TUNIT27 = 'Jy      '           / units of field
TTYPE28 = 'TCALTIME'           / label for field
TFORM28 = '16A     '           / format of field
TTYPE29 = 'AZIMUTH '           / label for field
TFORM29 = '1E      '           / format of field
TUNIT29 = 'deg     '           / units of field
TTYPE30 = 'ELEVATIO'           / label for field
TFORM30 = '1E      '           / format of field
TUNIT30 = 'deg     '           / units of field
TTYPE31 = 'PARANGLE'           / label for field
TFORM31 = '1E      '           / format of field
TUNIT31 = 'deg     '           / units of field
TTYPE32 = 'FOCUSAXI'           / label for field
TFORM32 = '1E      '           / format of field
TUNIT32 = 'm       '           / units of field
TTYPE33 = 'FOCUSTAN'           / label for field
TFORM33 = '1E      '           / format of field
TUNIT33 = 'm       '           / units of field
TTYPE34 = 'FOCUSROT'           / label for field
TFORM34 = '1E      '           / format of field
TUNIT34 = 'deg     '           / units of field
TTYPE35 = 'TAMBIENT'           / label for field
TFORM35 = '1E      '           / format of field
TUNIT35 = 'C       '           / units of field
TTYPE36 = 'PRESSURE'           / label for field
TFORM36 = '1E      '           / format of field
TUNIT36 = 'Pa      '           / units of field
TTYPE37 = 'HUMIDITY'           / label for field
TFORM37 = '1E      '           / format of field
TUNIT37 = '%       '           / units of field
TTYPE38 = 'WINDSPEE'           / label for field
TFORM38 = '1E      '           / format of field
TUNIT38 = 'm/s     '           / units of field
TTYPE39 = 'WINDDIRE'           / label for field
TFORM39 = '1E      '           / format of field
TUNIT39 = 'deg     '           / units of field

"""

# LIST_TTYPE = ["SCAN", "CYCLE", "DATE-OBS", "TIME",
#               "EXPOSURE", "OBJECT", "OBJ_RA", "OBJ_DEC",
#               "RESTFREQ", "OBSMODE", "BEAM", "_IF",
#               "FREQRES", "BANDWID", "CRPIX1", "CRVAL1",
#               "CDELT1", "CRVAL3", "CRVAL4", "SCANRATE",
#               "TSYS", "CALFCTR", # "DATA", "FLAGGED",
#               "XCALFCTR", "TCAL", "TCALTIME", "AZIMUTH",
#               "ELEVATIO", "PARANGLE", "FOCUSAXI", "FOCUSTAN",
#               "FOCUSROT", "TAMBIENT", "PRESSURE", "HUMIDITY",
#               "WINDSPEE", "WINDDIRE"]
#
# LIST_TFORM = ["I", "I", "10A", "D",
#               "E", "16A", "D", "D",
#               "D", "16A", "I", "I",
#               "D", "D", "E ", "D",
#               "D", "D", "D", "2E",
#               "2E", "2E", # tformat, tformat2,
#               "2E", "2E", "16A", "E",
#               "E", "E", "E", "E",
#               "E ", "E", "E", "E",
#               "E", "E"]
#
# LIST_TUNIT = [""] * len(LIST_TFORM)


def get_data_description_from_model_header(data_format=None):
    header = fits.Header.fromstring(model_header, sep='\n')
    num = []
    list_ttype = []
    list_tform = []
    list_tunit = []
    headerdict = dict(header)
    for k in header:
        if not k.startswith('TTYPE'):
            continue
        n = k.replace('TTYPE', '')
        unit = ""
        if 'TUNIT' + n in headerdict:
            unit = header['TUNIT' + n]
        tform = header['TFORM' + n]
        if header[k] in ['DATA', 'FLAGGED'] and data_format is not None:
            tform = "{}D".format(np.product(data_format))
        num.append(int(n))
        list_ttype.append(header[k])
        list_tform.append(tform)
        list_tunit.append(unit)
    num = np.asarray(num)
    list_ttype = np.asarray(list_ttype)
    list_tform = np.asarray(list_tform)
    list_tunit = np.asarray(list_tunit)
    order = np.argsort(num)
    list_ttype = list_ttype[order]
    list_tform = list_tform[order]
    list_tunit = list_tunit[order]
    return list_ttype, list_tform, list_tunit


def get_model_HDUlist(data_format, length=1, **kwargs):
    """Produce a model CLASS-compatible HDUlist."""
    cols = []

    list_ttype, list_tform, list_tunit = \
        get_data_description_from_model_header(data_format)

    for ttype, tform, tunit in zip(list_ttype, list_tform, list_tunit):
        newcol = fits.Column(name=ttype, format=tform, unit=tunit,
                             array=np.zeros(length))
        cols.append(newcol)

    coldefs = fits.ColDefs(cols)

    hdu = fits.BinTableHDU.from_columns(
        coldefs, header=fits.Header.fromstring(model_header, sep='\n'),
        name='SINGLE DISH', **kwargs)

    primary_hdu = fits.PrimaryHDU(
        header=fits.Header.fromstring(model_primary_header, sep='\n'))
    return fits.HDUList([primary_hdu, hdu])


class SDFITS_creator():
    """SDFITS converter"""
    def __init__(self, dirname, scandir=None, average=True, use_calon=False,
                 test=False):
        """Initialization.

        Initialization is easy. If scandir is given, the conversion is
        done right away.

        Parameters
        ----------
        dirname : str
            Output directory for products

        Other Parameters
        ----------------
        scandir : str
            Input data directory (to be clear, the directory containing a set
            of subscans plus a summary.fits file)
        average : bool, default True
            Average all spectra of a given configuration?
        use_calon : bool, default False
            If False, only the OFF + CAL is used for the calibration. If True,
            Also the ON + CAL is used and the calibration constant is averaged
            with that obtained through OFF + CAL.
        test : bool
            Only use for unit tests
        """
        self.dirname = dirname
        self.test = test
        mkdir_p(dirname)
        self.summary = {}
        self.tables = {}
        self.average = average
        if scandir is not None:
            self.get_scan(scandir, average=average)
            self.write_tables_to_disk()

    def fill_in_summary(self, summaryfile):
        """Fill in the information contained in the summary.fits file."""
        with fits.open(summaryfile) as hdul:
            self.summary.update(hdul[0].header)

    def get_scan(self, scandir, average=False):
        """Treat the data and produce the output, uncalibrated files.

        Fills in the `self.tables` attribute with a dictionary of HDU lists
        containing a primary header and a MATRIX extension in CLASS-compatible
        FITS format

        Parameters
        ----------
        scandir : str
            Input data directory (to be clear, the directory containing a set
            of subscans plus a summary.fits file)

        Other Parameters
        ----------------
        average : bool, default True
            Average all spectra of a given configuration?

        Returns
        -------
        tables
        """
        scandir = scandir.rstrip('/')
        fname = os.path.join(scandir, 'summary.fits')
        self.fill_in_summary(fname)
        for fname in sorted(glob.glob(os.path.join(scandir, '*.fits'))):
            if 'summary' in fname:
                continue
            subscan = read_data_fitszilla(fname)
            location = locations[subscan.meta['site']]
            times = Time(subscan['time'] * u.day, format='mjd', scale='utc',
                         location=location)
            date_col = [t.strftime('%d/%m/%y') for t in times.to_datetime()]

            # Different from CLASS converter - here we take seconds from the
            # First day (when data span multiple dats)
            ut_col = (times.mjd - np.floor(times.mjd[0])) * 86400

            allcolumns = get_chan_columns(subscan)
            channels = \
                [subscan[ch].meta['channels'] for ch in allcolumns]
            if not len(set(channels)) == 1:
                raise ValueError("Only files with the same number of spectral "
                                 "bins in each channel are supported. Please "
                                 "report")
            classif = classify_chan_columns(allcolumns)
            feeds = list(classif.keys())

            for f in feeds:
                azimuth = subscan['az'][:, f].to(u.deg).value
                elevation = subscan['el'][:, f].to(u.deg).value
                crval3 = subscan['ra'][:, f].to(u.deg).value
                crval4 = subscan['dec'][:, f].to(u.deg).value

                columns_allbase = [a for a in allcolumns
                                   if a.startswith('Feed{}'.format(f))]

                basebands = \
                    [interpret_chan_name(ch)[2] for ch in columns_allbase]

                for baseband in basebands:
                    if baseband is None:
                        baseband = 0
                        columns = columns_allbase
                    else:
                        columns = [ch for ch in columns_allbase if
                                   ch.endswith('{}'.format(baseband))]
                    ncol = len(columns)

                    data_matrix = np.stack([subscan[ch] for ch in columns],
                                            axis=1)
                    shape = data_matrix.shape
                    if len(shape) == 3:
                        data_matrix = \
                            data_matrix.reshape((shape[0], shape[1]*shape[2]))

                    array = subscan[columns[0]]

                    newhdu = \
                        get_model_HDUlist(data_format=(channels[0], ncol),
                                          length=len(array))

                    data = newhdu[1].data

                    nbin = subscan.meta['channels']

                    bandwidth = array.meta['bandwidth']
                    restfreq_label = 'RESTFREQ{}'.format(baseband + 1)
                    if restfreq_label not in self.summary:
                        restfreq_label = 'RESTFREQ1'
                    restfreq = self.summary[restfreq_label] * u.MHz
                    #
                    # data['RESTFREQ'] = restfreq.to(u.Hz).value
                    data['EXPOSURE'] = \
                        array.meta['integration_time'].value
                    data['TIME'] = ut_col

                    data['TSYS'] = 1
                    df = (bandwidth / nbin).to('Hz')
                    data['CDELT1'] = df
                    deltav = - df / restfreq * c.c
                    data['FREQRES'] = deltav.to('m/s').value

                    data['TDIM23'] = str(data_matrix[0].shape)
                    data['TDIM25'] = str(data_matrix[0].shape)
                    data['DATA'] = data_matrix

                    data['OBJECT'] = subscan.meta['SOURCE']
                    data['AZIMUTH'] = azimuth
                    data['ELEVATIO'] = elevation
                    data['CRPIX1'] = nbin // 2 + 1
                    data['CRVAL3'] = crval3
                    data['CRVAL4'] = crval4

                    data['PARANGLE'] = subscan['par_angle']
                    data['FOCUSROT'] = subscan['derot_angle']
                    data['CRVAL4'] = crval4
                    weather = subscan['weather']
                    data["HUMIDITY"] = weather[:, 0]
                    data["TAMBIENT"] = weather[:, 1]
                    data["PRESSURE"] = weather[:, 2]
                    data["BEAM"] = f

                    header = newhdu[1].header
                    header['TELESCOP'] = subscan.meta['site']
                    header['OBSERVER'] = subscan.meta['OBSERVER']
                    header['OBSGEO-X'] = \
                        locations[subscan.meta['site']].x.to('m').value
                    header['OBSGEO-Y'] = \
                        locations[subscan.meta['site']].y.to('m').value
                    header['OBSGEO-Z'] = \
                        locations[subscan.meta['site']].z.to('m').value

                    header['CTYPE1'] = "FREQ"
                    header['CRVAL'] = 0
                    header['CRVAL3'] = \
                        np.mean(subscan['ra'][:, f].to(u.deg).value)
                    header['CRVAL4'] = \
                        np.mean(subscan['dec'][:, f].to(u.deg).value)
                    header['LINE'] = subscan.meta['SOURCE']
                    header['DATE-OBS'] = date_col[0]

                    header['DATE-RED'] = \
                        Time.now().to_datetime().strftime('%d/%m/%y')
                    header['LINE'] = \
                        "FEED{}-{:3.3f}-MHz".format(f,
                                                    bandwidth.to('MHz').value)
                    header['CDELT1'] = df.to('Hz').value
                    header['CDELT3'] = \
                        subscan.meta["ra_offset"].to(u.deg).value
                    header['CDELT4'] = \
                        subscan.meta["dec_offset"].to(u.deg).value
                    header['RESTFREQ'] = restfreq.to(u.Hz).value
                    header['MAXIS1'] = channels[0]

                    filekey = \
                        os.path.basename(scandir) + \
                            '_all_feed{}_bband{}'.format(f, baseband)

                    if filekey in list(self.tables.keys()):
                        hdul = self.tables[filekey]
                        nrows1, nrows2 = len(hdul[1].data), len(data)
                        nrows = nrows1 + nrows2
                        newhdu = fits.BinTableHDU.from_columns(hdul[1].columns,
                                                               nrows=nrows)
                        for col in hdul[1].columns:
                            name = col.name
                            newhdu.data[name][:nrows1] = hdul[1].data[name]
                            newhdu.data[name][nrows1:] = data[name]
                        hdul[1].data = newhdu.data
                    else:
                        self.tables[filekey] = newhdu

        return self.tables

    def write_tables_to_disk(self):
        """Write all HDU lists produced until now in separate FITS files."""
        for (filekey, table) in self.tables.items():
            outfile = os.path.join(self.dirname, '{}.fits'.format(filekey))
            table.writeto(outfile, overwrite=True)
