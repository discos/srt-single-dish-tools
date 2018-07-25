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
XTENSION= 'BINTABLE'
BITPIX  =                    8         / Always 8.
NAXIS   =                    2         / Always 2: tables have 2 dimensions.
NAXIS1  =                 7275         / Number of bytes per row.
NAXIS2  =                    4         / Number of rows.
PCOUNT  =                    0         / Usually 0.
GCOUNT  =                    1         / Always 1.
TFIELDS =                   18         / Number of columns.
EXTNAME = 'SINGLE DISH'                   / Just a name, not important.
EXTVER  =                    1         / Always 1.
MAXIS   =                    4         / Number of dimensions in the data.
MAXIS1  =                 1793         / Dummy number of channels (see TTYPE1).
MAXIS2  =                    1         / 
MAXIS3  =                    1         / 
MAXIS4  =                    1         / 
CTYPE1  = 'FREQ    '                   / Dim1: freq => MAXIS1 = Nb channels.
CRVAL1  =  0.0000000000000E+00         / Frequency offset, always 0.
CDELT1  =  0.0000000000000E+00         / Frequency resolution [Hz].
CRPIX1  =  0.0000000000000E+00         / Dummy reference channel (see TTYPE18).
CTYPE3  = 'RA      '
CRVAL3  =  0.0000000000000E+00
CDELT3  =  0.0000000000000E+00
CRPIX3  =  0.0000000000000E+00
CTYPE4  = 'DEC     '
CRVAL4  =  0.0000000000000E+00
CDELT4  =  0.0000000000000E+00
CRPIX4  =  0.0000000000000E+00
SUBSCAN =                    1         / Subscan number.  Often 1.
LINE    = '            '               / Name of your line, up to 12 chars.
OBJECT  = '            '               / Name of your source, up to 12 chars.
TELESCOP= '            '               / Name of your telescope, up to 12 chars.
BANDWID =  0.000000000000               / 
DATE-OBS= '            '               / 
EXPOSURE=  0.000000000000               / 
TSYS    =  0.000000000000               / 
RESTFREQ=  0.0000000000000E+00         / Rest (signal) frequency at ref chan.
VELO-HEL=  0.0000000000000E+00         / Velocity at ref.  chan [m.s-1].
VELDEF  = 'RADI-LSR'                   / Type of velocity.
GAINIMAG=  0.0000000000000E+00         / Ratio Image/Signal.
BEAMEFF =  0.0000000000000E+00         / Beam efficiency.
FORWEFF =  0.0000000000000E+00         / Forward efficiency.
EPOCH   =  2.0000000000000E+03         / Epoch of coordinates.
DATE-RED= '15/07/97'                   / Date of reduction.
"""

LIST_TTYPE = ["SCAN", "CYCLE", "DATE-OBS", "TIME",
              "EXPOSURE", "OBJECT", "OBJ_RA", "OBJ_DEC",
              "RESTFREQ", "OBSMODE", "BEAM", "_IF",
              "FREQRES", "BANDWID", "CRPIX1", "CRVAL1",
              "CDELT1", "CRVAL3", "CRVAL4", "SCANRATE",
              "TSYS", "CALFCTR", # "DATA", "FLAGGED",
              "XCALFCTR", "TCAL", "TCALTIME", "AZIMUTH",
              "ELEVATIO", "PARANGLE", "FOCUSAXI", "FOCUSTAN",
              "FOCUSROT", "TAMBIENT", "PRESSURE", "HUMIDITY",
              "WINDSPEE", "WINDDIRE"]

LIST_TFORM = ["I", "I", "10A", "D",
              "E", "16A", "D", "D",
              "D", "16A", "I", "I",
              "D", "D", "E ", "D",
              "D", "D", "D", "2E",
              "2E", "2E", # tformat, tformat2,
              "2E", "2E", "16A", "E",
              "E", "E", "E", "E",
              "E ", "E", "E", "E",
              "E", "E"]

LIST_TUNIT = [""] * len(LIST_TFORM)


def get_model_HDUlist(additional_columns=None, **kwargs):
    """Produce a model CLASS-compatible HDUlist."""
    cols = []
    list_ttype = LIST_TTYPE
    list_tform = LIST_TFORM
    list_tunit = LIST_TUNIT

    for ttype, tform, tunit in zip(list_ttype, list_tform, list_tunit):
        newcol = fits.Column(name=ttype, format=tform, unit=tunit)
        cols.append(newcol)
    coldefs = fits.ColDefs(cols)
    if additional_columns is not None:
        coldefs += fits.ColDefs(additional_columns)

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

                    array = subscan[columns[0]]

                    newcol = fits.Column(array=data_matrix, name="DATA",
                                         unit="K",
                                         format="{}D".format(channels[0] * ncol))

                    newhdu = get_model_HDUlist(additional_columns=[newcol])

                    data = newhdu[1].data

                    nbin = subscan.meta['channels']

                    bandwidth = array.meta['bandwidth']
                    restfreq_label = 'RESTFREQ{}'.format(baseband + 1)
                    if restfreq_label not in self.summary:
                        restfreq_label = 'RESTFREQ1'
                    restfreq = self.summary[restfreq_label] * u.MHz

                    data['RESTFREQ'] = restfreq.to(u.Hz).value
                    data['EXPOSURE'] = \
                        array.meta['integration_time'].value
                    data['TIME'] = ut_col

                    data['TSYS'] = 1
                    df = (bandwidth / nbin).to('Hz')
                    data['CDELT1'] = df
                    deltav = - df / restfreq * c.c
                    data['FREQRES'] = deltav.to('m/s').value

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
                    header['CTYPE1'] = "FREQ"
                    header['CRVAL'] = 0
                    header['CRVAL3'] = \
                        np.mean(subscan['ra'][:, f].to(u.deg).value)
                    header['CRVAL4'] = \
                        np.mean(subscan['dec'][:, f].to(u.deg).value)
                    header['LINE'] = subscan.meta['SOURCE']
                    header['OBJECT'] = subscan.meta['SOURCE']
                    header['DATE-OBS'] = date_col[0]

                    header['SOURCE'] = subscan.meta['SOURCE']
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
