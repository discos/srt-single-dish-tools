from astropy.io import fits
from astropy.table import Table, vstack
from astropy.time import Time
import astropy.units as u
import astropy.constants as c
import os
import numpy as np
from srttools.io import mkdir_p, locations, read_data_fitszilla, \
    get_chan_columns, classify_chan_columns, interpret_chan_name
import glob
from ..utils import get_mH2O
from ..io import label_from_chan_name
from scipy.signal import medfilt
import copy
import warnings
import collections


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
EXTNAME = 'MATRIX  '                   / Just a name, not important.
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
CTYPE2  = 'RA      '
CRVAL2  =  0.0000000000000E+00
CDELT2  =  0.0000000000000E+00
CRPIX2  =  0.0000000000000E+00
CTYPE3  = 'DEC     '
CRVAL3  =  0.0000000000000E+00
CDELT3  =  0.0000000000000E+00
CRPIX3  =  0.0000000000000E+00
CTYPE4  = 'STOKES  '
CRVAL4  =  0.0000000000000E+00
CDELT4  =  0.0000000000000E+00
CRPIX4  =  0.0000000000000E+00
SUBSCAN =                    1         / Subscan number.  Often 1.
LINE    = '            '               / Name of your line, up to 12 chars.
OBJECT  = '            '               / Name of your source, up to 12 chars.
RESTFREQ=  0.0000000000000E+00         / Rest (signal) frequency at ref chan.
VELO-HEL=  0.0000000000000E+00         / Velocity at ref.  chan [m.s-1].
VELDEF  = 'RADI-LSR'                   / Type of velocity.
GAINIMAG=  0.0000000000000E+00         / Ratio Image/Signal.
BEAMEFF =  0.0000000000000E+00         / Beam efficiency.
FORWEFF =  0.0000000000000E+00         / Forward efficiency.
EPOCH   =  2.0000000000000E+03         / Epoch of coordinates.
DATE-RED= '15/07/97'                   / Date of reduction.
"""

LIST_TTYPE = \
    ["MJD",
     "MAXIS1", "SCAN", "TELESCOP", "TSYS",
     "IMAGFREQ", "DELTAV", "TAU-ATM", "MH2O",
     "TOUTSIDE", "PRESSURE", "CRVAL2", "CRVAL3",
     "ELEVATIO", "AZIMUTH", "DATE-OBS", "UT",
     "LST", "OBSTIME", "CRPIX1", "RESTFREQ",
     "OBJECT", "VELOCITY", "CDELT1", "CDELT2",
     "CDELT3", "LINE", "SIGNAL", "CAL_IS_ON"]

LIST_TFORM = \
    ["1D",
     "1J", "1J", "12A", "1E",
     "1E", "1E", "1E", "1E",
     "1E", "1E", "1E", "1E",
     "1E", "1E", "23A ", "1D",
     "1D", "1E", "1E", "1D",
     "12A", "1E", "1E", "1D",
     "1D", "12A", "1J", "1J"]

LIST_TUNIT = \
    ["d",
     " ", "", "", "K",
     "Hz", "m.s-1", "neper", "mm",
     "K", "hPa", "deg", "deg",
     "deg", "deg", "", "s",
     "s", "s", "", "Hz",
     "", "m.s-1", "Hz", "deg",
     "deg", "",  "", ""]


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
        name='MATRIX', **kwargs)

    primary_hdu = fits.PrimaryHDU(
        header=fits.Header.fromstring(model_primary_header, sep='\n'))
    return fits.HDUList([primary_hdu, hdu])


def create_variable_length_column(values, max_length=2048, name="SPECTRUM",
                                  unit="K"):
    """If we want to use variable length arrays, this is what we should do.

    Examples
    --------
    >>> col = create_variable_length_column([[1, 2]])
    >>> isinstance(col, fits.Column)
    True
    """
    format_str = "PJ({})".format(max_length)
    column = fits.Column(array=values, name=name, unit=unit,
                         format=format_str)
    return column


def on_or_off(subscan, feed):
    """Try to infer if a given subscan is ON or OFF."""
    is_on = False
    if 'SIGNAL' in subscan.meta and \
            subscan.meta['SIGNAL'] in ['SIGNAL', 'REFERENCE']:
        signal = subscan.meta['SIGNAL']
        # In nodding, feed 0 has
        if signal == 'SIGNAL' and feed == 0:
            is_on = True
        elif signal == 'REFERENCE' and feed != 0:
            is_on = True
    else:
        is_on = subscan.meta['az_offset'] > 1e-4 * u.rad
    return is_on


def cal_is_on(subscan):
    """Is the calibration mark on? Try to understand."""
    is_on = False
    if 'SUBSTYPE' in subscan.meta:
        is_on = subscan.meta['SUBSTYPE'] == 'CAL'
    elif 'flag_cal' in subscan.colnames and np.any(subscan['flag_cal'] == 1):
        is_on = subscan['flag_cal']

    return is_on


def find_cycles(table, list_of_keys):
    """Find cyclic patterns in table.

    Parameters
    ----------
    table : `astropy.table.Table` object, or compatible
        Input table
    list_of_keys : list
        List of keywords of the table to cycle on

    Examples
    --------
    >>> table = Table(data=[[0, 0, 1, 1, 0, 0, 1, 1],
    ...                     [0, 0, 0, 0, 0, 0, 0, 0]], names=['A', 'B'])
    >>> list_of_keys = ['A', 'B']
    >>> new_table = find_cycles(table, list_of_keys)
    >>> np.all(new_table['CYCLE'] == [0, 0, 0, 0, 1, 1, 1, 1])
    True
    >>> table = Table(data=[[0, 0, 1, 1, 0, 0, 1, 1],
    ...                     [0, 1, 0, 1, 0, 1, 0, 1]], names=['A', 'B'])
    >>> list_of_keys = ['A', 'B']
    >>> new_table = find_cycles(table, list_of_keys)
    >>> np.all(new_table['CYCLE'] == [0, 0, 0, 0, 1, 1, 1, 1])
    True
    >>> table = Table(data=[[0, 1, 0, 1], [1, 0, 1, 0]], names=['A', 'B'])
    >>> list_of_keys = ['A', 'B']
    >>> new_table = find_cycles(table, list_of_keys)
    >>> np.all(new_table['CYCLE'] == [0, 0, 1, 1])
    True
    """
    binary_values = 10 ** np.arange(len(list_of_keys), dtype=int)
    table['BINARY_COL'] = np.zeros(len(table), dtype=int)
    for i, k in enumerate(list_of_keys):
        good = (table[k] == 1)
        table['BINARY_COL'][good] += binary_values[i]

    table['CYCLE'] = np.zeros(len(table), dtype=int)
    cycle_counter = -1
    start_value = table['BINARY_COL'][0]
    last = -1
    for i, b in enumerate(table['BINARY_COL']):
        if b == start_value and b != last:
            cycle_counter += 1
        table['CYCLE'][i] = cycle_counter
        last = b
    return table


def normalize_on_off_cal(table, smooth=False, apply_cal=True, use_calon=False):
    """Do the actuall onoff/onoffcal calibration.

    The first passage is to combine the ON-source signal with the closest
    OFF-source signal (alternatively, SIGNAL/on-source vs REFERENCE/off-source)
    as (ON - OFF)/OFF.

    If a signal with calibration mark is present (CAL), and apply_cal is True,
    an additional calibration factor is applied. If CAL is applied only to
    the OFF/REFERENCE signal (let us call it OFFCAL), a single calibration
    factor is calculated: OFF/(OFFCAL - OFF). If the CAL signal is also applied
    to the ON signal, and `use_calon` is True, an additional factor
    ON/(ONCAL - ON) is calculated and averaged with the previous.

    Parameters
    ----------
    table : `astropy.table.Table` object
        The data table

    Other Parameters
    ----------------
    smooth : bool, default False
        Run a median filter on the reference data (the ones that go to the
        denominator) before calculating the ratio
    apply_cal : bool, default True
        If a CAL (calib. mark on) signal is present, use it!
    use_calon : bool, default False
        If False, only the OFF + CAL is used for the calibration. If True,
        Also the ON + CAL is used and the calibration constant is averaged
        with that obtained through OFF + CAL.
"""
    calibration_factor = 1
    unit = ""

    cal_on = table['CAL_IS_ON'] == 1
    onsource = table['SIGNAL'] == 1
    on_data = table[~cal_on & onsource]
    off_data = table[~cal_on & ~onsource]
    calon_data = table[cal_on & onsource]
    caloff_data = table[cal_on & ~onsource]

    newtable = Table(on_data)

    on = on_data['SPECTRUM']
    off = np.mean(off_data['SPECTRUM'], axis=0)
    calon = caloff = None

    if len(calon_data) > 0 and use_calon:
        calon = np.mean(calon_data['SPECTRUM'], axis=0)

    if len(caloff_data) > 0:
        caloff = np.mean(caloff_data['SPECTRUM'], axis=0)

    off_ref = off
    on_ref = on[0]
    if smooth:
        off_ref = medfilt(off_ref, 11)
        on_ref = medfilt(on_ref, 11)

    signal = copy.deepcopy(on)
    for i, o in enumerate(on):
        signal[i] = (o - off_ref) / off_ref

    if apply_cal:
        oncal = offcal = np.zeros_like(caloff)
        if caloff is not None:
            offcal = (caloff - off_ref) / off_ref
        if calon is not None:
            oncal = (calon - on_ref) / on_ref
        cal = np.array([oncal, offcal])
        good = cal != 0
        cal = cal[good]

        if len(cal) > 0:
            calibration_factor = 1 / np.mean(cal)
            unit = "K"
        else:
            return None, ""

    newtable['SPECTRUM'][:] = signal * calibration_factor

    return newtable, unit


class CLASSFITS_creator():
    """CLASS-compatible FITS creator obhject."""
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
            self.calibrate_all(use_calon)
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
            mjd_col = subscan['time']
            ut_col = (times.mjd - np.floor(times.mjd)) * 86400

            lsts = times.sidereal_time('apparent',
                                       locations[subscan.meta['site']].lon
                                       )
            lsts_col = lsts.value * u.hr

            mH2O = [get_mH2O(weather[1] + 273.15, weather[0])
                    for weather in subscan['weather']]
            lsts = lsts_col.to('s').value
            if average:
                date_col = date_col[0]
                ut_col = ut_col[0]
                mH2O = np.mean(mH2O)
                lsts = np.mean(lsts)
                mjd_col = mjd_col[0]

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
                crval2 = subscan['ra'][:, f].to(u.deg).value
                crval3 = subscan['dec'][:, f].to(u.deg).value

                columns = [a for a in allcolumns
                           if a.startswith('Feed{}'.format(f))]

                data_matrix = []
                for ch in columns:
                    array = subscan[ch]
                    if average:
                        array = Table(data=[[np.mean(array, axis=0)]],
                                      meta=array.meta)

                    data_matrix.extend(array)

                if average:
                    azimuth = np.mean(azimuth)
                    elevation = np.mean(elevation)
                    crval2 = np.mean(crval2)
                    crval3 = np.mean(crval3)

                newcol = fits.Column(array=data_matrix, name="SPECTRUM",
                                     unit="K",
                                     format="{}D".format(channels[0]))

                newhdu = get_model_HDUlist(additional_columns=[newcol])

                data = newhdu[1].data

                id0 = 0

                for ch in columns:
                    feed, polar, baseband = interpret_chan_name(ch)
                    if feed != f:
                        warnings.warn("Problem interpreting chan name: {} "
                                      "instead of {}".format(feed, f))
                    if baseband is None:
                        baseband = 1
                    array = subscan[ch]
                    if average:
                        length = len(array)
                        array = Table(data=[[np.mean(array, axis=0)]],
                                      meta=array.meta)
                        array.meta['integration_time'] *= length

                    length = len(array)
                    id1 = id0 + length
                    nbin = array.meta['channels']

                    bandwidth = array.meta['bandwidth']
                    restfreq_label = 'RESTFREQ{}'.format(baseband + 1)
                    if restfreq_label not in self.summary:
                        restfreq_label = 'RESTFREQ1'
                    restfreq = self.summary[restfreq_label] * u.MHz
                    data['MJD'][id0:id1] = mjd_col
                    data['RESTFREQ'][id0:id1] = restfreq.to(u.Hz).value
                    data['OBSTIME'][id0:id1] = \
                        array.meta['integration_time'].value
                    data['VELOCITY'][id0:id1] = \
                        subscan.meta['VLSR'].to("m/s").value
                    data['DATE-OBS'][id0:id1] = date_col
                    data['UT'][id0:id1] = ut_col
                    data['MH2O'][id0:id1] = mH2O

                    data['SIGNAL'][id0:id1] = on_or_off(subscan, f)
                    is_on = cal_is_on(subscan)
                    if isinstance(is_on, collections.Iterable) and average:
                        if len(list(set(is_on))) != 1:
                            raise ValueError('flag_cal is inconsistent '
                                             'in {}'.format(fname))
                        is_on = is_on[0]
                    data['CAL_IS_ON'][id0:id1] = is_on
                    label = 'SRT-{}-{}-{}'.format(subscan.meta['receiver'][0],
                                                  subscan.meta['backend'][:3],
                                                  label_from_chan_name(ch))
                    data['TELESCOP'][id0:id1] = label
                    data['TSYS'][id0:id1] = 1
                    df = (bandwidth / nbin).to('Hz')
                    data['CDELT1'][id0:id1] = df
                    data['CDELT2'][id0:id1] = \
                        subscan.meta["ra_offset"].to(u.deg).value
                    data['CDELT3'][id0:id1] = \
                        subscan.meta["dec_offset"].to(u.deg).value
                    deltav = - df / restfreq * c.c
                    data['DELTAV'][id0:id1] = deltav.to('m/s').value
                    data['LINE'][id0:id1] = \
                        "F{}-{:3.3f}-MHz".format(f, bandwidth.to('MHz').value)

                    data['OBJECT'][id0:id1] = subscan.meta['SOURCE']
                    data['AZIMUTH'][id0:id1] = azimuth
                    data['ELEVATIO'][id0:id1] = elevation
                    data['CRPIX1'][id0:id1] = nbin // 2 + 1
                    data['CRVAL2'][id0:id1] = crval2
                    data['CRVAL3'][id0:id1] = crval3
                    data['LST'][id0:id1] = lsts
                    data['MAXIS1'][id0:id1] = array.meta['channels']
                    id0 = id1

                header = newhdu[1].header
                header['CTYPE1'] = "FREQ"
                header['CRVAL'] = 0
                header['CRVAL2'] = np.mean(subscan['ra'][:, f].to(u.deg).value)
                header['CRVAL3'] = \
                    np.mean(subscan['dec'][:, f].to(u.deg).value)
                header['LINE'] = subscan.meta['SOURCE']
                header['OBJECT'] = subscan.meta['SOURCE']
                header['SOURCE'] = subscan.meta['SOURCE']
                header['DATE-RED'] = \
                    Time.now().to_datetime().strftime('%d/%m/%y')
                header['LINE'] = \
                    "FEED{}-{:3.3f}-MHz".format(f,
                                                bandwidth.to('MHz').value)
                header['CDELT1'] = df.to('Hz').value
                header['RESTFREQ'] = restfreq.to(u.Hz).value
                header['MAXIS1'] = channels[0]

                filekey = os.path.basename(scandir) + '_all_feed{}'.format(f)

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

    def calibrate_all(self, use_calon=False):
        """Calibrate the scan in all available ways.

        The basic calibration is `(on - off)/off`, where `on` and `off` are
        on-source and off-source spectra respectively.

        New HDU lists are produced and added to the existing, uncalibrated
        ones.

        If the calibration mark has been used in some scan, an additional
        calibration is applied and further HDU lists are produced.

        Other Parameters
        ----------------
        use_calon : bool, default False
            If False, only the OFF + CAL is used for the calibration. If True,
            Also the ON + CAL is used and the calibration constant is averaged
            with that obtained through OFF + CAL.

        """
        new_tables = {}
        for caltype in ["cal", "onoff"]:
            for (filekey, hdul) in self.tables.items():
                new_filekey = filekey.replace("_all", "_" + caltype)
                new_hdul = copy.deepcopy(hdul)

                table = Table(new_hdul[1].data)

                table.sort(['MJD', 'TELESCOP', 'LINE'])

                out_grouped = table.group_by(['TELESCOP', 'LINE'])
                new_rows = 0
                astropy_table_from_results = None

                apply_cal = caltype == "cal"

                for _, out_group in zip(out_grouped.groups.keys,
                                        out_grouped.groups):
                    out_group = find_cycles(out_group, ['SIGNAL', 'CAL_IS_ON'])

                    grouped = out_group.group_by(['CYCLE'])

                    for _, group in zip(grouped.groups.keys, grouped.groups):
                        # group = vstack([group, group])
                        results, _ = normalize_on_off_cal(group, smooth=False,
                                                          apply_cal=apply_cal,
                                                          use_calon=use_calon)
                        if results is None:
                            break
                        if astropy_table_from_results is None:
                            astropy_table_from_results = results
                        else:
                            astropy_table_from_results = \
                                vstack((astropy_table_from_results, results))
                        new_rows += 1
                if astropy_table_from_results is None:
                    continue

                astropy_table_from_results.remove_column('CYCLE')
                astropy_table_from_results.remove_column('SIGNAL')
                astropy_table_from_results.remove_column('CAL_IS_ON')
                astropy_table_from_results.remove_column('BINARY_COL')
                astropy_table_from_results.remove_column('MJD')

                dummy_hdu = \
                    fits.BinTableHDU(data=astropy_table_from_results)
                new_hdul[1].data = dummy_hdu.data

                new_tables[new_filekey] = new_hdul

        self.tables.update(new_tables)

    def write_tables_to_disk(self):
        """Write all HDU lists produced until now in separate FITS files."""
        for (filekey, table) in self.tables.items():
            outfile = os.path.join(self.dirname, '{}.fits'.format(filekey))
            table.writeto(outfile, overwrite=True)

            outscript = os.path.join(self.dirname,
                                     '{}_class_script.txt'.format(filekey))
            with open(outscript, 'w') as fobj:
                string = """
file out {file}.gdf multiple /overwrite
fits read {file}.fits
                """.format(file=filekey)
                print(string, file=fobj)
