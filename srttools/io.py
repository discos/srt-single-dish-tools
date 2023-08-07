"""Input/output functions."""
import astropy.io.fits as fits
from astropy.table import Table
import numpy as np
import astropy.units as u
from astropy.coordinates import (
    EarthLocation,
    AltAz,
    Angle,
    ICRS,
    GCRS,
    SkyCoord,
    get_sun,
)
import os
from astropy.time import Time
import warnings
import logging
import copy
import re
import glob
from collections.abc import Iterable
from scipy.interpolate import interp1d

from .utils import force_move_file

try:
    from sunpy.coordinates import frames, sun

    DEFAULT_SUN_FRAME = frames.Helioprojective
except ImportError:
    DEFAULT_SUN_FRAME = None


__all__ = [
    "mkdir_p",
    "detect_data_kind",
    "correct_offsets",
    "observing_angle",
    "get_rest_angle",
    "print_obs_info_fitszilla",
    "read_data_fitszilla",
    "read_data",
    "root_name",
    "get_chan_columns",
]


chan_re = re.compile(
    r"^Ch([0-9]+)$" r"|^Feed([0-9]+)_([a-zA-Z]+)$" r"|^Feed([0-9]+)_([a-zA-Z]+)_([0-9]+)$"
)


# 'srt': EarthLocation(4865182.7660, 791922.6890, 4035137.1740,
#                                   unit=u.m)
# EarthLocation(Angle("9:14:42.5764", u.deg),
#                                   Angle("39:29:34.93742", u.deg),
#                                   600 * u.meter) # not precise enough

locations = {
    "srt": EarthLocation(4865182.7660, 791922.6890, 4035137.1740, unit=u.m),
    "medicina": EarthLocation(Angle("11:38:49", u.deg), Angle("44:31:15", u.deg), 25 * u.meter),
    "greenwich": EarthLocation(lat=51.477 * u.deg, lon=0 * u.deg),
}


def interpret_chan_name(chan_name):
    """Get feed, polarization and baseband info from chan name.

    Examples
    >>> feed, polar, baseband = interpret_chan_name('blablabal')
    >>> feed  # None
    >>> polar  # None
    >>> baseband  # None
    >>> feed, polar, baseband = interpret_chan_name('Ch0')
    >>> feed
    0
    >>> polar  # None
    >>> baseband  # None
    >>> feed, polar, baseband = interpret_chan_name('Feed1_LCP')
    >>> feed
    1
    >>> polar
    'LCP'
    >>> baseband  # None
    >>> feed, polar, baseband = interpret_chan_name('Feed2_LCP_3')
    >>> feed
    2
    >>> polar
    'LCP'
    >>> baseband
    3
    """
    matchobj = chan_re.match(chan_name)
    if not matchobj:
        return None, None, None

    matches = [matchobj.group(i) for i in range(7)]
    polar, baseband = None, None
    if matches[6] is not None:
        baseband = int(matchobj.group(6))
        polar = matchobj.group(5)
        feed = int(matchobj.group(4))
    elif matches[3] is not None:
        polar = matchobj.group(3)
        feed = int(matchobj.group(2))
    else:
        feed = int(matchobj.group(1))

    return feed, polar, baseband


def classify_chan_columns(chans):
    """Classify the name of channels per feed, polarization, baseband.

    Examples
    --------
    >>> chans = ['Feed0_LCP_3', 'Feed0_RCP_3']
    >>> classif = classify_chan_columns(chans)
    >>> classif[0][3]['LCP']
    'Feed0_LCP_3'
    >>> classif[0][3]['RCP']
    'Feed0_RCP_3'
    >>> chans = ['Ch0']
    >>> classif = classify_chan_columns(chans)
    >>> classif[0][1]['N']
    'Ch0'
    >>> chans = ['Feed0_LCP']
    >>> classif = classify_chan_columns(chans)
    >>> classif[0][1]['LCP']
    'Feed0_LCP'
    """
    combinations = {}
    for ch in chans:
        feed, polar, baseband = interpret_chan_name(ch)
        if baseband is None:
            baseband = 1
        if polar is None:
            polar = "N"
        if feed not in combinations:
            combinations[feed] = {}

        if baseband not in combinations[feed]:
            combinations[feed][baseband] = {}

        combinations[feed][baseband][polar] = ch

    return combinations


def get_chan_columns(table):
    return np.array([i for i in table.columns if chan_re.match(i)])


def get_channel_feed(ch):
    if re.search("Feed?", ch):
        return int(ch[4])


def mkdir_p(path):
    """Safe mkdir function.

    Parameters
    ----------
    path : str
        Name of the directory/ies to create

    Notes
    -----
    Found at
    https://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
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
    if np.any(np.abs(derot_angle) > 2 * 360):
        return False
    return True


def detect_data_kind(fname):
    """Placeholder for function that recognizes data format."""
    if fname.endswith(".hdf5"):
        return "hdf5"
    elif "fits" in fname:
        return "fitszilla"
    else:
        warnings.warn("File {} is not in a known format".format(fname))
        return None


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
    sep = np.sqrt(xoffset**2.0 + yoffset**2.0)

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
    if not hasattr(rest_angle, "unit"):
        rest_angle *= u.rad
    if not hasattr(derot_angle, "unit"):
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
        return np.array([0] * len(xoffsets))
    xoffsets = np.asarray(xoffsets)
    yoffsets = np.asarray(yoffsets)
    n_lat_feeds = len(xoffsets) - 1
    rest_angle_default = _rest_angle_default(n_lat_feeds) * 2 * np.pi * u.rad
    w_0 = np.where((xoffsets[1:] > 0) & (yoffsets[1:] == 0.0))[0][0]
    return np.concatenate(([0], np.roll(rest_angle_default.to(u.rad).value, w_0))) * u.rad


def infer_skydip_from_elevation(elevation, azimuth=None):
    if azimuth is None:
        azimuth = np.array([0, 0])

    el_condition = np.max(elevation) - np.min(elevation) > np.pi / 3.0
    az_condition = np.max(azimuth) - np.min(azimuth) < 0.1 / 180.0 * np.pi
    return az_condition & el_condition


def get_sun_coords_from_radec(obstimes, ra, dec, sun_frame=None):
    if sun_frame is None:  # pragma: no cover
        sun_frame = DEFAULT_SUN_FRAME

    coords = GCRS(
        ra=Angle(ra),
        dec=Angle(dec),
        obstime=obstimes,
        distance=sun.earth_distance(obstimes),
    )

    coords_asec = coords.transform_to(sun_frame(obstime=obstimes, observer="earth"))

    lon = coords_asec.Tx
    lat = coords_asec.Ty
    dist = coords_asec.distance

    return lon.to(u.radian), lat.to(u.radian), dist.to(u.m).value


def update_table_with_sun_coords(new_table, feeds=None, sun_frame=None):
    lon_str, lat_str = "hpln", "hplt"

    if not ("dsun" in new_table.colnames):
        new_table[lon_str] = np.zeros_like(new_table["el"])
        new_table[lat_str] = np.zeros_like(new_table["az"])
        new_table["dsun"] = np.zeros(len(new_table["az"]))

    if feeds is None:
        feeds = np.arange(0, new_table["el"].shape[1], dtype=int)

    for i in feeds:
        obstimes = Time(new_table["time"] * u.day, format="mjd", scale="utc")

        lon, lat, dist = get_sun_coords_from_radec(
            obstimes,
            new_table["ra"][:, i],
            new_table["dec"][:, i],
            sun_frame=sun_frame,
        )
        new_table[lon_str][:, i] = lon
        new_table[lat_str][:, i] = lat
        if i == 0:
            new_table["dsun"][:] = dist

    return new_table


def get_coords_from_altaz_offset(obstimes, el, az, xoffs, yoffs, location, inplace=False):
    """"""
    # Calculate observing angle
    if not inplace:
        el = copy.deepcopy(el)
        az = copy.deepcopy(az)

    el += yoffs.to(u.rad).value
    az += xoffs.to(u.rad).value / np.cos(el)

    coords = AltAz(az=Angle(az), alt=Angle(el), location=location, obstime=obstimes)

    # According to line_profiler, coords.icrs is *by far* the longest
    # operation in this function, taking between 80 and 90% of the
    # execution time. Need to study a way to avoid this.
    coords_deg = coords.transform_to(ICRS())
    ra = np.radians(coords_deg.ra)
    dec = np.radians(coords_deg.dec)
    return ra, dec


def is_close_to_sun(ra, dec, obstime, tolerance=3 * u.deg):
    """Test if current source is close to the Sun.

    Examples
    --------
    >>> ra, dec = 131.13535699 * u.deg, 18.08202663 * u.deg
    >>> obstime = Time("2017-08-01")
    >>> is_close_to_sun(ra, dec, obstime, tolerance=3 * u.deg)
    True
    >>> is_close_to_sun(ra, dec + 4 * u.deg, obstime, tolerance=3 * u.deg)
    False
    """
    coords = SkyCoord(ra=ra, dec=dec, frame=GCRS(obstime=obstime))
    sun_position = get_sun(obstime).transform_to(GCRS(obstime=obstime))

    return (coords.separation(sun_position)).to(u.deg).value < tolerance.value


def update_table_with_offsets(
    new_table, xoffsets, yoffsets, rest_angles, feeds=None, inplace=False
):
    if not inplace:
        new_table = copy.deepcopy(new_table)

    lon_str, lat_str = "ra", "dec"

    if not (lon_str in new_table.colnames):
        new_table[lon_str] = np.zeros_like(new_table["el"])
        new_table[lat_str] = np.zeros_like(new_table["az"])

    if feeds is None:
        feeds = np.arange(0, new_table["el"].shape[1], dtype=int)

    for i in feeds:
        obs_angle = observing_angle(rest_angles[i], new_table["derot_angle"])

        # offsets < 0.001 arcseconds: don't correct (usually feed 0)
        if (
            np.abs(xoffsets[i]) < np.radians(0.001 / 60.0) * u.rad
            and np.abs(yoffsets[i]) < np.radians(0.001 / 60.0) * u.rad
        ):
            continue
        xoffs, yoffs = correct_offsets(obs_angle, xoffsets[i], yoffsets[i])
        obstimes = Time(new_table["time"] * u.day, format="mjd", scale="utc")

        location = locations[new_table.meta["site"]]
        lon, lat = get_coords_from_altaz_offset(
            obstimes,
            new_table["el"][:, i],
            new_table["az"][:, i],
            xoffs,
            yoffs,
            location=location,
            inplace=inplace,
        )
        new_table[lon_str][:, i] = lon
        new_table[lat_str][:, i] = lat

    return new_table


def print_obs_info_fitszilla(fname):
    """Placeholder for function that prints out oberving information."""
    with fits.open(fname, memmap=False) as lchdulist:
        section_table_data = lchdulist["SECTION TABLE"].data
        sample_rates = get_value_with_units(section_table_data, "sampleRate")

        print("Sample rates:", sample_rates)

        rf_input_data = lchdulist["RF INPUTS"].data
        print("Feeds          :", get_value_with_units(rf_input_data, "feed"))
        print("IFs            :", get_value_with_units(rf_input_data, "ifChain"))
        print(
            "Polarizations  :",
            get_value_with_units(rf_input_data, "polarization"),
        )
        print(
            "Frequencies    :",
            get_value_with_units(rf_input_data, "frequency"),
        )
        print(
            "Bandwidths     :",
            get_value_with_units(rf_input_data, "bandWidth"),
        )


def _chan_name(f, p, c=None):
    if c is not None:
        return "Feed{}_{}_{}".format(f, p, c)
    else:
        return "Feed{}_{}".format(f, p)


def read_data_fitszilla(fname):
    with fits.open(fname, memmap=False) as lchdulist:
        retval = _read_data_fitszilla(lchdulist)
    return retval


def get_value_with_units(fitsext, keyword, default=""):
    if isinstance(fitsext, fits.BinTableHDU):
        fitsext = fitsext.data
    unitstr = fitsext.columns[keyword].unit
    if unitstr is None:
        if default not in ["", None]:
            unit = u.Unit(default)
        else:
            unit = 1
    else:
        unit = u.Unit(unitstr)
    value = fitsext[keyword]
    is_string = isinstance(value, str)
    is_iterable = isinstance(value, Iterable)
    if is_string or (is_iterable and isinstance(value[0], str)):
        return value
    else:
        return value * unit


def adjust_temperature_size_rough(temp, comparison_array):
    """Adjust the size of the temperature array.

    Examples
    --------
    >>> temp = [1, 2, 3, 4]
    >>> adjust_temperature_size_rough(temp, [5, 6, 7])
    array([1, 2, 3])
    >>> adjust_temperature_size_rough(temp, [5, 6, 7, 5, 4])
    array([1, 2, 3, 4, 4])
    >>> adjust_temperature_size_rough(temp, [5, 6])
    array([2, 3])
    >>> adjust_temperature_size_rough(temp, [5, 6, 7, 5, 4, 6])
    array([1, 1, 2, 3, 4, 4])
    """
    import copy

    temp = np.asarray(temp)
    comparison_array = np.asarray(comparison_array)

    temp_save = copy.deepcopy(temp)

    sizediff = temp.size - comparison_array.size
    if sizediff > 0:
        temp = temp[sizediff // 2 : sizediff // 2 + comparison_array.size]
    elif sizediff < 0:
        # make it positive
        sizediff = -sizediff
        temp = np.zeros_like(comparison_array)
        temp[sizediff // 2 : sizediff // 2 + temp_save.size] = temp_save
        temp[: sizediff // 2] = temp_save[0]
        temp[sizediff // 2 + temp_save.size - 1 :] = temp_save[-1]

    return temp


def adjust_temperature_size(temp, comparison_array):
    """Adjust the size of the temperature array.

    Examples
    --------
    >>> temp = [1, 2, 3, 4]
    >>> np.allclose(adjust_temperature_size(temp, [5, 6]), [1.0, 4.0])
    True
    >>> temp = [1, 2, 3, 4]
    >>> np.allclose(adjust_temperature_size(temp, [5, 6, 4, 5]), temp)
    True
    """
    temp = np.asarray(temp)
    comparison_array = np.asarray(comparison_array)

    Ntemp = temp.shape[0]
    Ndata = comparison_array.shape[0]
    if Ntemp == Ndata:
        return temp

    temp_func = interp1d(np.linspace(0, 1, Ntemp), temp)

    newtemp = temp_func(np.linspace(0, 1, Ndata))
    return newtemp


# from memory_profiler import profile
# @profile
def _read_data_fitszilla(lchdulist):
    """Open a fitszilla FITS file and read all relevant information."""
    is_new_fitszilla = np.any(["coord" in i.name.lower() for i in lchdulist])

    # ----------- Extract generic observation information ------------------
    headerdict = dict(lchdulist[0].header.items())
    source = lchdulist[0].header["SOURCE"]
    site = lchdulist[0].header["ANTENNA"].lower()
    receiver = lchdulist[0].header["RECEIVER CODE"]

    ra = lchdulist[0].header["RIGHTASCENSION"] * u.rad
    dec = lchdulist[0].header["DECLINATION"] * u.rad
    ra_offset = dec_offset = az_offset = el_offset = 0 * u.rad
    if "RightAscension Offset" in lchdulist[0].header:
        ra_offset = lchdulist[0].header["RightAscension Offset"] * u.rad
    if "Declination Offset" in lchdulist[0].header:
        dec_offset = lchdulist[0].header["Declination Offset"] * u.rad
    if "Azimuth Offset" in lchdulist[0].header:
        az_offset = lchdulist[0].header["Azimuth Offset"] * u.rad
    if "Elevation Offset" in lchdulist[0].header:
        el_offset = lchdulist[0].header["Elevation Offset"] * u.rad

    # ----------- Read the list of channel ids ------------------
    section_table_data = lchdulist["SECTION TABLE"].data
    chan_ids = get_value_with_units(section_table_data, "id")
    nbin_per_chan = get_value_with_units(section_table_data, "bins")
    sample_rate = get_value_with_units(section_table_data, "sampleRate")
    try:
        bw_section = get_value_with_units(section_table_data, "bandWidth")
        fr_section = get_value_with_units(section_table_data, "frequency")
    except KeyError:
        bw_section = None
        fr_section = None
    integration_time = lchdulist["SECTION TABLE"].header["Integration"] * u.ms
    if len(list(set(nbin_per_chan))) > 1:
        raise ValueError(
            "Only datasets with the same nbin per channel are " "supported at the moment"
        )
    nbin_per_chan = list(set(nbin_per_chan))[0]
    types = get_value_with_units(section_table_data, "type")
    if "stokes" in types:
        is_polarized = True
    else:
        is_polarized = False

    # Check. If backend is not specified, use Total Power
    try:
        backend = lchdulist[0].header["BACKEND NAME"]
    except Exception:
        if "stokes" in types:
            if nbin_per_chan == 2048:
                backend = "XARCOS"
            else:
                backend = "SARDARA"
        elif "spectra" in types:
            backend = "SARDARA"
        else:
            backend = "TP"

    # ----------- Read the list of RF inputs, feeds, polarization, etc. --
    rf_input_data = lchdulist["RF INPUTS"].data
    feeds = get_value_with_units(rf_input_data, "feed")
    IFs = get_value_with_units(rf_input_data, "ifChain")
    polarizations = get_value_with_units(rf_input_data, "polarization")
    sections = get_value_with_units(rf_input_data, "section")
    frequencies_rf = get_value_with_units(rf_input_data, "frequency")
    bandwidths_rf = get_value_with_units(rf_input_data, "bandWidth")
    local_oscillator = get_value_with_units(rf_input_data, "localOscillator")
    try:
        cal_mark_temp = get_value_with_units(rf_input_data, "calibrationMark")
    except KeyError:
        # Old, stupid typo
        cal_mark_temp = get_value_with_units(rf_input_data, "calibratonMark")

    if bw_section is not None:
        bandwidths_section = [bw_section[i] for i in sections]
        frequencies_section = [fr_section[i] for i in sections]
        frequencies_section = [f + l for (f, l) in zip(frequencies_section, local_oscillator)]

    if backend == "TP" or bw_section is None:
        frequencies, bandwidths = frequencies_rf, bandwidths_rf
    else:
        frequencies, bandwidths = frequencies_section, bandwidths_section

    combinations = list(zip(frequencies, bandwidths))
    combination_idx = np.arange(len(combinations))

    # Solve stupid problem with old CCB data
    if receiver.lower() == "ccb":
        feeds[:] = 0

    if len(set(combinations)) > 1:
        chan_names = [_chan_name(f, p, c) for f, p, c in zip(feeds, polarizations, combination_idx)]
    else:
        chan_names = [_chan_name(f, p) for f, p in zip(feeds, polarizations)]

    # ----- Read the offsets of different feeds (nonzero only if multifeed)--
    feed_input_data = lchdulist["FEED TABLE"].data
    # Add management of historical offsets.
    # Note that we need to add the units by hand in this case.
    xoffsets = get_value_with_units(feed_input_data, "xOffset", default="rad")
    yoffsets = get_value_with_units(feed_input_data, "yOffset", default="rad")

    relpowers = get_value_with_units(feed_input_data, "relativePower")

    # -------------- Read data!-----------------------------------------
    datahdu = lchdulist["DATA TABLE"]
    # N.B.: there is an increase in memory usage here. This is just because
    # data are being read from the file at this point, not before.
    data_table_data = Table(datahdu.data)
    tempdata = Table(lchdulist["ANTENNA TEMP TABLE"].data)

    for col in data_table_data.colnames:
        if col == col.lower():
            continue
        data_table_data.rename_column(col, col.lower())
    for col in tempdata.colnames:
        if col == col.lower():
            continue
        tempdata.rename_column(col, col.lower())

    is_old_spectrum = "SPECTRUM" in list(datahdu.header.values())
    if is_old_spectrum:
        data_table_data.rename_column("spectrum", "ch0")
        sections = np.array([0, 0])

    unsupported_temperature = False
    if len(tempdata[tempdata.colnames[0]].shape) == 2:
        try:
            tempdata_new = Table()
            for i, (feed, ifnum) in enumerate(zip(feeds, IFs)):
                tempdata_new[f"ch{i}"] = tempdata[f"ch{feed}"][:, ifnum]
            tempdata = tempdata_new
        except Exception:  # pragma: no cover
            warnings.warn("Temperature format not supported", UserWarning)
            unsupported_temperature = True
            pass

    existing_columns = [chn for chn in data_table_data.colnames if chn.startswith("ch")]
    if existing_columns == []:
        raise ValueError("Invalid data")

    is_spectrum = nbin_per_chan > 1

    is_single_channel = len(set(combinations)) == 1

    good = np.ones(len(feeds), dtype=bool)

    for i, s in enumerate(sections):
        section_name = "ch{}".format(s)
        if section_name not in existing_columns:
            good[i] = False
    allfeeds = feeds
    feeds = allfeeds[good]
    IFs = IFs[good]
    polarizations = polarizations[good]
    sections = sections[good]

    unique_feeds = np.unique(feeds)

    rest_angles = get_rest_angle(xoffsets, yoffsets)

    if is_spectrum:
        nchan = len(chan_ids)

        sample_channel = existing_columns[0]

        _, nbins = data_table_data[sample_channel].shape

        # Development version of SARDARA -- will it remain the same?
        if nbin_per_chan == nbins:
            IFs = np.zeros_like(IFs)

        if nbin_per_chan * nchan * 2 == nbins and not is_polarized:
            warnings.warn(
                "Data appear to contain polarization information "
                "but are classified as simple, not stokes, in the "
                "Section table."
            )
            is_polarized = True

        if (
            nbin_per_chan != nbins
            and nbin_per_chan * nchan != nbins
            and nbin_per_chan * nchan * 2 != nbins
            and not is_polarized
        ):
            raise ValueError(
                "Something wrong with channel subdivision: "
                "{} bins/channel, {} channels, "
                "{} total bins".format(nbin_per_chan, nchan, nbins)
            )

        for f, ic, p, s in zip(feeds, IFs, polarizations, sections):
            c = s
            if is_single_channel:
                c = None
            section_name = "ch{}".format(s)
            ch = _chan_name(f, p, c)
            start, end = ic * nbin_per_chan, (ic + 1) * nbin_per_chan
            data_table_data[ch] = data_table_data[section_name][:, start:end]

        if is_polarized:
            # for f, ic, p, s in zip(feeds, IFs, polarizations, sections):
            for s in list(set(sections)):
                f = feeds[sections == s][0]
                c = s
                if is_single_channel:
                    c = None

                section_name = "ch{}".format(s)
                qname, uname = _chan_name(f, "Q", c), _chan_name(f, "U", c)
                qstart, qend = 2 * nbin_per_chan, 3 * nbin_per_chan
                ustart, uend = 3 * nbin_per_chan, 4 * nbin_per_chan
                data_table_data[qname] = data_table_data[section_name][:, qstart:qend]
                data_table_data[uname] = data_table_data[section_name][:, ustart:uend]

                chan_names += [qname, uname]

        for f, ic, p, s in zip(feeds, IFs, polarizations, sections):
            section_name = "ch{}".format(s)
            if section_name in data_table_data.colnames:
                data_table_data.remove_column(section_name)
    else:
        for ic, ch in enumerate(chan_names):
            data_table_data[ch] = data_table_data["ch{}".format(chan_ids[ic])]

    # ----------- Read temperature data, if possible ----------------
    for ic, ch in enumerate(chan_names):
        data_table_data[ch + "-Temp"] = 0.0
        if unsupported_temperature:
            continue

        if len(chan_ids) <= ic:
            continue
        ch_string = f"ch{chan_ids[ic]}"
        if ch_string not in tempdata.colnames:
            continue

        td = np.asarray(tempdata[ch_string])
        data_table_data[ch + "-Temp"] = adjust_temperature_size(td, data_table_data[ch + "-Temp"])

    info_to_retrieve = [
        "time",
        "derot_angle",
        "weather",
        "par_angle",
        "flag_track",
        "flag_cal",
    ] + [ch + "-Temp" for ch in chan_names]

    new_table = Table()

    new_table.meta.update(headerdict)
    new_table.meta["SOURCE"] = source
    new_table.meta["site"] = site
    new_table.meta["backend"] = backend
    new_table.meta["receiver"] = receiver
    new_table.meta["RA"] = ra
    new_table.meta["Dec"] = dec
    new_table.meta["channels"] = nbin_per_chan
    new_table.meta["VLSR"] = new_table.meta["VLSR"] * u.Unit("km/s")

    for i, off in zip(
        "ra,dec,el,az".split(","),
        [ra_offset, dec_offset, el_offset, az_offset],
    ):
        new_table.meta[i + "_offset"] = off

    for info in info_to_retrieve:
        new_table[info] = data_table_data[info]

    if not _check_derotator(new_table["derot_angle"]):
        logging.debug("Derotator angle looks weird. Setting to 0")
        new_table["derot_angle"][:] = 0

    # Duplicate raj and decj columns (in order to be corrected later)
    Nfeeds = np.max(allfeeds) + 1

    for newcol, oldcol in [("ra", "raj2000"), ("dec", "decj2000"), ("el", "el"), ("az", "az")]:
        new_table[newcol] = np.tile(data_table_data[oldcol], (Nfeeds, 1)).transpose()

    new_table.meta["is_skydip"] = infer_skydip_from_elevation(
        data_table_data["el"], data_table_data["az"]
    )

    for info in ["ra", "dec", "az", "el", "derot_angle"]:
        new_table[info].unit = u.radian

    if not is_new_fitszilla:
        update_table_with_offsets(
            new_table, xoffsets, yoffsets, rest_angles, feeds=unique_feeds, inplace=True
        )
    else:
        for i in range(len(xoffsets)):
            try:
                ext = lchdulist["Coord{}".format(i)]
                extdata = ext.data
                ra, dec = extdata["raj2000"], extdata["decj2000"]
                el, az = extdata["el"], extdata["az"]
            except KeyError:
                ra, dec = new_table["ra"][:, 0], new_table["dec"][:, 0]
                el, az = new_table["el"][:, 0], new_table["az"][:, 0]

            new_table["ra"][:, i] = ra
            new_table["dec"][:, i] = dec
            new_table["el"][:, i] = el
            new_table["az"][:, i] = az

    # Don't know if better euristics is needed
    obstime = Time(np.mean(new_table["time"]) * u.day, format="mjd", scale="utc")
    if is_close_to_sun(
        new_table.meta["RA"],
        new_table.meta["Dec"],
        obstime,
        tolerance=3 * u.deg,
    ):
        if DEFAULT_SUN_FRAME is None:
            raise ValueError("You need Sunpy to process Sun observations.")
        update_table_with_sun_coords(
            new_table,
            feeds=unique_feeds,
            sun_frame=DEFAULT_SUN_FRAME,
        )

    lchdulist.close()

    # So ugly. But it works
    filtered_frequencies = [f for (f, g) in zip(frequencies, good) if g]

    for i, fr in enumerate(filtered_frequencies):
        f = feeds[i]
        s = sections[i]
        ic = IFs[i]
        p = polarizations[i]
        b = bandwidths[i]
        lo = local_oscillator[i]
        cal = cal_mark_temp[i]

        c = s
        if is_single_channel:
            c = None
        chan_name = _chan_name(f, p, c)
        if bandwidths[ic] < 0:
            frequencies[ic] -= bandwidths[ic]
            bandwidths[ic] *= -1
            for i in range(data_table_data[chan_name].shape[0]):
                data_table_data[chan_name][f, :] = data_table_data[chan_name][f, ::-1]

        new_table[chan_name] = data_table_data[chan_name] * relpowers[feeds[ic]]

        new_table[chan_name + "-filt"] = np.ones(len(data_table_data[chan_name]), dtype=bool)
        data_table_data.remove_column(chan_name)

        newmeta = {
            "polarization": polarizations[ic],
            "feed": int(f),
            "IF": int(ic),
            "frequency": fr.to("MHz"),
            "bandwidth": b.to("MHz"),
            "sample_rate": sample_rate[s],
            "sample_time": (1 / (sample_rate[s].to(u.Hz))).to("s"),
            "local_oscillator": lo.to("MHz"),
            "cal_mark_temp": cal.to("K"),
            "integration_time": integration_time.to("s"),
            "xoffset": xoffsets[f].to(u.rad),
            "yoffset": yoffsets[f].to(u.rad),
            "relpower": float(relpowers[f]),
        }
        new_table[chan_name].meta.update(headerdict)
        new_table[chan_name].meta.update(new_table.meta)
        new_table[chan_name].meta.update(newmeta)

    if is_polarized:
        for s in list(set(sections)):
            feed = feeds[sections == s][0]
            c = s
            if is_single_channel:
                c = None
            for stokes_par in "QU":
                chan_name = _chan_name(feed, stokes_par, c)
                try:
                    new_table[chan_name] = data_table_data[chan_name]
                except KeyError:
                    continue
                sample_time = 1 / (sample_rate[s].to(u.Hz))

                newmeta = {
                    "polarization": stokes_par,
                    "feed": int(feed),
                    "IF": -1,
                    # There are two IFs for each section
                    "frequency": frequencies[2 * s].to("MHz"),
                    "bandwidth": bandwidths[2 * s].to("MHz"),
                    "sample_rate": sample_rate[s],
                    "sample_time": sample_time.to("s"),
                    "local_oscillator": local_oscillator[2 * s].to("MHz"),
                    "cal_mark_temp": cal_mark_temp[2 * s].to("K"),
                    "integration_time": integration_time.to("s"),
                    "xoffset": xoffsets[feed].to(u.rad),
                    "yoffset": yoffsets[feed].to(u.rad),
                    "relpower": 1.0,
                }
                new_table[chan_name].meta.update(headerdict)
                new_table[chan_name].meta.update(new_table.meta)
                new_table[chan_name].meta.update(newmeta)

                new_table[chan_name + "-filt"] = np.ones(
                    len(data_table_data[chan_name]), dtype=bool
                )
                data_table_data.remove_column(chan_name)

    return new_table


def read_data(fname):
    """Read the data, whatever the format, and return them."""
    kind = detect_data_kind(fname)
    if kind == "fitszilla":
        return read_data_fitszilla(fname)
    elif kind == "hdf5":
        return Table.read(fname)
    else:
        return None


def root_name(fname):
    """Return the file name without extension."""
    fn, ext = os.path.splitext(fname)
    if "fits" in ext and not ext.endswith("fits"):
        fn += ext.replace("fits", "").replace(".", "")
    return fn


def _try_type(value, dtype):
    """
    Examples
    --------
    >>> _try_type("1", int)
    1
    >>> _try_type(1.0, int)
    1
    >>> _try_type("ab", float)
    'ab'
    """
    try:
        return dtype(value)
    except ValueError:
        return value


def label_from_chan_name(ch):
    """
    Examples
    --------
    >>> label_from_chan_name('Feed0_LCP_1')
    'LL'
    >>> label_from_chan_name('Feed0_Q_2')
    'LR'
    >>> label_from_chan_name('Feed3_RCP_1')
    'RR'
    >>> label_from_chan_name('Feed2_U_3')
    'RL'
    """
    _, polar, _ = interpret_chan_name(ch)

    if polar.startswith("L"):
        return "LL"
    elif polar.startswith("R"):
        return "RR"
    elif polar.startswith("Q"):
        return "LR"
    elif polar.startswith("U"):
        return "RL"
    else:
        raise ValueError("Unrecognized polarization")


def bulk_change(file, path, value):
    """Bulk change keyword or column values in FITS file.

    Parameters
    ----------
    file : str
        Input file
    path : str
        it has to be formatted as EXT,data,COLUMN or EXT,header,KEY depending
        on what is being changed (a data column or a header key resp.). Ex.
        1,TIME to change the values of column TIME in ext. n. 1
    value : any acceptable type
        Value to be filled in
    """

    with fits.open(file, memmap=False) as hdul:
        ext, attr, key = path.split(",")
        ext = _try_type(ext, int)

        data = getattr(hdul[ext], attr)
        data[key] = value
        setattr(hdul[ext], attr, data)

        hdul.writeto("tmp.fits", overwrite=True)
    force_move_file("tmp.fits", file)


def main_bulk_change(args=None):
    """Preprocess the data."""
    import argparse

    description = "Change all values of a given column or header keyword in " "fits files"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "files",
        nargs="*",
        help="Single files to preprocess",
        default=None,
        type=str,
    )
    parser.add_argument(
        "-k",
        "--key",
        type=str,
        default=None,
        help="Path to key or data column. E.g. "
        '"EXT,header,KEY" to change key KEY in the header'
        "in extension EXT; EXT,data,COL to change column"
        "COL in the data of extension EXT",
    )
    parser.add_argument("-v", "--value", default=None, type=str, help="Value to be written")
    parser.add_argument(
        "--apply-cal-mark",
        action="store_true",
        default=False,
        help='Short for -k "DATA TABLE,data,flag_cal" -v 1',
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        default=False,
        help="Look for file in up to two subdirectories",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="Plot stuff and be verbose",
    )

    args = parser.parse_args(args)

    if args.apply_cal_mark:
        args.key = "DATA TABLE,data,flag_cal"
        args.value = 1

    if args.key is None:
        raise ValueError(
            "What should I do? Please specify either key and " "value, or apply-cal-mark"
        )

    fnames = []
    for fname in args.files:
        if args.recursive:
            if not fname == os.path.basename(fname):
                raise ValueError(
                    "Options recursive requires a file name, not " "a full path: {}".format(fname)
                )

            fs = glob.glob(os.path.join("**", fname), recursive=True)

            fnames.extend(fs)
        else:
            fnames.append(fname)

    for fname in fnames:
        print("Updating", fname, "...", end="")

        bulk_change(fname, args.key, args.value)
        print(fname, " Done.")
