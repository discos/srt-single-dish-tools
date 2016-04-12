"""Functions to clean up images by fitting linear trends to the initial scans."""

from srttools.core.scan import contiguous_regions

from numba import jit, vectorize
from .histograms import histogram2d
import numpy as np


@vectorize ('(float64(float64,float64,float64,float64))', nopython=True)
def _align_fast(x, scan, m, q):
    """Align ``scan`` to a linear function."""
    return scan - x * m - q


XBUFFER = None
YBUFFER = None


def _get_coords(xedges, yedges):
    """Get coordinates given the edges of the histogram."""
    global XBUFFER, YBUFFER
    if XBUFFER is None:
        xcenters = (xedges[:-1] + xedges[1:]) / 2
        ycenters = (yedges[:-1] + yedges[1:]) / 2

        X, Y = np.meshgrid(xcenters, ycenters)
        XBUFFER = X
        YBUFFER = Y
    return XBUFFER, YBUFFER


EXPOMAP=None


def _calculate_image(x, y, counts, bx, by, nsamp):
    """Calculate the image."""
    global EXPOMAP

    if EXPOMAP is None:
        EXPOMAP, xedges, yedges = histogram2d(x, y, bins=(bx, by),
                                              weights=nsamp)

    histograms, xedges, yedges = histogram2d(x, y, bins=(bx, by),
                                             weights=[counts * nsamp, (counts) ** 2 * nsamp])

    img, img_var = histograms
    X, Y = _get_coords(xedges, yedges)

    good = EXPOMAP > 0
    mean = img.copy()

    mean[good] /= EXPOMAP[good]

    img_var[good] = img_var[good] / EXPOMAP[good] - mean[good] ** 2

    return X, Y, mean.T, img_var.T


@jit(nopython=True)
def _align_all(newd_t, newd_c, data_idx, par):
    ms = np.zeros_like(newd_c, dtype=np.float64)
    qs = np.zeros_like(newd_c, dtype=np.float64)

    for i_p in range(0, len(par), 2):
        i0, i1 = data_idx[i_p // 2]
        if i0 == i1:
            continue
        sliceobj = slice(i0, i1)
        ms[sliceobj] = par[i_p]
        qs[sliceobj] = par[i_p + 1]

    return _align_fast(newd_t, newd_c, ms, qs)


def _obj_fun(par, data, data_idx, excluded, bx, by):
    """
    This is the function we have to minimize.

    Parameters
    ----------
    par : array([m0, q0, m1, q1, ...])
        linear baseline parameters for the image.
    data : [times, idxs, x, y, counts]
        All five quantities are ``numpy`` ``array``s; ``time`` is time
        from the start of the scan; ``x``, ``y`` are the image coordinates,
        ``idx`` corresponds to the scan number and `counts` to the scan
        values at those coordinates.
    excluded : [[centerx0, centery0, radius0]]
        list of circular regions to exclude from fitting (e.g. strong sources
        that might alter the total rms)

    """

    newd_t, newd_i, newd_x, newd_y, newd_c, newd_e = data

    newd_c_new = _align_all(newd_t, newd_c, data_idx, par)
    X, Y, img, img_var = _calculate_image(newd_x, newd_y, newd_c_new, bx, by, newd_e)

    good = np.ones_like(img, dtype=bool)
    if excluded is not None:
        for e in excluded:
            centerx, centery, radius = e
            filt = (X - centerx) ** 2 + (Y - centery) ** 2 < radius ** 2
            good[filt] = 0

    stat = np.sum(img_var[good]) + np.var(img[good]) * img[good].size
    return stat


def _resample_scans(data):
    """Resample all scans to match the pixels of the image."""
    t, idx, x, y, c = data

    xmax, xmin = np.max(x), np.min(x)
    ymax, ymin = np.max(y), np.min(y)

    x_range = xmax - xmin
    y_range = ymax - ymin

    bx = np.linspace(xmin, xmax, int(x_range) + 1)
    by = np.linspace(ymin, ymax, int(y_range) + 1)

    newt = np.array([], dtype=np.float64)
    newi = np.array([], dtype=int)
    newx = np.array([], dtype=np.float64)
    newy = np.array([], dtype=np.float64)
    newc = np.array([], dtype=np.float64)
    newe = np.array([], dtype=np.float64)
    for i in list(set(idx)):
        good = idx == i
        x_filt = x[good]
        n = len(x_filt)
        if n == 0:
            continue
        y_filt = y[good]
        c_filt = c[good]
        t_filt = t[good]
        t_filt -= t_filt[0]

        hists, xbins, ybins = \
            histogram2d(x_filt, y_filt, bins=(bx, by),
                        weights=[np.ones(n), t_filt, x_filt, y_filt, c_filt])
        expo, time, X, Y, counts = hists
        good = expo > 0

        goodexpo = expo[good]
        tdum = np.ndarray.flatten(time[good] / goodexpo)
        cdum = np.ndarray.flatten(counts[good] / goodexpo)
        idum = np.ndarray.flatten(i + np.zeros(len(goodexpo), dtype=int))
        xdum = np.ndarray.flatten(X[good] / goodexpo)
        ydum = np.ndarray.flatten(Y[good] / goodexpo)
        edum = np.ndarray.flatten(goodexpo)

        newt = np.append(newt, tdum)
        newc = np.append(newc, cdum)
        newi = np.append(newi, idum)
        newx = np.append(newx, xdum)
        newy = np.append(newy, ydum)
        newe = np.append(newe, edum)

    return [newt, newi, newx, newy, newc, newe], bx, by


def _get_data_idx(par, idx):
    """Get the index in the data arrays corresponding to different scans."""
    data_idx = []

    for i_p, p in enumerate(list(zip(par[:-1:2], par[1::2]))):
        good = idx == i_p
        if not np.any(good):
            continue
        data_idx.append(contiguous_regions(good)[0])

    data_idx = np.array(data_idx, dtype=int)
    return data_idx


def fit_full_image(scanset, chan=0, excluded=None, par=None):
    """Get a clean image by subtracting linear trends from the initial scans.

    Parameters
    ----------
    scanset : a :class:``ScanSet`` instance
        The scanset to be fit

    Other parameters
    ----------------
    chan : int
        channel of the scanset to be fit. Defaults to 0
    excluded : [[centerx0, centery0, radius0]]
        List of circular regions to exclude from fitting (e.g. strong sources
        that might alter the total rms)
    par : [m0, q0, m1, q1, ...] or None
        Initial parameters -- slope and intercept for linear trends to be
        subtracted from the scans

    """
    from scipy.optimize import minimize
    global EXPOMAP, XBUFFER, YBUFFER
    EXPOMAP = None
    XBUFFER = None
    YBUFFER = None

    X = np.array(scanset['x'][:, 0])
    Y = np.array(scanset['y'][:, 0])

    times = np.array(scanset['time'], dtype=np.float64)
    times -= times[0]

    data_to_fit = [times,
                   np.array(scanset['Scan_id'], dtype=int),
                   np.array(X, dtype=np.float64),
                   np.array(Y, dtype=np.float64),
                   np.array(scanset['Ch0'], dtype=np.float64)]

    data, bx, by = _resample_scans(data_to_fit)

    t, i, x, y, c, expo = data

    data_idx = _get_data_idx(par, i)

    res = minimize(_obj_fun, par, args=(data, data_idx, excluded, bx, by),
                   method="SLSQP")

    return res