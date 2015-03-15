from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
from scipy.optimize import curve_fit
import numpy as np


def linear_fun(x, m, q):
    '''A linear function'''
    return m * x + q


def linear_fit(time, lc, start_pars, return_err=False):
    '''A linear fit with any set of data. Return the parameters'''
    par, pcov = curve_fit(linear_fun, time, lc, start_pars,
                          maxfev=6000)
    if return_err:
        pass
    else:
        return par


def rough_baseline_sub(time, lc, start_pars=None):
    '''Rough function to subtract the baseline'''

    if start_pars is None:
        m0 = 0
        q0 = min(lc)
        start_pars = [m0, q0]

    # only consider start and end quarters of image
    nbin = len(time)
    #    bins = np.arange(nbin, dtype=int)

    for percentage in [0.8, 0.15]:
        sorted_els = np.argsort(lc)
        # Select the lowest half elements
        good = sorted_els[: int(nbin * percentage)]
        #    good = np.logical_or(bins <= nbin / 4, bins >= nbin / 4 * 3)

        time_filt = time[good]
        lc_filt = lc[good]
        back_in_order = np.argsort(time_filt)
        lc_filt = lc_filt[back_in_order]
        time_filt = time_filt[back_in_order]

        par = linear_fit(time_filt, lc_filt, start_pars)

        lc -= linear_fun(time, *par)

    return lc


