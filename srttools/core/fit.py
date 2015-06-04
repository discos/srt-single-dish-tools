from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
from scipy.optimize import curve_fit
import numpy as np


def linear_fun(x, q, m):
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
        start_pars = [q0, m0]

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


def function_to_minimize_align(xs, ys, params):
    '''Calculate the total variance of a series of scans, after subtracting
    a linear function from each of them (excluding the first one)'''

    params = np.array(params).flatten()
    qs = params[:len(xs) - 1]
    ms = params[len(xs) - 1:]

    x = xs[0].copy()
    y = ys[0].copy()
    for i in range(1, len(xs)):
        x = np.append(x, xs[i])
        scaled_y = ys[i] - (xs[i] * ms[i - 1] + qs[i - 1])
        y = np.append(y, scaled_y)

    order = np.argsort(x)

    x = x[order]
    y = y[order]

    x_range = [np.min(x), np.max(x)]

    xints = np.linspace(x_range[0], x_range[1], len(x) / 20)

    values = np.array([np.var(y[(x >= xints[k]) & (x < xints[k+1])])
                      for k in range(len(xints[:-1]))])
    good = values == values
    value = np.mean(values[good])

    return value


def objective_function(params, args):
    '''Put the parameters in the right order to use with scipy's minimize'''
    return function_to_minimize_align(args[0], args[1], params)


def align(xs, ys):
    '''Given the first scan, it aligns all the others to that'''
    from scipy.optimize import minimize

    qs = np.zeros(len(xs) - 1)
    ms = np.zeros(len(xs) - 1)

    result = minimize(objective_function, [qs, ms], args=[xs, ys],
                      options={'disp': True})

    print(result)

    qs = result.x[:len(xs) - 1]
    ms = result.x[len(xs) - 1:]

    return qs, ms


def test_function_to_minimize_align():
    import matplotlib.pyplot as plt

    shape = lambda x: 100 * np.exp(-(x - 50) ** 2 / 3)
    x1 = np.arange(0, 100, 0.1)
    y1 = np.random.poisson(100, len(x1)) + shape(x1)
    x2 = np.arange(0.02, 100, 0.1)
    y2 = np.random.poisson(100, len(x2)) + shape(x2)
    x3 = np.arange(0.053, 98.34, 0.1)
    y3 = np.random.poisson(100, len(x3)) + shape(x3)

    xs = [x1, x2, x3]
    ys = [y1, y2, y3]

    qs = [20, -60, 60]
    ms = [0, 0.3, -0.8]

    plt.figure('Original')
    for ix, x in enumerate(xs):
        ys[ix] = ys[ix] + qs[ix] + ms[ix] * xs[ix]

        plt.plot(xs[ix], ys[ix], label=ix)

    plt.legend()

    qs, ms = align(xs, ys)

    plt.figure('After fit')
    for ix, x in enumerate(xs):
        if ix == 0:
            plt.plot(xs[0], ys[0], label=ix)
        else:
            ys[ix] = ys[ix] - (qs[ix-1] + ms[ix-1] * xs[ix])

            plt.plot(xs[ix], ys[ix], label=ix)

    plt.legend()
    plt.show()
