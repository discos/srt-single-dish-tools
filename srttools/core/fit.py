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


def minimize_align(xs, ys, params):
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
    return minimize_align(args[0], args[1], params)


def align(xs, ys):
    '''Given the first scan, it aligns all the others to that'''
    from scipy.optimize import minimize

    qs = np.zeros(len(xs) - 1)
    ms = np.zeros(len(xs) - 1)

    result = minimize(objective_function, [qs, ms], args=[xs, ys],
                      options={'disp': True})

    qs = result.x[:len(xs) - 1]
    ms = result.x[len(xs) - 1:]

    return qs, ms


def fit_baseline_plus_bell(x, y, ye=None, kind='gauss'):
    """Fit a function composed of a linear baseline plus a bell function.

    kind:     'gauss' or 'lorentz'
    """
    assert kind in ['gauss', 'lorentz'], \
        'kind has to be one of: gauss, lorentz'
    from astropy.modeling import models, fitting

    base = models.Linear1D(slope=0, intercept=np.min(y), name='Baseline')

    xrange = np.max(x) - np.min(x)
    yrange = np.max(y) - np.min(y)

    if kind == 'gauss':
        bell = models.Gaussian1D(mean=np.mean(x), stddev=xrange / 20,
                                 amplitude=yrange, name='Bell')
        bell.amplitude.bounds = (0, None)
        bell.mean.bounds = (0, None)
        bell.stddev.bounds = (0, None)
        # max_name = 'mean'
    elif kind == 'lorentz':
        bell = models.Lorentz1D(x_0=np.mean(x), fwhm=xrange / 20,
                                amplitude=yrange, name='Bell')
        bell.amplitude.bounds = (0, None)
        bell.x_0.bounds = (0, None)
        bell.fwhm.bounds = (0, None)
        # max_name = 'x_0'

    mod_init = base + bell

    fit = fitting.LevMarLSQFitter()

    mod_out = fit(mod_init, x, y)

    return mod_out, fit.fit_info


def _test_shape(x):
    return 1000 * np.exp(-(x - 50) ** 2 / 3)


def test_fit_baseline_plus_bell():
    import matplotlib.pyplot as plt

    x = np.arange(0, 100, 0.1)
    y = np.random.poisson(1000, len(x)) + _test_shape(x) + x * 6 + 20

    model, _ = fit_baseline_plus_bell(x, y, ye=10, kind='gauss')

    plt.figure('After fit')
    plt.errorbar(x, y, yerr=10)

    plt.plot(x, model(x))

    print('Amplitude: {}'.format(model.amplitude_1.value))
    print('Stddev: {}'.format(model.stddev_1.value))
    print('Mean: {}'.format(model.mean_1.value))
    plt.show()


def test_minimize_align():
    import matplotlib.pyplot as plt

    x1 = np.arange(0, 100, 0.1)
    y1 = np.random.poisson(100, len(x1)) + _test_shape(x1)
    x2 = np.arange(0.02, 100, 0.1)
    y2 = np.random.poisson(100, len(x2)) + _test_shape(x2)
    x3 = np.arange(0.053, 98.34, 0.1)
    y3 = np.random.poisson(100, len(x3)) + _test_shape(x3)

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
