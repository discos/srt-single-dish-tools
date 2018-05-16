from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit
from .utils import HAS_MPL
if HAS_MPL:
    import matplotlib.pyplot as plt


def exptau(airmass, tatm, tau, t0):
    """Function to fit to the T vs airmass data."""
    bx = np.exp(-tau * airmass)
    return tatm * (1 - bx) + t0


def calculate_opacity(file, plot=True, tatm=None, tau0=None, t0=None):
    """Calculate opacity from a skydip scan.

    Atmosphere temperature is fixed, from Buffa et al.'s calculations.

    Parameters
    ----------
    file : str
        File name of the skydip scan in Fits format
    plot : bool
        Plot diagnostics about the fit

    Other parameters
    ----------------
    tatm : float
        Atmospheric temperature (fixed in the fit). The default value is
        calculated from an empyrical formula
    tau0 : float
        Initial opacity in the fit. The default value is
        np.log(2 / (1 + np.sqrt(1 - 4 * (t30 - t90) / tatm))), where
        t30 and t90 are the Tsys values calculated at 30 and 90 degrees of
        elevation respectively.
    t0 : float
        Initial value for Tsys in the fit.

    Returns
    -------
    opacities : dict
        Dictionary containing the opacities calculated for each channel, plus
        the time in the middle of the observation.
    """
    with fits.open(file) as hdulist:
        data = hdulist['DATA TABLE'].data
        tempdata = hdulist['ANTENNA TEMP TABLE'].data
        rfdata = hdulist['RF INPUTS'].data

    time = np.mean(data['Time'])

    freq = (rfdata['frequency'] + rfdata['bandwidth'] / 2)[0]

    elevation = data['el']
    airmass = 1 / np.sin(elevation)
    if tatm is None:
        airtemp = np.median(data['weather'][:, 1])
        tatm = 0.683 * (airtemp + 273.15) + 78

    el30 = np.argmin(np.abs(elevation - np.radians(30)))
    el90 = np.argmin(np.abs(elevation - np.radians(90)))

    results = {'time': time}
    for ch in ['Ch0', 'Ch1']:
        temp = tempdata[ch]
        if plot and HAS_MPL:
            fig = plt.figure(ch)
            plt.scatter(airmass, temp, c='k')

        if tau0 is None:
            t90 = temp[el90]
            t30 = temp[el30]
            tau0 = np.log(2 / (1 + np.sqrt(1 - 4 * (t30 - t90) / tatm)))

        if t0 is None:
            t0 = freq / 1e3

        init_par = [tatm, tau0, t0]

        epsilon = 1.e-5
        par, _ = curve_fit(exptau, airmass, temp, p0=init_par,
                           maxfev=10000000,
                           bounds=([tatm - epsilon, -np.inf, -np.inf],
                                   [tatm + epsilon, np.inf, np.inf]))

        print('The opacity for channel {} is {}'.format(ch, par[1]))
        if plot and HAS_MPL:
            plt.plot(airmass, exptau(airmass, *par), color='r', zorder=10)
            plt.xlabel('Airmass')
            plt.ylabel('T (K)')
            plt.title('T_atm: {:.2f}; tau: {:.4f}; t0: {:.2f}'.format(*par))
            plt.savefig(file.replace('.fits', '_fit_{}.png'.format(ch)))
            plt.close(fig)

        results[ch] = par[1]

    return results


def main_opacity(args=None):
    import argparse

    description = ('Calculate opacity from a skydip scan and plot the fit '
                   'results')
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("files", nargs='+',
                        help="File to inspect",
                        default=None, type=str)

    parser.add_argument("--tatm", type=float, default=None,
                        help='Atmospheric temperature')

    parser.add_argument("--tau0", type=float, default=None,
                        help='Initial value for tau (to be fit)')

    parser.add_argument("--t0", type=float, default=None,
                        help='Initial value for Tsys (to be fitted)')

    args = parser.parse_args(args)

    for f in args.files:
        _ = calculate_opacity(f, tatm=args.tatm, tau0=args.tau0, t0=args.t0)
