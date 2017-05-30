"""Scan class."""
from __future__ import (absolute_import, division,
                        print_function)

from .io import read_data, root_name
import glob
from .read_config import read_config, get_config_file
from .fit import ref_mad, contiguous_regions
import os
import numpy as np
from astropy.table import Table, Column
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from .fit import baseline_rough, baseline_als, linear_fun, _rolling_window
from .interactive_filter import select_data
from astropy import stats

import re
import warnings
import logging
import traceback


def _split_freq_splat(freqsplat):
    freqmin, freqmax = \
        [float(f) for f in freqsplat.split(':')]
    return freqmin, freqmax


def interpret_frequency_range(freqsplat, bandwidth, nbin):
    """Interpret the frequency range specified in freqsplat.

    Examples
    --------
    >>> interpret_frequency_range(None, 1024, 512)
    (102.4, 921.6, 51, 459)
    >>> interpret_frequency_range('default', 1024, 512)
    (102.4, 921.6, 51, 459)
    >>> interpret_frequency_range(':', 1024, 512)
    (0, 1024, 0, 511)
    >>> interpret_frequency_range('all', 1024, 512)
    (0, 1024, 0, 511)
    >>> interpret_frequency_range('200:800', 1024, 512)
    (200.0, 800.0, 100, 399)
    """

    if freqsplat is None or freqsplat == 'default':
        freqmin, freqmax = bandwidth / 10, bandwidth * 0.9
    elif freqsplat in ['all', ':']:
        freqmin, freqmax = 0, bandwidth
    else:
        freqmin, freqmax = _split_freq_splat(freqsplat)

    binmin = int(nbin * freqmin / bandwidth)
    binmax = int(nbin * freqmax / bandwidth) - 1

    return freqmin, freqmax, binmin, binmax


def _clean_dyn_spec(dynamical_spectrum, bad_intervals):
    cleaned_dynamical_spectrum = dynamical_spectrum.copy()
    for b in bad_intervals:
        if b[0] == 0:
            fill_lc = np.array(dynamical_spectrum[:, b[1]])
        elif b[1] >= dynamical_spectrum.shape[1]:
            fill_lc = np.array(dynamical_spectrum[:, b[0]])
        else:
            previous = np.array(dynamical_spectrum[:, b[0] - 1])
            next = np.array(dynamical_spectrum[:, b[1]])
            fill_lc = (previous + next) / 2

        for bsub in range(b[0], np.min([b[1], dynamical_spectrum.shape[1]])):
            cleaned_dynamical_spectrum[:, bsub] = fill_lc
        if b[0] == 0 or b[1] >= dynamical_spectrum.shape[1]:
            continue
    return cleaned_dynamical_spectrum


def _clean_scan_using_variability(dynamical_spectrum, length, bandwidth,
                                  good_mask=None, freqsplat=None,
                                  noise_threshold=5, debug=True, nofilt=False,
                                  outfile="out", label=""):
    if len(dynamical_spectrum.shape) == 1:
        return None
    dynspec_len, nbin = dynamical_spectrum.shape

    # Calculate first light curve

    times = length * np.arange(dynspec_len) / dynspec_len
    lc = np.sum(dynamical_spectrum, axis=1)
    lc = baseline_als(times, lc)
    lcbins = np.arange(len(lc))

    # Calculate spectral variability curve

    meanspec = np.sum(dynamical_spectrum, axis=0) / dynspec_len
    spectral_var = \
        np.sqrt(np.sum((dynamical_spectrum - meanspec) ** 2,
                       axis=0) / dynspec_len) / meanspec

    df = bandwidth / len(meanspec)
    allbins = np.arange(len(meanspec)) * df

    # Mask frequencies -- avoid those excluded from splat

    freqmask = np.ones(len(meanspec), dtype=bool)
    freqmin, freqmax, binmin, binmax = \
        interpret_frequency_range(freqsplat, bandwidth, nbin)
    freqmask[0:binmin] = False
    freqmask[binmax:] = False

    # Calculate the variability image

    varimg = np.sqrt((dynamical_spectrum - meanspec) ** 2) / meanspec

    # Set up corrected spectral var

    mod_spectral_var = spectral_var.copy()
    mod_spectral_var[0:binmin] = spectral_var[binmin]
    mod_spectral_var[binmax:] = spectral_var[binmax]

    # Some statistical information on spectral var

    median_spectral_var = np.median(mod_spectral_var[freqmask])
    stdref = ref_mad(mod_spectral_var[freqmask], 20)

    # Calculate baseline of spectral var ---------------
    # Empyrical formula, with no physical meaning
    lam = 10**(-6.2 + np.log2(nbin) * 1.2)
    _, baseline = baseline_als(np.arange(binmax - binmin),
                               np.array(mod_spectral_var[binmin:binmax]),
                               return_baseline=True,
                               lam=lam,
                               p=0.001, offset_correction=False,
                               outlier_purging=(False, True), niter=30)

    baseline = \
        np.concatenate((np.zeros(binmin) + baseline[0],
                        baseline,
                        np.zeros(nbin - binmax) + baseline[-1]
                        ))

    # Set threshold

    if not nofilt:
        threshold = baseline + 2 * noise_threshold * stdref
    else:
        threshold = np.zeros_like(baseline) + 1e32

    # Set mask

    mask = spectral_var < threshold
    wholemask = freqmask & mask & np.logical_not(good_mask)

    # Calculate frequency-masked lc
    lc_masked = np.sum(dynamical_spectrum[:, freqmask], axis=1)
    lc_masked = baseline_als(times, lc_masked, outlier_purging=False)

    bad_intervals = contiguous_regions(np.logical_not(wholemask))

    # Calculate cleaned dynamical spectrum

    cleaned_dynamical_spectrum = \
        _clean_dyn_spec(dynamical_spectrum, bad_intervals)

    cleaned_meanspec = \
        np.sum(cleaned_dynamical_spectrum,
               axis=0) / len(cleaned_dynamical_spectrum)
    cleaned_varimg = \
        np.sqrt((cleaned_dynamical_spectrum - cleaned_meanspec) ** 2 /
                cleaned_meanspec ** 2)
    cleaned_spectral_var = \
        np.sqrt(np.sum((cleaned_dynamical_spectrum - cleaned_meanspec) ** 2,
                       axis=0) / dynspec_len) / cleaned_meanspec

    mean_varimg = np.mean(cleaned_varimg[:, freqmask])
    std_varimg = np.std(cleaned_varimg[:, freqmask])

    lc_corr = np.sum(cleaned_dynamical_spectrum[:, freqmask], axis=1)
    lc_corr = baseline_als(times, lc_corr, outlier_purging=False)

    results = type('test', (), {})()  # create empty object
    results.lc = lc_corr
    results.freqmin = freqmin
    results.freqmax = freqmax

    if not debug:
        return results

    # Now, PLOT IT ALL --------------------------------
    # Prepare subplots
    fig = plt.figure("{}_{}".format(outfile, label), figsize=(15, 15))
    gs = GridSpec(4, 3, hspace=0, wspace=0,
                  height_ratios=(1.5, 1.5, 1.5, 1.5),
                  width_ratios=(3, 0.3, 1.2))
    ax_meanspec = plt.subplot(gs[0, 0])
    ax_dynspec = plt.subplot(gs[1, 0], sharex=ax_meanspec)
    ax_cleanspec = plt.subplot(gs[2, 0], sharex=ax_meanspec)
    ax_lc = plt.subplot(gs[1, 2], sharey=ax_dynspec)
    ax_cleanlc = plt.subplot(gs[2, 2], sharey=ax_dynspec, sharex=ax_lc)
    ax_var = plt.subplot(gs[3, 0], sharex=ax_meanspec)
    ax_varhist = plt.subplot(gs[3, 1], sharey=ax_var)
    ax_meanspec.set_ylabel('Counts')
    ax_dynspec.set_ylabel('Sample')
    ax_cleanspec.set_ylabel('Sample')
    ax_var.set_ylabel('r.m.s.')
    ax_var.set_xlabel('Frequency (MHz)')
    ax_cleanlc.set_xlabel('Counts')

    # Plot mean spectrum

    ax_meanspec.plot(allbins[1:], meanspec[1:], label="Unfiltered")
    # ax_meanspec.plot(allbins[1:], meanspec[1:], label="Whitelist applied")
    ax_meanspec.plot(allbins[wholemask], meanspec[wholemask],
                     label="Final mask")
    ax_meanspec.set_ylim([np.min(cleaned_meanspec),
                          np.max(cleaned_meanspec)])

    try:
        cmap = plt.get_cmap("magma")
    except:
        cmap = plt.get_cmap("gnuplot2")
    ax_dynspec.imshow(varimg, origin="lower", aspect='auto',
                      cmap=cmap,
                      vmin=mean_varimg - 5 * std_varimg,
                      vmax=mean_varimg + 5 * std_varimg,
                      extent=(0, bandwidth,
                              0, varimg.shape[0]), interpolation='none')

    ax_cleanspec.imshow(cleaned_varimg, origin="lower", aspect='auto',
                        cmap=cmap,
                        vmin=mean_varimg - 5 * std_varimg,
                        vmax=mean_varimg + 5 * std_varimg,
                        extent=(0, bandwidth,
                                0, varimg.shape[0]), interpolation='none')

    # Plot variability

    ax_var.plot(allbins[1:], spectral_var[1:], label="Spectral rms")
    ax_var.plot(allbins[mask], spectral_var[mask])
    ax_var.plot(allbins, cleaned_spectral_var,
                zorder=10, color="k")
    ax_var.plot(allbins[1:], baseline[1:])
    ax_var.plot(allbins[1:], baseline[1:] + 2 * noise_threshold * stdref)
    minb = np.min(baseline[1:])
    ax_var.set_ylim([minb, median_spectral_var + 10 * stdref])

    # Plot light curves

    ax_lc.plot(lc, lcbins, color="grey")
    ax_lc.plot(lc_masked, lcbins, color="b")
    ax_cleanlc.plot(lc_masked, lcbins, color="grey")
    ax_cleanlc.plot(lc_corr, lcbins, color="k")
    dlc = max(lc_corr) - min(lc_corr)
    ax_lc.set_xlim([np.min(lc_corr) - dlc / 10, max(lc_corr) + dlc / 10])

    # Indicate bad intervals

    for b in bad_intervals:
        maxsp = np.max(meanspec)
        ax_meanspec.plot(b * df, [maxsp] * 2, color='k', lw=2)
        middleimg = [varimg.shape[0] / 2]
        ax_dynspec.plot(b * df, [middleimg] * 2, color='k', lw=2)
        maxsp = np.max(spectral_var)
        ax_var.plot(b * df, [maxsp] * 2, color='k', lw=2)

    # Indicate freqmin and freqmax

    ax_dynspec.axvline(freqmin)
    ax_dynspec.axvline(freqmax)
    ax_var.axvline(freqmin)
    ax_var.axvline(freqmax)
    ax_meanspec.axvline(freqmin)
    ax_meanspec.axvline(freqmax)

    plt.savefig(
        "{}_{}.pdf".format(outfile, label))
    plt.close(fig)
    return results


chan_re = re.compile(r'^Ch[0-9]+$')


def list_scans(datadir, dirlist):
    """List all scans contained in the directory listed in config."""
    scan_list = []

    for d in dirlist:
        for f in glob.glob(os.path.join(datadir, d, '*.fits')):
            scan_list.append(f)
    return scan_list


class Scan(Table):
    """Class containing a single scan."""

    def __init__(self, data=None, config_file=None, norefilt=True,
                 interactive=False, nosave=False, debug=False,
                 freqsplat=None, nofilt=False, nosub=False, **kwargs):
        """Initialize a Scan object.

        Freqsplat is a string, freqmin:freqmax, and gives the limiting
        frequencies of the interval to splat in a single channel.
        """
        if config_file is None:
            config_file = get_config_file()

        if isinstance(data, Table):
            Table.__init__(self, data, **kwargs)
        elif data is None:
            Table.__init__(self, **kwargs)
            self.meta['config_file'] = config_file
            self.meta.update(read_config(self.meta['config_file']))
        else:  # if data is a filename
            if os.path.exists(root_name(data) + '.hdf5') and norefilt:
                data = root_name(data) + '.hdf5'
            if debug:
                logging.info('Loading file {}'.format(data))
            table = read_data(data)
            Table.__init__(self, table, masked=True, **kwargs)
            if not data.endswith('hdf5'):
                self.meta['filename'] = os.path.abspath(data)
            self.meta['config_file'] = config_file

            self.meta.update(read_config(self.meta['config_file']))

            self.check_order()

            self.clean_and_splat(freqsplat=freqsplat, nofilt=nofilt,
                                 noise_threshold=self.meta['noise_threshold'],
                                 debug=debug)

            if interactive:
                self.interactive_filter()

            if (('backsub' not in self.meta.keys() or
                    not self.meta['backsub'])) and not nosub:
                logging.info('Subtracting the baseline')
                self.baseline_subtract()

            if not nosave:
                self.save()

    def interpret_frequency_range(self, freqsplat, bandwidth, nbin):
        """Interpret the frequency range specified in freqsplat."""
        return interpret_frequency_range(freqsplat, bandwidth, nbin)

    def make_single_channel(self, freqsplat, masks=None):
        """Transform a spectrum into a single-channel count rate."""
        for ic, ch in enumerate(self.chan_columns()):
            if len(self[ch].shape) == 1:
                continue

            _, nbin = self[ch].shape

            freqmin, freqmax, binmin, binmax = \
                self.interpret_frequency_range(freqsplat,
                                               self[ch].meta['bandwidth'],
                                               nbin)

            if masks is not None:
                self[ch][:, np.logical_not(masks[ch])] = 0

            self[ch + 'TEMP'] = \
                Column(np.sum(self[ch][:, binmin:binmax], axis=1))

            self[ch + 'TEMP'].meta.update(self[ch].meta)
            self.remove_column(ch)
            self[ch + 'TEMP'].name = ch
            self[ch].meta['bandwidth'] = freqmax - freqmin

    def chan_columns(self):
        """List columns containing samples."""
        return np.array([i for i in self.columns
                         if chan_re.match(i)])

    def clean_and_splat(self, good_mask=None, freqsplat=None,
                        noise_threshold=5, debug=True,
                        save_spectrum=False, nofilt=False):
        """Clean from RFI.

        Very rough now, it will become complicated eventually.

        Parameters
        ----------
        good_mask : boolean array
            this mask specifies intervals that should never be discarded as
            RFI, for example because they contain spectral lines
        noise_threshold : float
            The threshold, in sigmas, over which a given channel is
            considered noisy
        freqsplat : str
            Specification of frequency interval to merge into a single channel

        Returns
        -------
        masks : dictionary of boolean arrays
            this dictionary contains, for each detector/polarization, True
            values for good spectral channels, and False for bad channels.

        Other parameters
        ----------------
        save_spectrum : bool, default False
            Save the spectrum into a 'ChX_spec' column
        debug : bool, default True
            Save images with quicklook information on single scans
        """
        logging.debug("Noise threshold:", noise_threshold)

        if self.meta['filtering_factor'] > 0.5:
            warnings.warn("Don't use filtering factors > 0.5. Skipping.")
            return

        chans = self.chan_columns()
        for ic, ch in enumerate(chans):
            results = \
                _clean_scan_using_variability(
                    self[ch], self['time'],
                    self[ch].meta['bandwidth'],
                    good_mask=good_mask,
                    freqsplat=freqsplat,
                    noise_threshold=noise_threshold,
                    debug=debug, nofilt=nofilt,
                    outfile=root_name(self.meta['filename']),
                    label="{}".format(ic))

            if results is None:
                continue
            lc_corr = results.lc
            freqmin, freqmax = results.freqmin, results.freqmax

            self[ch + 'TEMP'] = Column(lc_corr)

            self[ch + 'TEMP'].meta.update(self[ch].meta)
            if save_spectrum:
                self[ch].name = ch + "_spec"
            else:
                self.remove_column(ch)
            self[ch + 'TEMP'].name = ch
            self[ch].meta['bandwidth'] = freqmax - freqmin

    def baseline_subtract(self, kind='als', plot=False):
        """Subtract the baseline."""
        for ch in self.chan_columns():
            if plot:
                fig = plt.figure("Sub" + ch)
                plt.plot(self['time'], self[ch] - np.min(self[ch]),
                         alpha=0.5)

            if kind == 'als':
                self[ch] = baseline_als(self['time'], self[ch])
            elif kind == 'rough':
                self[ch] = baseline_rough(self['time'], self[ch])

            if plot:
                plt.plot(self['time'], self[ch])
                out = self.meta['filename'].replace('.fits',
                                                    '_{}.png'.format(ch))
                plt.savefig(out)
                plt.close(fig)
        self.meta['backsub'] = True

    def zap_birdies(self):
        """Zap bad intervals."""
        pass

    def __repr__(self):
        """Give the print() function something to print."""
        reprstring = \
            '\n\n----Scan from file {0} ----\n'.format(self.meta['filename'])
        reprstring += repr(Table(self))
        return reprstring

    def write(self, fname, **kwargs):
        """Set default path and call Table.write."""
        logging.info('Saving to {}'.format(fname))

        Table.write(self, fname, path='scan', serialize_meta=True, **kwargs)

    def check_order(self):
        """Check that times in a scan are monotonically increasing."""
        assert np.all(self['time'] == np.sort(self['time'])), \
            'The order of times in the table is wrong'

    def interactive_filter(self, save=True):
        """Run the interactive filter."""
        for ch in self.chan_columns():
            # Temporary, waiting for AstroPy's metadata handling improvements
            feed = self[ch + '_feed'][0]

            selection = self['ra'][:, feed]

            ravar = np.abs(selection[-1] -
                           selection[0])

            selection = self['dec'][:, feed]
            decvar = np.abs(selection[-1] -
                            selection[0])

            # Choose if plotting by R.A. or Dec.
            if ravar > decvar:
                dim = 'ra'
            else:
                dim = 'dec'

            # ------- CALL INTERACTIVE FITTER ---------
            info = select_data(self[dim][:, feed], self[ch],
                               xlabel=dim)

            # -----------------------------------------

            # Treat zapped intervals
            xs = info['Ch']['zap'].xs
            good = np.ones(len(self[dim]), dtype=bool)
            if len(xs) >= 2:
                intervals = list(zip(xs[:-1:2], xs[1::2]))
                for i in intervals:
                    good[np.logical_and(self[dim][:, feed] >= i[0],
                                        self[dim][:, feed] <= i[1])] = False
            self['{}-filt'.format(ch)] = good

            if len(info['Ch']['fitpars']) > 1:
                self[ch] -= linear_fun(self[dim][:, feed],
                                       *info['Ch']['fitpars'])
            # TODO: make it channel-independent
                self.meta['backsub'] = True

            # TODO: make it channel-independent
            if info['Ch']['FLAG']:
                self.meta['FLAG'] = True
        if save:
            self.save()
        self.meta['ifilt'] = True

    def save(self, fname=None):
        """Call self.write with a default filename, or specify it."""
        if fname is None:
            fname = root_name(self.meta['filename']) + '.hdf5'
        self.write(fname, overwrite=True)
