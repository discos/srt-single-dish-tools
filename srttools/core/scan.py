"""Scan class."""
from __future__ import (absolute_import, division,
                        print_function)

from .io import read_data, root_name
import glob
from .read_config import read_config, get_config_file
from .fit import ref_std
import os
import numpy as np
from astropy.table import Table, Column
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from .fit import baseline_rough, baseline_als, linear_fun, _rolling_window
from .interactive_filter import select_data

import re
import warnings
import logging
import traceback


def contiguous_regions(condition):
    """Find contiguous True regions of the boolean array "condition".

    Return a 2D array where the first column is the start index of the region
    and the second column is the end index.

    Parameters
    ----------
    condition : boolean array

    Returns
    -------
    idx : [[i0_0, i0_1], [i1_0, i1_1], ...]
        A list of integer couples, with the start and end of each True blocks
        in the original array

    Notes
    -----
    From http://stackoverflow.com/questions/4494404/find-large-number-of-consecutive-values-fulfilling-condition-in-a-numpy-array
    """  # NOQA
    # Find the indicies of changes in "condition"
    diff = np.diff(condition)
    idx, = diff.nonzero()
    # We need to start things after the change in "condition". Therefore,
    # we'll shift the index by 1 to the right.
    idx += 1
    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = np.r_[0, idx]
    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = np.r_[idx, condition.size]
    # Reshape the result into two columns
    idx.shape = (-1, 2)
    return idx


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

    def __init__(self, data=None, config_file=None, norefilt=False,
                 interactive=False, nosave=False, verbose=True,
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
            if os.path.exists(root_name(data) + '.hdf5'):
                data = root_name(data) + '.hdf5'
            if verbose:
                logging.info('Loading file {}'.format(data))
            table = read_data(data)
            Table.__init__(self, table, masked=True, **kwargs)
            self.meta['filename'] = os.path.abspath(data)

            self.meta['config_file'] = config_file

            self.meta.update(read_config(self.meta['config_file']))

            self.check_order()

            self.clean_and_splat(freqsplat=freqsplat, nofilt=nofilt,
                                 noise_threshold=self.meta['noise_threshold'])

            if interactive:
                self.interactive_filter()

            if (('backsub' not in self.meta.keys() or
                    not self.meta['backsub']) \
                    or not norefilt) and not nosub:
                logging.info('Subtracting the baseline')
                self.baseline_subtract()

            if not nosave:
                self.save()

    def interpret_frequency_range(self, freqsplat, bandwidth, nbin):
        """Interpret the frequency range specified in freqsplat."""
        try:
            freqmin, freqmax = \
                [float(f) for f in freqsplat.split(':')]
        except:
            freqsplat = ":"

        if freqsplat == ":" or freqsplat == "all" or freqsplat is None:
            freqmin = 0
            freqmax = bandwidth

        binmin = int(nbin * freqmin / bandwidth)
        binmax = int(nbin * freqmax / bandwidth) - 1

        return freqmin, freqmax, binmin, binmax

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

    def clean_and_splat(self, good_mask=None, freqsplat=None, noise_threshold=5,
                        debug=True,
                        save_spectrum=False, nofilt=False):
        """Clean from RFI.

        Very rough now, it will become complicated eventually.

        Parameters
        ----------
        good_mask : boolean array
            this mask specifies intervals that should never be discarded as
            RFI, for example because they contain spectral lines
        noise_threshold : float
            The threshold, in sigmas, over which a given channel is considered noisy
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
            if len(self[ch].shape) == 1:
                break
            _, nbin = self[ch].shape

            lc = np.sum(self[ch], axis=1)
            lc = baseline_als(self['time'], lc)
            lcbins = np.arange(len(lc))
            total_spec = np.sum(self[ch], axis=0) / len(self[ch])
            spectral_var = \
                np.sqrt(np.sum((self[ch] - total_spec) ** 2 / total_spec ** 2,
                        axis=0))

            df = self[ch].meta['bandwidth'] / len(total_spec)
            allbins = np.arange(len(total_spec)) * df

            freqmask = np.ones(len(total_spec), dtype=bool)

            freqmin, freqmax, binmin, binmax = \
                self.interpret_frequency_range(freqsplat,
                                               self[ch].meta['bandwidth'],
                                               nbin)
            freqmask[0:binmin] = False
            freqmask[binmax:] = False

            if debug:
                fig = plt.figure("{}_{}".format(self.meta['filename'], ic))
                gs = GridSpec(3, 2, hspace=0, height_ratios=(1.5, 3, 1.5),
                              width_ratios=(3, 1.5))
                ax1 = plt.subplot(gs[0, 0])
                ax2 = plt.subplot(gs[1, 0], sharex=ax1)
                ax3 = plt.subplot(gs[1, 1], sharey=ax2)
                ax4 = plt.subplot(gs[2, 0], sharex=ax1)
                ax1.plot(allbins, total_spec, label="Unfiltered")
                ax4.plot(allbins, spectral_var, label="Spectral rms")

            if good_mask is not None:
                total_spec[good_mask] = 0

            varimg = np.sqrt((self[ch] - total_spec) ** 2 / total_spec ** 2)
            mean_varimg = np.mean(varimg[:, freqmask])
            std_varimg = np.std(varimg[:, freqmask])

            stdref = ref_std(spectral_var[freqmask], np.max([nbin // 20, 20]))

            mod_spectral_var = spectral_var.copy()
            mod_spectral_var[0:binmin] = spectral_var[binmin]
            mod_spectral_var[binmax:] = spectral_var[binmax]


            _, baseline = baseline_als(np.arange(len(spectral_var)),
                                       mod_spectral_var, return_baseline=True,
                                       lam=1000, p=0.001, offset_correction=False,
                                       outlier_purging=False)
            if not nofilt:
                threshold = baseline + noise_threshold * stdref
            else:
                threshold = np.zeros_like(baseline) + 1e32

            mask = spectral_var < threshold

            wholemask = freqmask & mask
            lc_corr = np.sum(self[ch][:, wholemask], axis=1)
            lc_corr = baseline_als(self['time'], lc_corr)

            if debug:
                ax1.plot(allbins, total_spec, label="Whitelist applied")
                ax1.axvline(binmin)
                ax1.axvline(binmax)
                ax1.plot(allbins[mask], total_spec[mask],
                         label="Final mask")
                # ax1.legend()

                ax2.imshow(varimg, origin="lower", aspect='auto',
                           cmap=plt.get_cmap("magma"),
                           vmin=mean_varimg - 5 * std_varimg,
                           vmax=mean_varimg + 5 * std_varimg,
                           extent=(0, self[ch].meta['bandwidth'],
                                   0, varimg.shape[0]))

                bad_intervals = contiguous_regions(np.logical_not(wholemask))
                for b in bad_intervals:
                    maxsp = np.max(total_spec)
                    ax1.plot(b*df, [maxsp]*2, color='k')
                    middleimg = [varimg.shape[0] / 2]
                    ax2.plot(b*df, [middleimg]*2, color='k')
                    maxsp = np.max(spectral_var)
                    ax4.plot(b*df, [maxsp]*2, color='k')


                ax2.axvline(binmin)
                ax2.axvline(binmax)

                ax3.plot(lc, lcbins)
                ax3.plot(lc_corr, lcbins)
                ax3.set_xlim([np.min(lc), max(lc)])
                ax4.axvline(freqmin)
                ax4.axvline(freqmax)
                ax4.plot(allbins[mask], spectral_var[mask])
                ax4.plot(allbins, baseline)

                plt.savefig(
                    "{}_{}.pdf".format(
                        root_name(self.meta['filename']), ic))
                plt.close(fig)

            self[ch + 'TEMP'] = Column(lc_corr)

            self[ch + 'TEMP'].meta.update(self[ch].meta)
            if save_spectrum:
                self[ch].name = ch + "_spec"
            else:
                self.remove_column(ch)
            self[ch + 'TEMP'].name = ch
            self[ch].meta['bandwidth'] = freqmax - freqmin

    def baseline_subtract(self, kind='als'):
        """Subtract the baseline."""
        if kind == 'als':
            for col in self.chan_columns():
                self[col] = baseline_als(self['time'], self[col])
        elif kind == 'rough':
            for col in self.chan_columns():
                self[col] = baseline_rough(self['time'], self[col])

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
        t = Table(self)
        t.write(fname, path='scan', **kwargs)

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
