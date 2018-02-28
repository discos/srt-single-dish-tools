"""Interactive operations."""
from __future__ import (absolute_import, division,
                        print_function)

import copy

try:
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

import numpy as np
from .fit import linear_fit, linear_fun, align
import warnings
from .utils import compare_anything


__all__ = ["TestWarning", "PlotWarning", "mask", "intervals", "DataSelector",
           "ImageSelector", "create_empty_info", "select_data"]


class TestWarning(UserWarning):
    pass


class PlotWarning(UserWarning):
    pass


def create_empty_info(keys):
    info = {}
    for key in keys:
        info[key] = {}
        info[key]['FLAG'] = None
        info[key]['zap'] = intervals()
        info[key]['base'] = intervals()
        info[key]['fitpars'] = np.array([0, 0])

    return info


def mask(xs, border_xs, invert=False):
    """Create mask from a list of interval borders.

    Parameters
    ----------
    xs : array
        the array of values to filter
    border_xs : array
        the list of borders. Should be an even number of positions

    Returns
    -------
    mask : array
        Array of boolean values, that work as a mask to xs

    Other Parameters
    ----------------
    invert : bool
        Mask value is False if invert = False, and vice versa.
        E.g. for zapped intervals, invert = False. For baseline fit selections,
        invert = True
    """
    good = np.ones(len(xs), dtype=bool)
    if len(border_xs) >= 2:
        intervals = list(zip(border_xs[:-1:2], border_xs[1::2]))
        for i in intervals:
            good[np.logical_and(xs >= i[0],
                                xs <= i[1])] = False
    if invert:
        good = np.logical_not(good)

    return good


class intervals():
    """A list of xs and ys of the points taken during interactive selection."""

    def __init__(self):
        """Initialize."""
        self.xs = []
        self.ys = []

    def clear(self):
        """Clear."""
        self.xs = []
        self.ys = []

    def add(self, point):
        """Append points."""
        self.xs.append(point[0])
        self.ys.append(point[1])

    def __eq__(self, other):
        if isinstance(self, other.__class__):
            return self.__dict__ == other.__dict__
        return False

    def __ne__(self, other):
        return not self.__eq__(other)


class DataSelector:
    """Plot and process scans interactively."""

    def __init__(self, xs, ys, ax1, ax2, masks=None, xlabel=None, title=None,
                 test=False):
        """Initialize."""

        self.instructions = """
-------------------------------------------------------------

Interactive plotter.

-------------------------------------------------------------

Choose line to fit: Click on the line

Interval selection: Point mouse + <key>
    z     create zap intervals
    b     suggest intervals to use for baseline fit

Flagging actions:
    x     flag as bad;
    v     Remove flags and all masks from data;

Actions:
    P     print current zap list and fit parameters
    A     align all scans w.r.t. the selected one
    u     update plots with new selections
    B     subtract the baseline;
    r     reset baseline and zapping intervals, and fit parameters;
    q     quit

-------------------------------------------------------------
    """
        if not HAS_MPL:
            raise ImportError("matplotlib not installed")
        self.xs = xs
        self.ys = ys
        self.test = test
        if masks is None:
            masks = dict(list(zip(self.xs.keys(),
                                  [np.ones(len(self.xs[k]), dtype=bool)
                                   for k in self.xs.keys()])))
        self.masks = masks
        self.ax1 = ax1
        self.ax2 = ax2

        self.xlabel = xlabel
        self.title = title
        self.starting_info = create_empty_info(self.xs.keys())
        self.info = copy.deepcopy(self.starting_info)
        self.lines = []
        if not test:
            self.print_instructions()
        self.current = None

        if not test:
            ax1.figure.canvas.mpl_connect('button_press_event', self.on_click)
            ax1.figure.canvas.mpl_connect('key_press_event', self.on_key)
            ax1.figure.canvas.mpl_connect('pick_event', self.on_pick)
            ax2.figure.canvas.mpl_connect('button_press_event', self.on_click)
            ax2.figure.canvas.mpl_connect('key_press_event', self.on_key)
            ax2.figure.canvas.mpl_connect('pick_event', self.on_pick)

        self.plot_all()
        self.zcounter = 0
        self.bcounter = 0
        if not test:
            plt.show()

    def on_click(self, event):
        """Dummy function, in case I want to do something with a click."""
        pass

    def zap(self, event):
        """Create a zap interval."""
        key = self.current
        if key is None:
            return
        self.info[key]['zap'].add([event.xdata, event.ydata])
        self.zcounter += 1
        color = 'r'
        if self.zcounter % 2 == 1:
            ls = '-'
        else:
            ls = '--'
        line = self.ax1.axvline(event.xdata, color=color, ls=ls)
        line = self.ax2.axvline(event.xdata, color=color, ls=ls)
        self.lines.append(line)
        plt.draw()
        if self.test:
            warnings.warn("I select a zap interval at {}".format(event.xdata),
                          TestWarning)

    def base(self, event):
        """Add an interval to the ones that will be used by baseline sub."""
        key = self.current
        if key is None:
            return
        self.info[key]['base'].add([event.xdata, event.ydata])
        self.bcounter += 1
        color = 'b'
        if self.bcounter % 2 == 1:
            ls = '-'
        else:
            ls = '--'
        line = self.ax1.axvline(event.xdata, color=color, ls=ls)
        line = self.ax2.axvline(event.xdata, color=color, ls=ls)
        self.lines.append(line)
        plt.draw()
        if self.test:
            warnings.warn("I put a baseline mark at {}".format(event.xdata),
                          TestWarning)

    def on_key(self, event):
        """Do something when the keyboard is used."""
        if event.key == 'z':
            self.zap(event)
        elif event.key == 'h':
            self.print_instructions()
        elif event.key == 'b':
            self.base(event)
        elif event.key == 'B':
            self.subtract_baseline()
        elif event.key == 'u':
            self.plot_all()
        elif event.key == 'x':
            self.flag()
        elif event.key == 'P':
            self.print_info()
        elif event.key == 'A':
            self.align_all()
        elif event.key == 'v':
            self.flag(value=False)
        elif event.key == 'r':
            self.reset()
        elif event.key == 'q':
            self.quit()
        else:
            pass

    def flag(self, value=True):
        self.info[self.current]['FLAG'] = value
        print('Scan was {}flagged'.format("un" if not value else ""))

    def reset(self):
        for l in self.lines:
            l.remove()
        for current in self.xs.keys():
            self.lines = []
            self.info[current]['zap'].clear()
            self.info[current]['base'].clear()
            self.info[current]['fitpars'] = np.array([0, 0])
            self.info[current]['FLAG'] = None
        self.plot_all(silent=True)

    def quit(self):
        print("Closing all figures and quitting.")
        for key in self.info.keys():
            if compare_anything(self.info[key], self.starting_info[key]):
                self.info.pop(key)
        plt.close(self.ax1.figure)

    def subtract_baseline(self):
        """Subtract the baseline based on the selected intervals."""
        key = self.current
        if len(self.info[key]['base'].xs) < 2:
            self.info[key]['fitpars'] = np.array([np.min(self.ys[key]), 0])
        else:
            base_xs = self.info[key]['base'].xs
            good = mask(self.xs[key], base_xs, invert=True)

            self.info[key]['fitpars'] = linear_fit(self.xs[key][good],
                                                   self.ys[key][good],
                                                   self.info[key]['fitpars'])

        self.plot_all(silent=True)
        if self.test:
            warnings.warn("I subtracted the baseline", TestWarning)

    def subtract_model(self, channel):
        """Subtract the model from the scan."""
        fitpars = list(self.info[channel]['fitpars'])
        return self.ys[channel] - linear_fun(self.xs[channel], *fitpars)

    def align_all(self):
        """Given the selected scan, aligns all the others to that."""

        # During tests, let's suppress all unwanted warnings
        if self.test:
            warnings.filterwarnings("ignore")
        reference = self.current

        x = np.array(self.xs[reference])
        y = np.array(self.subtract_model(reference))
        zap_xs = self.info[reference]['zap'].xs

        good = mask(x, zap_xs)

        xs = [x[good]]
        ys = [y[good]]
        keys = [reference]

        for key in self.xs.keys():
            if key == reference:
                continue

            x = np.array(self.xs[key].copy())
            y = np.array(self.ys[key].copy())

            zap_xs = self.info[key]['zap'].xs

            good = mask(x, zap_xs)

            good = good * self.masks[key]
            if len(x[good]) == 0:
                continue

            xs.append(x[good])
            ys.append(y[good])
            keys.append(key)

        # ------- Make FIT!!! -----
        qs, ms = align(xs, ys)
        # -------------------------

        for ik, key in enumerate(keys):
            if ik == 0:
                continue
            self.info[key]['fitpars'] = np.array([qs[ik - 1], ms[ik - 1]])

        self.plot_all(silent=True)
        if self.test:
            warnings.filterwarnings("default")
            warnings.warn("I aligned all", TestWarning)

    def on_pick(self, event):
        """Do this when I pick a line in the plot."""
        thisline = event.artist

        self.current = (thisline._label)
        self.plot_all(silent=True)

    def plot_all(self, silent=False):
        """Plot everything."""
        update_limits = False
        if self.lines:
            xlim_save = self.ax1.get_xlim()
            ylim_save = self.ax1.get_ylim()
            update_limits = True
        for l in self.lines:
            l.remove()
        self.lines = []
        self.ax1.cla()
        plt.setp(self.ax1.get_xticklabels(), visible=False)
        good = {}
        model = {}
        if self.current is not None:
            self.ax1.plot(self.xs[self.current], self.ys[self.current],
                          color='g', lw=3, zorder=10,
                          rasterized=True)
        for key in self.xs.keys():
            self.ax1.plot(self.xs[key], self.ys[key], color='k', picker=True,
                          label=key, lw=1, rasterized=True)

            zap_xs = self.info[key]['zap'].xs

            # Eliminate zapped intervals
            plt.draw()
            good[key] = mask(self.xs[key], zap_xs)

            if self.info[key]['FLAG'] is True:
                good[key][:] = 0
            elif self.info[key]['FLAG'] is False:
                # "v" eliminates all flags!
                good[key][:] = 1
                self.masks[key][:] = 1

            good[key] = good[key] * self.masks[key]

            fitpars = list(self.info[key]['fitpars'])

            if len(fitpars) >= 2:
                model[key] = linear_fun(self.xs[key], *fitpars)
                self.ax1.plot(self.xs[key], model[key], color='b',
                              rasterized=True)
            else:
                model[key] = np.zeros(len(self.xs[key])) + np.min(self.ys[key])

        self.ax2.cla()
        self.ax2.axhline(0, ls='--', color='k')
        for key in self.xs.keys():
            self.ax2.plot(self.xs[key], self.ys[key] - model[key],
                          color='grey', ls='-', picker=True,
                          label=key,
                          rasterized=True)
            self.ax2.plot(self.xs[key][good[key]],
                          self.ys[key][good[key]] - model[key][good[key]],
                          'k-', lw=1,
                          rasterized=True)

        if self.current is not None:
            print("Current scan is {}".format(self.current))
            key = self.current
            self.ax2.plot(self.xs[key][good[key]],
                          self.ys[key][good[key]] - model[key][good[key]],
                          color='g', lw=3, zorder=10,
                          rasterized=True)
        if self.xlabel is not None:
            self.ax2.set_xlabel(self.xlabel)

        if update_limits:
            self.ax1.set_xlim(xlim_save)
            self.ax1.set_ylim(ylim_save)
        plt.draw()
        if self.test and not silent:
            warnings.warn("I plotted all", PlotWarning)

    def print_instructions(self):
        """Print to terminal some instructions for the interactive window."""
        print(self.instructions)

    def print_info(self):
        """Print info on the current scan.

        Info includes zapped intervals and fit parameters.
        """
        for key in self.info.keys():
            print(key + ':')
            if len(self.info[key]['zap'].xs) >= 2:
                print('  Zap intervals: ',
                      list(zip(self.info[key]['zap'].xs[:-1:2],
                               self.info[key]['zap'].xs[1::2])))

            print('  Fit pars:      ', self.info[key]['fitpars'])


def select_data(xs, ys, masks=None, title=None, xlabel=None, test=False):
    """Open a DataSelector window."""
    if not HAS_MPL:
        raise ImportError("matplotlib not installed")
    try:
        xs.keys()
    except Exception:
        xs = {'Ch': xs}
        ys = {'Ch': ys}

    if title is None:
        title = 'Data selector (press "h" for help)'

    plt.figure(title)
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 2], hspace=0)

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1], sharex=ax1)

    datasel = DataSelector(xs, ys, ax1, ax2, masks=masks, title=title,
                           xlabel=xlabel, test=test)

    return datasel.info


class ImageSelector():
    """Return xs and ys of the image, and the key that was pressed.

    Attributes
    ----------
    img : array
        the image
    ax : pyplot.axis instance
        the axis where the image will be plotted
    fun : function
        the function to call when a key is pressed. It must accept three
        arguments: `x`, `y` and `key`
    """

    def __init__(self, data, ax, fun=None, test=False):
        """
        Initialize an ImageSelector class.

        Parameters
        ----------
        data : array
            the image
        ax : pyplot.axis instance
            the axis where the image will be plotted
        fun : function, optional
            the function to call when a key is pressed. It must accept three
            arguments: `x`, `y` and `key`
        """
        if not HAS_MPL:
            raise ImportError("matplotlib not installed")
        self.img = data
        self.ax = ax
        self.fun = fun
        self.plot_img()
        if not test:
            ax.figure.canvas.mpl_connect('key_press_event', self.on_key)
            plt.show()

    def on_key(self, event):
        """Do this when the keyboard is pressed."""
        x, y = event.xdata, event.ydata
        key = event.key
        print(x, y, key)

        if key == 'q':
            plt.close(self.ax.figure)
        elif x is None or y is None or x != x or y != y:
            print("Invalid choice. Is the window under focus?")
            return
        elif self.fun is not None:
            plt.close(self.ax.figure)
            self.fun(x, y, key)

        return x, y, key

    def plot_img(self):
        """Plot the image on the interactive display."""
        from .utils import ds9_like_log_scale
        img_to_plot = ds9_like_log_scale(self.img)
        self.ax.imshow(img_to_plot, origin='lower',
                       vmin=np.percentile(img_to_plot, 20),
                       interpolation="nearest", cmap="gnuplot2",
                       extent=[0, self.img.shape[1], 0, self.img.shape[0]])
