"""Interactive operations."""
from __future__ import (absolute_import, division,
                        print_function)

from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np
from .fit import linear_fit, linear_fun, align


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


class DataSelector:
    """Plot and process scans interactively."""

    def __init__(self, xs, ys, ax1, ax2, masks=None, xlabel=None, title=None):
        """Initialize."""
        self.xs = xs
        self.ys = ys
        if masks is None:
            masks = dict(list(zip(self.xs.keys(),
                                  [np.ones(len(self.xs[k]), dtype=bool)
                                   for k in self.xs.keys()])))
        self.masks = masks
        self.ax1 = ax1
        self.ax2 = ax2

        self.xlabel = xlabel
        self.title = title
        self.info = {}
        for key in self.xs.keys():
            self.info[key] = {}
            self.info[key]['FLAG'] = False
            self.info[key]['zap'] = intervals()
            self.info[key]['base'] = intervals()
            self.info[key]['fitpars'] = np.array([0, 0])
        self.lines = []
        self.print_instructions()
        self.current = None
        ax1.figure.canvas.mpl_connect('button_press_event', self.on_click)
        ax1.figure.canvas.mpl_connect('key_press_event', self.on_key)
        ax1.figure.canvas.mpl_connect('pick_event', self.on_pick)
        ax2.figure.canvas.mpl_connect('button_press_event', self.on_click)
        ax2.figure.canvas.mpl_connect('key_press_event', self.on_key)
        ax2.figure.canvas.mpl_connect('pick_event', self.on_pick)
        self.plot_all()
        self.zcounter = 0
        self.bcounter = 0
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

    def on_key(self, event):
        """Do something when the keyboard is used."""
        current = self.current
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
            self.info[current]['FLAG'] = True
            print('Marked as flagged')
        elif event.key == 'P':
            self.print_info()
        elif event.key == 'A':
            self.align_all()
        elif event.key == 'v':
            if self.info[current]['FLAG']:
                self.info[current]['FLAG'] = False
                print('Removed flag ()')
        elif event.key == 'r':
            for l in self.lines:
                l.remove()
            for current in self.xs.keys():
                self.lines = []
                self.info[current]['zap'].clear()
                self.info[current]['base'].clear()
                self.info[current]['fitpars'] = np.array([0, 0])
            self.plot_all()
        elif event.key == 'q':
            plt.close(self.ax1.figure)
        else:
            pass

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

        self.plot_all()

    def subtract_model(self, channel):
        """Subtract the model from the scan."""
        fitpars = list(self.info[channel]['fitpars'])
        return self.ys[channel] - linear_fun(self.xs[channel], *fitpars)

    def align_all(self):
        """Given the selected scan, aligns all the others to that."""
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

        self.plot_all()

    def on_pick(self, event):
        """Do this when I pick a line in the plot."""
        thisline = event.artist

        self.current = (thisline._label)
        self.plot_all()

    def plot_all(self):
        """Plot everything."""
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
                          label=key,
                          rasterized=True)

            zap_xs = self.info[key]['zap'].xs

            # Eliminate zapped intervals
            plt.draw()
            good[key] = mask(self.xs[key], zap_xs)
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
                          'k-', lw=2,
                          rasterized=True)

        if self.current is not None:
            key = self.current
            self.ax2.plot(self.xs[key][good[key]],
                          self.ys[key][good[key]] - model[key][good[key]],
                          color='g', lw=3, zorder=10,
                          rasterized=True)
        if self.xlabel is not None:
            self.ax2.set_xlabel(self.xlabel)
        plt.draw()

    def print_instructions(self):
        """Print to terminal some instructions for the interactive window."""
        instructions = """
-------------------------------------------------------------

Interactive plotter.

-------------------------------------------------------------

Choose line to fit: Click on the line

Interval selection: Point mouse + <key>
    z     create zap intervals
    b     suggest intervals to use for baseline fit

Flagging actions:
    x     flag as bad;
    v     Remove flag;

Actions:
    P     print current zap list and fit parameters
    A     align all scans w.r.t. the selected one
    u     update plots with new selections
    B     subtract the baseline;
    r     reset baseline and zapping intervals, and fit parameters;
    q     quit

-------------------------------------------------------------
    """

        print(instructions)

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


def select_data(xs, ys, masks=None, title=None, xlabel=None):
    """Open a DataSelector window."""
    try:
        xs.keys()
    except Exception:
        xs = {'Ch': xs}
        ys = {'Ch': ys}

    if title is None:
        title = 'Data selector'

    plt.figure(title)
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 2], hspace=0)

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1], sharex=ax1)

    datasel = DataSelector(xs, ys, ax1, ax2, masks=masks, title=title,
                           xlabel=xlabel)

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

    def __init__(self, data, ax, fun=None):
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
        self.img = data
        self.ax = ax
        self.fun = fun
        self.plot_img()
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
        self.ax.imshow(np.log10(self.img), origin='lower',
                       vmin=np.percentile(np.log10(self.img), 20),
                       interpolation="nearest", cmap="gnuplot2",
                       extent=[0, self.img.shape[1], 0, self.img.shape[0]])
