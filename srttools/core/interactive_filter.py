from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np
from .fit import linear_fit, linear_fun


class intervals():
    '''A list of xs and ys of the points taken during interactive selection'''
    def __init__(self):
        self.xs = []
        self.ys = []

    def clear(self):
        self.xs = []
        self.ys = []

    def add(self, point):
        self.xs.append(point[0])
        self.ys.append(point[1])


class DataSelector:
    def __init__(self, xs, ys, ax1, ax2, xlabel=None, title=None):
        self.xs = xs
        self.ys = ys
        self.ax1 = ax1
        self.ax2 = ax2
        self.xlabel = xlabel
        self.title = title
        self.info = {}
        self.info['FLAG'] = False
        self.info['zap'] = intervals()
        self.info['base'] = intervals()
        self.info['fitpars'] = []
        self.lines = []
        self.print_instructions()
        self.plot_all()

        ax1.figure.canvas.mpl_connect('button_press_event', self.on_click)
        ax1.figure.canvas.mpl_connect('key_press_event', self.on_key)
        self.zcounter = 0
        self.bcounter = 0

    def on_click(self, event):
        pass

    def zap(self, event):
        '''Create a zap interval'''
        self.info['zap'].add([event.xdata, event.ydata])
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
        '''Adds an interval to the ones that will be used by baseline sub.'''
        self.info['base'].add([event.xdata, event.ydata])
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
        '''What to do when the keyboard is used'''

        if event.key == 'z':
            self.zap(event)
        if event.key == 'h':
            self.print_instructions()
        if event.key == 'b':
            self.base(event)
        if event.key == 'B':
            self.subtract_baseline()
        if event.key == 'u':
            self.plot_all()
        if event.key == 'x':
            self.info['FLAG'] = True
            print('Marked as flagged')
        if event.key == 'v':
            if self.info['FLAG']:
                self.info['FLAG'] = False
                print('Removed flag ()')
        elif event.key == 'r':
            for l in self.lines:
                l.remove()
            self.lines = []
            self.info['zap'].clear()
            self.info['base'].clear()
            self.info['fitpars'] = []
            self.plot_all()
        elif event.key == 'q':
            plt.close(self.ax1.figure)
        else:
            pass

    def subtract_baseline(self):
        '''Subtract the baseline based on the selected intervals'''
        if len(self.info['base'].xs) < 2:
            self.info['fitpars'] = [np.min(self.ys)]
        else:
            base_xs = self.info['base'].xs
            good = np.zeros(len(self.xs), dtype=bool)
            intervals = list(zip(base_xs[:-1:2], base_xs[1::2]))
            for i in intervals:
                good[np.logical_and(self.xs >= i[0],
                                    self.xs <= i[1])] = True
            self.info['fitpars'] = linear_fit(self.xs[good], self.ys[good],
                                              [np.min(self.ys), 0])

        self.plot_all()

    def plot_all(self):
        for l in self.lines:
            l.remove()
        self.lines = []
        self.ax1.cla()
        self.ax1.plot(self.xs, self.ys, color='k')

        plt.setp(self.ax1.get_xticklabels(), visible=False)

        zap_xs = self.info['zap'].xs

        # Eliminate zapped intervals
        plt.draw()
        good = np.ones(len(self.xs), dtype=bool)
        if len(zap_xs) >= 2:
            intervals = list(zip(zap_xs[:-1:2], zap_xs[1::2]))
            for i in intervals:
                good[np.logical_and(self.xs >= i[0],
                                    self.xs <= i[1])] = False
        fitpars = self.info['fitpars']
        if len(self.info['fitpars']) >= 2:
            model = linear_fun(self.xs, *fitpars)
            self.ax1.plot(self.xs, model, color='b')
        else:
            model = np.zeros(len(self.xs)) + np.min(self.ys)

        self.ax2.cla()
        self.ax2.axhline(0, ls='--', color='k')
        self.ax2.plot(self.xs, self.ys - model, color='grey', ls='-')
        self.ax2.plot(self.xs[good], self.ys[good] - model[good], 'k-', lw=2)

        if self.xlabel is not None:
            self.ax2.set_xlabel(self.xlabel)
        plt.draw()

    def print_instructions(self):
            instructions = '''
-------------------------------------------------------------

Interactive plotter.

-------------------------------------------------------------

Interval selection: Point mouse + <key>
    z     create zap intervals
    b     suggest intervals to use for baseline fit

Flagging actions:
    x     flag as bad;
    v     Remove flag;

Actions:
    u     update plots with new selections
    B     subtract the baseline;
    r     reset baseline and zapping intervals, and fit parameters;
    q     quit

-------------------------------------------------------------
    '''

            print(instructions)


def select_data(xs, ys, title=None, xlabel=None):

    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 2], hspace=0)

    if title is not None:
        plt.title(title)

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1], sharex=ax1)

    datasel = DataSelector(xs, ys, ax1, ax2, title=title, xlabel=xlabel)

    plt.show()

    return datasel.info


class ImageSelector():
    '''Return xs and ys of the image, and the key that was pressed.
    Inputs:
        data:       the image
        ax:         a pyplot.axis instance where the image will be plotted
        fun:        (optional) the function to call when a key is pressed.
                    it must accept three arguments: `x`, `y` and `key`
    '''
    def __init__(self, data, ax, fun=None):
        self.img = data
        self.ax = ax
        self.fun = fun
        self.plot_img()
        ax.figure.canvas.mpl_connect('key_press_event', self.on_key)

    def on_key(self, event):
        x, y = event.xdata, event.ydata
        key = event.key

        if x is None or y is None or x != x or y != y:
            print("Invalid choice. Is the window under focus?")
            return
        if self.fun is not None:
            self.fun(x, y, key)
        else:
            print(x, y, key)
        return x, y, key

    def plot_img(self):
        self.ax.imshow(self.img)


def test_select_data():
    times = np.arange(0, 100, 0.1)
    lc = np.random.normal(0, 0.3, len(times)) + times * 0.5 +\
        10 * np.exp(-(times - 30) ** 2 / 20 ** 2) + \
        3 * np.exp(-(times - 55) ** 2 / 6 ** 2) + \
        10 * np.exp(-(times - 80) ** 2 / 0.1 ** 2)
    select_data(times, lc)


def test_imageselector():
    img = np.random.normal(100, 100, (100, 100))
    fig = plt.figure()
    ax = fig.add_subplot(111)

    imagesel = ImageSelector(img, ax)

    plt.show()


