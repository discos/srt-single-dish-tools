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
        for key in self.xs.keys():
            self.info[key] = {}
            self.info[key]['FLAG'] = False
            self.info[key]['zap'] = intervals()
            self.info[key]['base'] = intervals()
            self.info[key]['fitpars'] = []
        self.lines = []
        self.print_instructions()
        self.current = None
        self.plot_all()

        ax1.figure.canvas.mpl_connect('button_press_event', self.on_click)
        ax1.figure.canvas.mpl_connect('key_press_event', self.on_key)
        ax1.figure.canvas.mpl_connect('pick_event', self.on_pick)
        self.zcounter = 0
        self.bcounter = 0

    def on_click(self, event):
        pass

    def zap(self, event):
        '''Create a zap interval'''
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
        '''Adds an interval to the ones that will be used by baseline sub.'''
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
        '''What to do when the keyboard is used'''

        current = self.current
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
            self.info[current]['FLAG'] = True
            print('Marked as flagged')
        if event.key == 'P':
            self.print_info()
        if event.key == 'v':
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
                self.info[current]['fitpars'] = []
            self.plot_all()
        elif event.key == 'q':
            plt.close(self.ax1.figure)
        else:
            pass

    def subtract_baseline(self):
        '''Subtract the baseline based on the selected intervals'''
        key = self.current
        if len(self.info[key]['base'].xs) < 2:
            self.info[key]['fitpars'] = [np.min(self.ys[key])]
        else:
            base_xs = self.info[key]['base'].xs
            good = np.zeros(len(self.xs[key]), dtype=bool)
            intervals = list(zip(base_xs[:-1:2], base_xs[1::2]))
            for i in intervals:
                good[np.logical_and(self.xs[key] >= i[0],
                                    self.xs[key] <= i[1])] = True
            self.info[key]['fitpars'] = linear_fit(self.xs[key][good],
                                              self.ys[key][good],
                                              [np.min(self.ys[key]), 0])

        self.plot_all()

    def on_pick(self, event):
        thisline = event.artist
        self.current=(thisline._label)
        self.plot_all()

    def plot_all(self):
        for l in self.lines:
            l.remove()
        self.lines = []
        self.ax1.cla()
        plt.setp(self.ax1.get_xticklabels(), visible=False)
        good = {}
        model = {}

        if self.current is not None:
            self.ax1.plot(self.xs[self.current], self.ys[self.current], color='g', lw=3)
        for key in self.xs.keys():
            self.ax1.plot(self.xs[key], self.ys[key], color='k', picker=True,
                          label=key)


            zap_xs = self.info[key]['zap'].xs

            # Eliminate zapped intervals
            plt.draw()
            good[key] = np.ones(len(self.xs[key]), dtype=bool)
            if len(zap_xs) >= 2:
                intervals = list(zip(zap_xs[:-1:2], zap_xs[1::2]))
                for i in intervals:
                    good[key][np.logical_and(self.xs[key] >= i[0],
                                             self.xs[key] <= i[1])] = False
            fitpars = self.info[key]['fitpars']
            if len(self.info[key]['fitpars']) >= 2:
                model[key] = linear_fun(self.xs[key], *fitpars)
                self.ax1.plot(self.xs[key], model[key], color='b')
            else:
                model[key] = np.zeros(len(self.xs[key])) + np.min(self.ys[key])

        self.ax2.cla()
        self.ax2.axhline(0, ls='--', color='k')
        for key in self.xs.keys():
            self.ax2.plot(self.xs[key], self.ys[key] - model[key], color='grey', ls='-')
            self.ax2.plot(self.xs[key][good[key]], self.ys[key][good[key]] - model[key][good[key]], 'k-', lw=2)

        if self.xlabel is not None:
            self.ax2.set_xlabel(self.xlabel)
        plt.draw()

    def print_instructions(self):
        instructions = '''
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
    u     update plots with new selections
    B     subtract the baseline;
    r     reset baseline and zapping intervals, and fit parameters;
    q     quit

-------------------------------------------------------------
    '''

        print(instructions)

    def print_info(self):
        for key in self.info.keys():
            print(key + ':')
            if len(self.info[key]['zap'].xs) >= 2:
                print('  Zap intervals: ', list(zip(self.info[key]['zap'].xs[:-1:2],
                                                    self.info[key]['zap'].xs[1::2])))

            print('  Fit pars:      ', self.info[key]['fitpars'])

def select_data(xs, ys, title=None, xlabel=None):
    try:
        xs.keys()
    except:
        xs = {'Ch': xs}
        ys = {'Ch': ys}

    if title is None:
        title = 'Data selector'

    fig = plt.figure(title)
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 2], hspace=0)

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1], sharex=ax1)

    datasel = DataSelector(xs, ys, ax1, ax2, title=title, xlabel=xlabel)

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
            plt.close(self.ax.figure)
            self.fun(x, y, key)
        else:
            print(x, y, key)
        return x, y, key

    def plot_img(self):
        self.ax.imshow(self.img, origin='lower')


def test_select_data():
    times = np.arange(0, 100, 0.1)
    lc1 = np.random.normal(0, 0.3, len(times)) + times * 0.5 +\
        10 * np.exp(-(times - 30) ** 2 / 20 ** 2) + \
        3 * np.exp(-(times - 55) ** 2 / 6 ** 2) + \
        10 * np.exp(-(times - 80) ** 2 / 0.1 ** 2)
    lc2 = np.random.normal(0, 0.3, len(times)) + times * 0.3 +\
        10 * np.exp(-(times - 30) ** 2 / 20 ** 2) + \
        3 * np.exp(-(times - 55) ** 2 / 6 ** 2) + \
        10 * np.exp(-(times - 80) ** 2 / 0.1 ** 2)
    xs = {'Ch0': times, 'Ch1': times}
    ys = {'Ch0': lc1, 'Ch1': lc2}
    select_data(xs, ys)


def test_imageselector():
    img = np.random.normal(100, 100, (100, 100))
    fig = plt.figure()
    ax = fig.add_subplot(111)

    imagesel = ImageSelector(img, ax)

    plt.show()


