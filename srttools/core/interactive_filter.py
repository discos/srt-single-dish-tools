from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np


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
    def __init__(self, xs, ys, ax1, ax2):
        self.xs = xs
        self.ys = ys
        self.ax1 = ax1
        self.ax2 = ax2
        self.info = {}
        self.info['zap'] = intervals()
        self.lines = []
        self.plot_all()

        ax1.figure.canvas.mpl_connect('button_press_event', self.on_click)
        ax1.figure.canvas.mpl_connect('key_press_event', self.on_key)
        self.counter = 0

    def on_click(self, event):
        pass

    def zap(self, event):
        '''Create a zap interval'''
        self.info['zap'].add([event.xdata, event.ydata])
        self.counter += 1
        color = 'r'
        if self.counter % 2 == 1:
            ls = '-'
        else:
            ls = '--'
        line = self.ax1.axvline(event.xdata, color=color, ls=ls)
        self.lines.append(line)
        plt.draw()

    def on_key(self, event):
        '''What to do when the keyboard is used'''

        if event.key == 'z':
            self.zap(event)
        if event.key == 'u':
            self.plot_all()
        elif event.key == 'r':
            for l in self.lines:
                l.remove()
            self.lines = []
            self.info['zap'].clear()
            self.plot_all()
        elif event.key == 'q':
            plt.close(self.ax1.figure)
        else:
            pass

    def plot_all(self):
        for l in self.lines:
            l.remove()
        self.lines = []
        self.ax1.cla()
        self.ax1.plot(self.xs, self.ys)

        zap_xs = self.info['zap'].xs

        plt.setp(self.ax1.get_xticklabels(), visible=False)
        plt.draw()
        good = np.ones(len(self.xs), dtype=bool)
        if len(zap_xs) >= 2:
            intervals = list(zip(zap_xs[:-1:2], zap_xs[1::2]))
            for i in intervals:
                good[np.logical_and(self.xs >= i[0],
                                    self.xs <= i[1])] = False
        self.ax2.cla()
        self.ax2.plot(self.xs, self.ys, 'c-')
        self.ax2.plot(self.xs[good], self.ys[good], 'b-', lw=3)
        plt.draw()


def select_data(xs, ys):
    instructions = '''Interactive plotter.


    Point and z: create zap intervals

    r:           reset;

    q:           quit

    u:           update plots with new selections
    '''

    print(instructions)

    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 2], hspace=0)

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1], sharex=ax1)

    datasel = DataSelector(xs, ys, ax1, ax2)

    plt.show()

    return datasel.info


def test_select_data():
    times = np.arange(0, 100, 0.1)
    lc = np.random.poisson(10, len(times))
    info = select_data(times, lc)
