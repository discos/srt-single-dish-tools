from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from matplotlib import pyplot as plt
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
    def __init__(self, fig):
        self.info = {}
        self.info['zap'] = intervals()

        self.cid = fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.cid = fig.canvas.mpl_connect('key_press_event', self.on_key)
        self.counter = 0
        self.lines = []
        self.fig = fig

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
        line = plt.gca().axvline(event.xdata, color=color, ls=ls)
        self.lines.append(line)
        plt.draw()

    def on_key(self, event):
        '''What to do when the keyboard is used'''

        if event.key == 'z':
            self.zap(event)
        elif event.key == 'r':
            for l in self.lines:
                l.remove()
            self.lines = []
            self.info['zap'].clear()
            plt.draw()
        elif event.key == 'q':
            plt.close(self.fig)
        else:
            pass


def select_data(xs, ys):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(xs, ys)
    ax.set_title('Point and z to create zap intervals; r to reset; q to quit')

    datasel = DataSelector(fig)

    plt.show()

    return datasel.info


def test_select_data():
    times = np.arange(0, 100, 0.1)
    lc = np.random.poisson(10, len(times))
    info = select_data(times, lc)

    x = info['zap'].xs

    good = np.ones(len(times), dtype=bool)
    if len(x) >= 2:
        intervals = list(zip(x[:-1:2], x[1::2]))
        for i in intervals:
            good[np.logical_and(times >= i[0],
                                times <= i[1])] = False
    plt.plot(times, lc, 'c-')
    plt.plot(times[good], lc[good], 'b-', lw=3)
    plt.show()
