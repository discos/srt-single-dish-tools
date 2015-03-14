from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from matplotlib import pyplot as plt
import numpy as np


class DataSelector:
    def __init__(self, fig):
        self.xs = []
        self.ys = []
        self.cid = fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.cid = fig.canvas.mpl_connect('key_press_event', self.on_key)
        self.counter = 0
        self.lines = []
        self.fig = fig

    def on_click(self, event):
        self.xs.append(event.xdata)
        self.ys.append(event.ydata)
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
        if event.key == 'r':
            for l in self.lines:
                l.remove()
            self.lines = []
            self.xs = []
            self.ys = []

            plt.draw()
        if event.key == 'q':
            plt.close(self.fig)


def select_data(xs, ys):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(xs, ys)
    ax.set_title('right click to get x intervals; r to reset; q to quit')

    datasel = DataSelector(fig)

    plt.show()

    return datasel.xs, datasel.ys


def test_select_data():
    times = np.arange(0, 100, 0.1)
    lc = np.random.poisson(10, len(times))
    x, y = select_data(times, lc)

    good = np.ones(len(times), dtype=bool)
    if len(x) >= 2:
        intervals = list(zip(x[:-1:2], x[1::2]))
        for i in intervals:
            good[np.logical_and(times >= i[0],
                                times <= i[1])] = False
    plt.plot(times, lc, 'c-')
    plt.plot(times[good], lc[good], 'b-', lw=3)
    plt.show()
