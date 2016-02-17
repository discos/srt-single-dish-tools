# -*- coding: utf-8 -*-

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from matplotlib import pyplot as plt
import numpy as np
from .interactive_filter import select_data, ImageSelector


def test_select_data():
    plt.ioff()
    times = np.arange(0, 100, 0.1)
    lc1 = np.random.normal(0, 0.3, len(times)) + times * 0.5 +\
        10 * np.exp(-(times - 30) ** 2 / 20 ** 2) + \
        3 * np.exp(-(times - 55) ** 2 / 6 ** 2) + \
        10 * np.exp(-(times - 80) ** 2 / 0.1 ** 2)
    lc2 = np.random.normal(0, 0.3, len(times)) + times * 0.3 +\
        10 * np.exp(-(times - 30) ** 2 / 20 ** 2) + \
        3 * np.exp(-(times - 55) ** 2 / 6 ** 2) + \
        10 * np.exp(-(times - 80) ** 2 / 0.1 ** 2)
    lc3 = np.random.normal(0, 0.3, len(times)) + times * (-0.3) +\
        10 * np.exp(-(times - 30) ** 2 / 20 ** 2) + \
        3 * np.exp(-(times - 55) ** 2 / 6 ** 2) + \
        10 * np.exp(-(times - 80) ** 2 / 0.1 ** 2)
    xs = {'Ch0': times, 'Ch1': times, 'Ch2': times}
    ys = {'Ch0': lc1, 'Ch1': lc2, 'Ch2': lc3}
    select_data(xs, ys)
    plt.show()


def test_imageselector():
    plt.ioff()
    img = np.random.poisson(100, (100, 100))
    fig = plt.figure()
    ax = fig.add_subplot(111)

    imagesel = ImageSelector(img, ax)

    plt.show()
