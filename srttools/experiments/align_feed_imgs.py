from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy.random as ra
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from scipy.ndimage.interpolation import shift
import re
channel_re = re.compile(r'.*CH([0-9]+)$')


class FeedAligner():
    '''Align images from different feeds, interactively
    Inputs:
        imgs:       a dictionary containing images to align,
                    like this: {'Ch0': img0, 'Ch1': img1, ...}
        ax:         a pyplot.axis instance where the image will be plotted
    '''
    def __init__(self, imgs, ax, xextent=None, yextent=None, feeds=None,
                 radius=None, angle=None):

        self.imgs = imgs
        self.ax = ax
        if xextent is None:
            xextent = [-1, 1]
        if yextent is None:
            yextent = [-1, 1]
        self.xextent = xextent
        self.yextent = yextent
        keys = imgs.keys()
        if feeds is None:
            idxs = [int(channel_re.match(k).group(1)) for k in keys]
            feeds = dict(zip(keys, idxs))
        self.feeds = feeds
        self.xoffsets = dict(zip(keys, np.zeros(len(keys))))
        self.yoffsets = dict(zip(keys, np.zeros(len(keys))))
        if radius is None:
            radius = 0
        self.radius = radius
        if angle is None:
            angle = 0
        self.angle = 0
        self.mode = 'ROUGH'
        self.variable_to_change = 'A'
        self.radius_step = (xextent[1] - xextent[0]) / 50
        self.angle_step = 5.
        self.step_factor = 1.
        self.ax.figure.canvas.mpl_connect('key_press_event', self.on_key)
        self.plotted = None
        self.plot_imgs()

    def on_key(self, event):
        x, y = event.xdata, event.ydata
        key = event.key

        if key == 'm':
            self.switch_mode()
        if key == 'r':
            self.variable_to_change = 'R'
        if key == 'a':
            self.variable_to_change = 'A'
        if key == 'down':
            self.makestep(-1)
        if key == 'up':
            self.makestep(1)

        self.print_variables()
        return x, y, key

    def makestep(self, sign):
        if self.variable_to_change == 'R':
            self.radius += sign * self.step_factor * self.radius_step
        if self.variable_to_change == 'A':
            self.angle += sign * self.step_factor * self.angle_step
        self.recalculate_offsets()

    def print_variables(self):
        print('Radius: {}; Angle: {}'.format(self.radius, self.angle))

    def recalculate_offsets(self):
        xoffsets, yoffsets = \
            standard_offsets(self.radius, -self.angle,
                             len(self.imgs.keys()))
        for k in self.imgs.keys():
            feed = self.feeds[k]

            self.xoffsets[k] = xoffsets[feed]
            self.yoffsets[k] = yoffsets[feed]

        self.plot_imgs()

    def switch_mode(self):
        if self.mode == 'ROUGH':
            self.mode = 'FINE'
            self.step_factor = 0.05
        else:
            self.mode = 'ROUGH'
            self.step_factor = 1

    def plot_imgs(self):
        self.ax.cla()
        total_image = 0
        keys = list(self.imgs.keys())
#        self.fig = plt.figure('Increase')
#        gs = GridSpec(4, len(keys) // 2)

#        half = len(keys) // 2

        for ik, k in enumerate(keys):
#            col = ik % half
#            row = ik // half
            img = self.imgs[k]
            xextent = self.xextent
            yextent = self.yextent
            xoff = \
                self.xoffsets[k] / (xextent[1] - xextent[0]) * len(img[0, :])
            yoff = \
                self.yoffsets[k] / (yextent[1] - yextent[0]) * len(img[:, 0])

            shifted = shift(img, [-xoff, -yoff])
            total_image += shifted
#            inc_ax = plt.subplot(gs[row * 2, col])
#            inc_ax.imshow(total_image,
#                          extent=xextent + yextent,
#                          vmin=np.percentile(total_image, 40))
#            inc_ax = plt.subplot(gs[row * 2 + 1, col])
#            inc_ax.imshow(shifted,
#                          extent=xextent + yextent,
#                          vmin=np.percentile(self.imgs[keys[0]], 40))
        self.ax.imshow(total_image,
                       extent=xextent + yextent,
                       vmin=np.percentile(self.imgs[keys[0]], 40)*len(keys))
        plt.draw()


def standard_offsets(radius=3., angle=0, nfeeds=7):
    if nfeeds == 1:
        return [0], [0]

    # 0 for feed 0, radius for the other six
    radii = np.array([0] + [radius]*(nfeeds - 1))
    # Feeds 1--6 are at angles -60, -120, etc. Here I use angle 0 for
    # convenience for feed 0, but it has no effect since radii[0] is 0
    feed_angles = \
        -np.arange(0, nfeeds, 1) * np.pi * 2/(nfeeds - 1) + np.radians(angle)

    xoffsets = radii * np.cos(feed_angles)
    yoffsets = radii * np.sin(feed_angles)
    return xoffsets, yoffsets


def test_with_point_source():
        # make the colormaps
    colors = 'white,red,green,blue,magenta,cyan,yellow'.split(',')
    cmaps = []
    for i in range(7):
        cmap = mpl.colors.LinearSegmentedColormap.from_list('cmap{}'.format(i),
                                                            ['black', colors[i]],
                                                            256)
        cmaps.append(cmap)

    real_xoffsets, real_yoffsets = standard_offsets(2.4, 45)

    imgs = {}
    xbins = np.linspace(-30, 30, 101)
    ybins = np.linspace(-30, 30, 101)
    for i in range(7):
        xoff = real_xoffsets[i]
        yoff = real_yoffsets[i]

        mean = [-xoff - 3, -yoff - 1]
        mean2 = [-xoff + 4, -yoff + 5]
        cov = [[1, 0], [0, 1]]  # diagonal covariance, points lie on x or y-axis

        x, y = ra.multivariate_normal(mean, cov, 50000).T
        x2, y2 = ra.multivariate_normal(mean2, cov, 50000).T

        x3 = ra.uniform(-30, 30, 50000)
        y3 = ra.uniform(-30, 30, 50000)

        hist, _, _ = np.histogram2d(x, y, bins=[xbins, ybins])
        hist2, _, _ = np.histogram2d(x2, y2, bins=[xbins, ybins])
        hist3, _, _ = np.histogram2d(x3, y3, bins=[xbins, ybins])

        img = hist.T + hist2.T + hist3.T
        channel = 'CH{}'.format(i)
        imgs[channel] = img[::-1]
        plt.imshow(img[::-1], extent=[-30, 30, -30, 30], alpha=1/7,
                   cmap=cmaps[i])
        plt.scatter([-xoff], [-yoff])

    fig = plt.figure('Interactive')
    ax = fig.add_subplot(111)
    FeedAligner(imgs, ax, [-30, 30], [-30, 30])
    plt.show()


if __name__ == '__main__':
    import astropy.io.fits as pf
    hdulist = pf.open('try_2.fits')

    imgs = {}
    feeds = {}
    for h in hdulist[1:]:
        name = h.name
        img = h.data
        imgs[name] = img[::-1]
        chnum = int(channel_re.match(name).group(1))
        feeds[name] = chnum // 2

    nx = hdulist['IMGCH0'].header['NAXIS1']
    ny = hdulist['IMGCH0'].header['NAXIS2']

    dx = np.abs(hdulist['IMGCH0'].header['CDELT1'])
    dy = np.abs(hdulist['IMGCH0'].header['CDELT2'])

    cx = hdulist['IMGCH0'].header['CRPIX1']
    cy = hdulist['IMGCH0'].header['CRPIX2']

    xextent = [-cx * dx, (nx - cx) * dx]
    yextent = [-cy * dy, (ny - cy) * dy]

    print(xextent, yextent)
    fig = plt.figure('Interactive')
    ax = fig.add_subplot(111)
    FeedAligner(imgs, ax, xextent, yextent, feeds)
    plt.show()
