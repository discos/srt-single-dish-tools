from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy.random as ra
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib as mpl


class FeedAligner():
    '''Align images from different feeds, interactively
    Inputs:
        data:       a dictionary containing images to align and extents,
                    like this: {'img0': [img0, xextent0, yextent0],
                                'img1': [img1, xextent1, yextent1], ...}
        ax:         a pyplot.axis instance where the image will be plotted
    '''
    def __init__(self, data, ax):
        self.data = data
        self.ax = ax
        self.xoffsets = dict(zip(data[0].keys(),
                                 np.zeros(len(data[0].keys()))))
        self.yoffsets = dict(zip(data[0].keys(),
                                 np.zeros(len(data[0].keys()))))
        self.radius = 0
        self.angle = 0
        self.mode = 'ROUGH'
        self.variable_to_change = 'A'
        self.step = 1.
        self.ax.figure.canvas.mpl_connect('key_press_event', self.on_key)
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
            self.radius += sign * self.step
        if self.variable_to_change == 'A':
            self.angle += sign * self.step
        self.recalculate_offsets()

    def print_variables(self):
        print('Radius: {}; Angle: {}'.format(self.radius, self.angle))

    def recalculate_offsets(self):
        import re
        channel_re = re.compile(r'^Ch([0-9]+)$')
        img, xextent, yextent = self.data
        xoffsets, yoffsets = \
            standard_offsets(self.radius, self.angle,
                             len(img.keys()))
        for k in img.keys():
            feed = int(channel_re.match(k).group(1))

            self.xoffsets[k] = xoffsets[feed]
            self.yoffsets[k] = yoffsets[feed]

        self.plot_imgs()

    def switch_mode(self):
        if self.mode == 'ROUGH':
            self.mode = 'FINE'
            self.step = 0.01
        else:
            self.mode = 'ROUGH'
            self.step = 1.

    def plot_imgs(self):
        self.ax.cla()
        img, xextent, yextent = self.data
        for ik, k in enumerate(img.keys()):
            xoff = self.xoffsets[k]
            yoff = self.yoffsets[k]
            self.ax.imshow(img[k], origin='lower',
                           extent=[xextent[0] + xoff, xextent[1] + xoff,
                                   yextent[0] + yoff, yextent[1] + yoff],
                           alpha=1/7)
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

#    print(np.degrees(feed_angles))
    xoffsets = radii * np.cos(feed_angles)
    yoffsets = radii * np.sin(feed_angles)
    return xoffsets, yoffsets


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

    mean = [-xoff, -yoff]
    cov = [[1, 0], [0, 1]]  # diagonal covariance, points lie on x or y-axis

    x, y = ra.multivariate_normal(mean, cov, 5000).T

    hist, _, _ = np.histogram2d(x, y, bins=[xbins, ybins])

    img = hist.T
    channel = 'Ch{}'.format(i)
    imgs[channel] = img
    plt.imshow(img, origin='lower', extent=[-30, 30, -30, 30], alpha=1/7,
               cmap=cmaps[i])
    plt.scatter([-xoff], [-yoff])


#fig = plt.figure('Gridspec', figsize=(9,9))
#gs = GridSpec(5, 5)
#
#radii = np.linspace(1, 5, 5)
#angles = np.linspace(-20, 20, 5)
#
#for ir, r in enumerate(radii):
#    for ia, a in enumerate(angles):
#        ax = plt.subplot(gs[ir, ia])
#        ax.set_title('Radius: {}, angle: {}'.format(r, a))
#        xoffsets, yoffsets = standard_offsets(r, a)
#        for i in range(7):
#            xoff = xoffsets[i]
#            yoff = yoffsets[i]
#            channel = 'Ch{}'.format(i)
#            img = imgs[channel]
#            ax.imshow(img, origin='lower',
#                      extent=[-15 + xoff, 15 + xoff, -15 + yoff, 15 + yoff],
#                      alpha=1/7,
#                      cmap=cmaps[i])
#            ax.scatter([-real_xoffsets[i] + xoff],
#                       [-real_yoffsets[i] + yoff], color=colors[i])
#
#plt.show()

fig = plt.figure('Interactive')
ax = fig.add_subplot(111)
FeedAligner([imgs, [-30, 30], [-30, 30]], ax)
plt.show()