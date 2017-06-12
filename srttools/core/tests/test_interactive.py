from __future__ import (absolute_import, division,
                        print_function)
import numpy as np
from ..interactive_filter import ImageSelector, DataSelector, select_data
import warnings
import pytest


np.random.seed(1241347)


class TestImageSelector(object):
    @classmethod
    def setup_class(klass):
        import matplotlib.pyplot as plt
        klass.data = np.zeros((100, 100))
        klass.ax = plt.subplot()

        def fun(x, y, key):
            warnings.warn("It is working: {}, {}, {}".format(x, y, key),
                          UserWarning)
        klass.selector = ImageSelector(klass.data, klass.ax, test=True,
                                       fun=fun)

    def test_interactive_valid_data(self):
        fake_event = type('event', (), {})()
        fake_event.key = 'q'
        fake_event.xdata, fake_event.ydata = (130, 30)

        retval = self.selector.on_key(fake_event)
        assert retval == (130, 30, 'q')

    def test_interactive_invalid_data(self):
        fake_event = type('event', (), {})()
        fake_event.key = 'b'
        fake_event.xdata, fake_event.ydata = (None, 30)

        retval = self.selector.on_key(fake_event)
        assert retval is None

    def test_interactive_fun(self):
        fake_event = type('event', (), {})()
        fake_event.key = 'b'
        fake_event.xdata, fake_event.ydata = (130, 30)

        with pytest.warns(UserWarning) as record:
            retval = self.selector.on_key(fake_event)
        assert "It is working: 130, 30, b" in record[0].message.args[0]
        assert retval == (130, 30, 'b')


class TestDataSelector(object):
    @classmethod
    def setup_class(klass):
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        chans = ['scan1.fits', 'scan2.fits']
        klass.xs = {c: np.arange(30) for c in chans}
        klass.ys = {c: -1 ** i * np.random.normal(klass.xs[c] * 0.1, 0.1) + i
                    for i, c in enumerate(chans)}

        gs = mpl.gridspec.GridSpec(2, 1)

        klass.ax0 = plt.subplot(gs[0])
        klass.ax1 = plt.subplot(gs[1])

        klass.selector = DataSelector(klass.xs, klass.ys, klass.ax0, klass.ax1,
                                      test=True)
        klass.selector.current = 'scan1.fits'

    def test_interactive_zap(self):
        fake_event = type('event', (), {})()
        fake_event.key = 'z'
        fake_event.xdata, fake_event.ydata = (1, 3)
        with pytest.warns(UserWarning) as record:
            self.selector.on_key(fake_event)
        assert "I select a zap interval at 1" in record[0].message.args[0]
        assert self.selector.info['scan1.fits']['zap'].xs == \
            [fake_event.xdata]
        assert self.selector.info['scan1.fits']['zap'].ys == \
            [fake_event.ydata]
        assert self.selector.zcounter == 1

    def test_interactive_base(self):
        fake_event = type('event', (), {})()
        fake_event.key = 'b'
        fake_event.xdata, fake_event.ydata = (1, 3)
        with pytest.warns(UserWarning) as record:
            self.selector.on_key(fake_event)
        assert "I put a baseline mark at 1" in record[0].message.args[0]
        assert self.selector.info['scan1.fits']['base'].xs == \
            [fake_event.xdata]
        assert self.selector.info['scan1.fits']['base'].ys == \
            [fake_event.ydata]
        assert self.selector.bcounter == 1

    def test_subtract_baseline(self):
        fake_event = type('event', (), {})()
        fake_event.key = 'B'
        fake_event.xdata, fake_event.ydata = (1, 3)
        with pytest.warns(UserWarning) as record:
            self.selector.on_key(fake_event)
        assert "I subtract" in record[0].message.args[0]

    def test_align_all(self):
        fake_event = type('event', (), {})()
        fake_event.key = 'A'
        fake_event.xdata, fake_event.ydata = (1, 3)
        with pytest.warns(UserWarning) as record:
            self.selector.on_key(fake_event)
        assert "I aligned all" in record[0].message.args[0]

    def test_select_data(self):
        info = select_data(self.xs, self.ys, test=True)
        assert info['scan1.fits']['zap'].xs == []
        assert info['scan1.fits']['base'].xs == []
