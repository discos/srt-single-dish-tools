from __future__ import (absolute_import, division,
                        print_function)
import numpy as np
from srttools.interactive_filter import ImageSelector, DataSelector
from srttools.interactive_filter import select_data, HAS_MPL
from srttools.interactive_filter import TestWarning, PlotWarning
import warnings
import pytest

np.random.seed(1241347)


@pytest.mark.skipif('not HAS_MPL')
class TestImageSelector(object):
    @classmethod
    def setup_class(klass):
        import matplotlib.pyplot as plt
        klass.data = np.zeros((100, 100))
        klass.ax = plt.subplot()

        def fun(x, y, key):
            warnings.warn("It is working: {}, {}, {}".format(x, y, key),
                          TestWarning)
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

        with pytest.warns(TestWarning) as record:
            retval = self.selector.on_key(fake_event)
        assert "It is working: 130, 30, b" in record[0].message.args[0]
        assert retval == (130, 30, 'b')


@pytest.mark.skipif('not HAS_MPL')
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

    def test_interactive_zap_and_print_info(self, capsys):
        fake_event = type('event', (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ('z', 1, 3)
        with pytest.warns(TestWarning) as record:
            self.selector.on_key(fake_event)
        fake_event = type('event', (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ('z', 4, 3)
        with pytest.warns(TestWarning) as record:
            self.selector.on_key(fake_event)
        assert "I select a zap interval at 4" in record[0].message.args[0]
        assert self.selector.info['scan1.fits']['zap'].xs == [1, 4]
        assert self.selector.info['scan1.fits']['zap'].ys == [3, 3]
        assert self.selector.zcounter == 2
        fake_event = type('event', (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ('P', 1, 3)
        self.selector.on_key(fake_event)
        out, err = capsys.readouterr()
        assert "scan1.fits" + ":" in out

        assert "Zap intervals:" in out
        assert "[(1, 4)]" in out

    def test_interactive_base(self):
        fake_event = type('event', (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ('b', 1, 3)
        with pytest.warns(TestWarning) as record:
            self.selector.on_key(fake_event)
        fake_event = type('event', (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ('b', 4, 3)
        with pytest.warns(TestWarning) as record:
            self.selector.on_key(fake_event)
        assert "I put a baseline mark at 4" in record[0].message.args[0]
        assert self.selector.info['scan1.fits']['base'].xs == [1, 4]
        assert self.selector.info['scan1.fits']['base'].ys == [3, 3]
        assert self.selector.bcounter == 2

    def test_print_instructions(self, capsys):
        fake_event = type('event', (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ('h', 1, 3)
        self.selector.on_key(fake_event)
        out, err = capsys.readouterr()
        assert "Interactive plotter." in out
        assert 'z     create zap intervals' in out

    def test_update(self):
        fake_event = type('event', (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ('u', 1, 3)
        with pytest.warns(PlotWarning) as record:
            self.selector.on_key(fake_event)
        assert "I plotted all" in record[0].message.args[0]

    def test_flag(self, capsys):
        assert self.selector.info['scan1.fits']['FLAG'] is None
        fake_event = type('event', (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ('x', 1, 3)
        self.selector.on_key(fake_event)
        out, err = capsys.readouterr()
        assert "Scan was flagged" in out
        assert self.selector.info['scan1.fits']['FLAG'] is True

    def test_flag_otherscan(self, capsys):
        self.selector.current = 'scan2.fits'
        assert self.selector.info['scan2.fits']['FLAG'] is None
        fake_event = type('event', (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ('x', 1, 3)
        with pytest.warns(PlotWarning) as record:
            self.selector.on_key(fake_event)
            out, err = capsys.readouterr()
            assert "Scan was flagged" in out
            fake_event.key, fake_event.xdata, fake_event.ydata = ('u', 1, 3)
            self.selector.on_key(fake_event)
        assert "I plotted all" in record[0].message.args[0]
        assert self.selector.info['scan2.fits']['FLAG'] is True
        self.selector.current = 'scan1.fits'

    def test_unflag(self, capsys):
        fake_event = type('event', (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ('v', 1, 3)
        with pytest.warns(PlotWarning) as record:
            self.selector.on_key(fake_event)
            out, err = capsys.readouterr()
            assert "Scan was unflagged" in out
            fake_event.key, fake_event.xdata, fake_event.ydata = ('u', 1, 3)
            self.selector.on_key(fake_event)
        assert "I plotted all" in record[0].message.args[0]
        assert self.selector.info['scan1.fits']['FLAG'] is False

    def test_reset(self):
        fake_event = type('event', (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ('b', 1, 3)
        self.selector.on_key(fake_event)
        self.selector.current = 'scan2.fits'
        fake_event = type('event', (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ('z', 1, 3)
        self.selector.on_key(fake_event)
        self.selector.current = 'scan1.fits'
        fake_event = type('event', (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ('r', 1, 3)
        self.selector.on_key(fake_event)
        assert self.selector.info['scan2.fits']['FLAG'] is None
        assert self.selector.info['scan1.fits']['base'].xs == []
        assert self.selector.info['scan1.fits']['zap'].xs == []
        assert self.selector.info['scan1.fits']['fitpars'][0] == 0

    def test_subtract_baseline_one_interval(self):
        fake_event = type('event', (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ('b', 1, 3)
        self.selector.on_key(fake_event)
        fake_event = type('event', (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ('b', 4, 3)
        self.selector.on_key(fake_event)

        fake_event = type('event', (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ('B', 1, 3)
        with pytest.warns(TestWarning) as record:
            self.selector.on_key(fake_event)
        assert "I subtract" in record[0].message.args[0]
        assert self.selector.info['scan1.fits']['fitpars'][1] != 0

    def test_subtract_baseline_no_interval(self):
        # Reset all
        fake_event = type('event', (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ('r', 1, 3)
        self.selector.on_key(fake_event)
        # Then fit
        fake_event = type('event', (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ('B', 1, 3)
        with pytest.warns(TestWarning) as record:
            self.selector.on_key(fake_event)
        assert "I subtract" in record[0].message.args[0]
        assert self.selector.info['scan1.fits']['fitpars'][1] == 0

    def test_align_all(self):
        fake_event = type('event', (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ('A', 1, 3)
        with pytest.warns(TestWarning) as record:
            self.selector.on_key(fake_event)
        assert "I aligned all" in record[0].message.args[0]

    def test_quit(self, capsys):
        fake_event = type('event', (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ('q', 1, 3)
        self.selector.on_key(fake_event)
        out, err = capsys.readouterr()
        assert "Closing all figures and quitting." in out

    def test_select_data(self):
        info = select_data(self.xs, self.ys, test=True)
        assert info['scan1.fits']['zap'].xs == []
        assert info['scan1.fits']['base'].xs == []
