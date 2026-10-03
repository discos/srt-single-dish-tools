import logging
import warnings

import numpy as np
import pytest

from srttools.interactive_filter import (
    HAS_MPL,
    DataSelector,
    ImageSelector,
    PlotWarning,
    TestWarning,
    _disable_keymaps,
    _restore_keymaps,
    mask,
    select_data,
)

np.random.seed(1241347)


@pytest.mark.skipif("not HAS_MPL")
class TestImageSelector:
    @classmethod
    def setup_class(klass):
        import matplotlib.pyplot as plt

        klass.data = np.zeros((100, 100))
        klass.fig = plt.figure()
        klass.ax = plt.subplot()

        def fun(x, y, key):
            warnings.warn(f"It is working: {x}, {y}, {key}", TestWarning)

        klass.selector = ImageSelector(klass.data, klass.ax, test=True, fun=fun)

    def test_interactive_valid_data(self):
        fake_event = type("event", (), {})()
        fake_event.key = "q"
        fake_event.xdata, fake_event.ydata = (130, 30)

        retval = self.selector.on_key(fake_event)
        assert retval == (130, 30, "q")

    def test_interactive_invalid_data(self):
        fake_event = type("event", (), {})()
        fake_event.key = "b"
        fake_event.xdata, fake_event.ydata = (None, 30)

        retval = self.selector.on_key(fake_event)
        assert retval is None

    def test_interactive_fun(self):
        fake_event = type("event", (), {})()
        fake_event.key = "b"
        fake_event.xdata, fake_event.ydata = (130, 30)

        with pytest.warns(TestWarning) as record:
            retval = self.selector.on_key(fake_event)
        assert np.any(["It is working: 130, 30, b" in r.message.args[0] for r in record])
        assert retval == (130, 30, "b")

    @classmethod
    def teardown_class(klass):
        import matplotlib.pyplot as plt

        plt.close(klass.fig)


@pytest.mark.skipif("not HAS_MPL")
class TestDataSelector:
    @classmethod
    def setup_class(klass):
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        chans = ["scan1.fits", "scan2.fits"]
        klass.xs = {c: np.arange(30) for c in chans}
        klass.ys = {
            c: -(1**i) * np.random.normal(klass.xs[c] * 0.1, 0.1) + i for i, c in enumerate(chans)
        }

        gs = mpl.gridspec.GridSpec(2, 1)

        klass.fig = plt.figure()
        klass.ax0 = plt.subplot(gs[0])
        klass.ax1 = plt.subplot(gs[1])

        klass.selector = DataSelector(klass.xs, klass.ys, klass.ax0, klass.ax1, test=True)
        klass.selector.current = "scan1.fits"

    def test_interactive_zap_and_print_info(self, capsys):
        fake_event = type("event", (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ("z", 1, 3)
        with pytest.warns(TestWarning) as record:
            self.selector.on_key(fake_event)
        fake_event = type("event", (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ("z", 4, 3)
        with pytest.warns(TestWarning) as record:
            self.selector.on_key(fake_event)
        assert np.any(["I select a zap interval at 4" in r.message.args[0] for r in record])
        assert self.selector.info["scan1.fits"]["zap"].xs == [1, 4]
        assert self.selector.info["scan1.fits"]["zap"].ys == [3, 3]
        assert self.selector.zcounter == 2
        fake_event = type("event", (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ("P", 1, 3)
        self.selector.on_key(fake_event)
        out, err = capsys.readouterr()
        assert "scan1.fits" + ":" in out

        assert "Zap intervals:" in out
        assert "[(1, 4)]" in out

    def test_interactive_base(self):
        fake_event = type("event", (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ("b", 1, 3)
        with pytest.warns(TestWarning) as record:
            self.selector.on_key(fake_event)
        fake_event = type("event", (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ("b", 4, 3)
        with pytest.warns(TestWarning) as record:
            self.selector.on_key(fake_event)
        assert np.any(["I put a baseline mark at 4" in r.message.args[0] for r in record])
        assert self.selector.info["scan1.fits"]["base"].xs == [1, 4]
        assert self.selector.info["scan1.fits"]["base"].ys == [3, 3]
        assert self.selector.bcounter == 2

    def test_print_instructions(self, capsys):
        fake_event = type("event", (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ("h", 1, 3)
        self.selector.on_key(fake_event)
        out, err = capsys.readouterr()
        assert "Interactive plotter." in out
        assert "z     create zap intervals" in out

    def test_update(self):
        fake_event = type("event", (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ("u", 1, 3)
        with pytest.warns(PlotWarning) as record:
            self.selector.on_key(fake_event)
        assert np.any(["I plotted all" in r.message.args[0] for r in record])

    def test_flag(self, caplog):
        assert self.selector.info["scan1.fits"]["FLAG"] is None
        fake_event = type("event", (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ("x", 1, 3)
        self.selector.on_key(fake_event)
        assert "Scan was flagged" in caplog.text
        assert self.selector.info["scan1.fits"]["FLAG"] is True

    def test_flag_otherscan(self, caplog):
        self.selector.current = "scan2.fits"
        assert self.selector.info["scan2.fits"]["FLAG"] is None
        fake_event = type("event", (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ("x", 1, 3)
        with pytest.warns(PlotWarning) as record:
            self.selector.on_key(fake_event)
            assert "Scan was flagged" in caplog.text
            fake_event.key, fake_event.xdata, fake_event.ydata = ("u", 1, 3)
            self.selector.on_key(fake_event)
        assert np.any(["I plotted all" in r.message.args[0] for r in record])
        assert self.selector.info["scan2.fits"]["FLAG"] is True
        self.selector.current = "scan1.fits"

    def test_unflag(self, caplog):
        fake_event = type("event", (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ("v", 1, 3)
        with pytest.warns(PlotWarning) as record:
            self.selector.on_key(fake_event)
            assert "Scan was unflagged" in caplog.text
            fake_event.key, fake_event.xdata, fake_event.ydata = ("u", 1, 3)
            self.selector.on_key(fake_event)
        assert np.any(["I plotted all" in r.message.args[0] for r in record])
        assert self.selector.info["scan1.fits"]["FLAG"] is False

    def test_reset(self):
        fake_event = type("event", (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ("b", 1, 3)
        self.selector.on_key(fake_event)
        self.selector.current = "scan2.fits"
        fake_event = type("event", (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ("z", 1, 3)
        self.selector.on_key(fake_event)
        self.selector.current = "scan1.fits"
        fake_event = type("event", (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ("r", 1, 3)
        self.selector.on_key(fake_event)
        assert self.selector.info["scan2.fits"]["FLAG"] is None
        assert self.selector.info["scan1.fits"]["base"].xs == []
        assert self.selector.info["scan1.fits"]["zap"].xs == []
        assert self.selector.info["scan1.fits"]["fitpars"][0] == 0

    def test_subtract_baseline_one_interval(self):
        fake_event = type("event", (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ("b", 1, 3)
        self.selector.on_key(fake_event)
        fake_event = type("event", (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ("b", 4, 3)
        self.selector.on_key(fake_event)

        fake_event = type("event", (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ("B", 1, 3)
        with pytest.warns(TestWarning) as record:
            self.selector.on_key(fake_event)
        assert np.any(["I subtract" in r.message.args[0] for r in record])
        assert self.selector.info["scan1.fits"]["fitpars"][1] != 0

    def test_subtract_baseline_no_interval(self):
        # Reset all
        fake_event = type("event", (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ("r", 1, 3)
        self.selector.on_key(fake_event)
        # Then fit
        fake_event = type("event", (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ("B", 1, 3)
        with pytest.warns(TestWarning) as record:
            self.selector.on_key(fake_event)
        assert np.any(["I subtract" in r.message.args[0] for r in record])
        assert self.selector.info["scan1.fits"]["fitpars"][1] == 0

    def test_align_all(self):
        fake_event = type("event", (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ("A", 1, 3)
        with pytest.warns(TestWarning) as record:
            self.selector.on_key(fake_event)
        assert np.any("I aligned all" in rec.message.args[0] for rec in record)

    def test_quit(self, caplog):
        fake_event = type("event", (), {})()
        fake_event.key, fake_event.xdata, fake_event.ydata = ("q", 1, 3)
        self.selector.on_key(fake_event)
        assert "Closing all figures and quitting." in caplog.text

    def test_select_data(self):
        info = select_data(self.xs, self.ys, test=True)
        assert info["scan1.fits"]["zap"].xs == []
        assert info["scan1.fits"]["base"].xs == []

    @classmethod
    def teardown_class(klass):
        import matplotlib.pyplot as plt

        plt.close(klass.fig)


def _send_key(fig, ax, key, xdata, ydata):
    """Send a real key press at data coordinates through the figure canvas."""
    from matplotlib.backend_bases import KeyEvent

    x, y = ax.transData.transform((xdata, ydata))
    fig.canvas.callbacks.process(
        "key_press_event", KeyEvent("key_press_event", fig.canvas, key, x, y)
    )


def _send_click(fig, ax, xdata, ydata):
    """Send a real mouse click at data coordinates (it also triggers picking)."""
    from matplotlib.backend_bases import MouseButton, MouseEvent

    x, y = ax.transData.transform((xdata, ydata))
    event = MouseEvent("button_press_event", fig.canvas, x, y, button=MouseButton.LEFT)
    fig.canvas.callbacks.process("button_press_event", event)


@pytest.fixture
def selector():
    """A DataSelector on two scans, driven by real matplotlib events."""
    import matplotlib.pyplot as plt

    xs = {"scan1.fits": np.arange(30.0), "scan2.fits": np.arange(30.0)}
    ys = {"scan1.fits": np.zeros(30), "scan2.fits": np.full(30, 10.0)}
    fig = plt.figure()
    ax1, ax2 = fig.subplots(2, 1, sharex=True)
    sel = DataSelector(xs, ys, ax1, ax2, test=True)
    fig.canvas.draw()
    yield sel
    plt.close(fig)


@pytest.mark.skipif("not HAS_MPL")
class TestDataSelectorEvents:
    def test_pick_selects_scan(self, selector):
        """Clicking on a scan in the upper panel selects it."""
        _send_click(selector.ax1.figure, selector.ax1, 15, 10)
        assert selector.current == "scan2.fits"

    def test_zap_with_keys(self, selector):
        """Each key press adds one point; intervals work in any order."""
        fig, ax = selector.ax1.figure, selector.ax1
        _send_click(fig, ax, 15, 0)
        with pytest.warns(TestWarning):
            _send_key(fig, ax, "z", 20, 0)
            _send_key(fig, ax, "z", 10, 0)
        zap = selector.info["scan1.fits"]["zap"].xs
        assert np.allclose(zap, [20, 10], atol=0.5)
        good = mask(selector.xs["scan1.fits"], zap)
        assert not np.any(good[11:20]) and np.all(good[:10]) and np.all(good[21:])

    def test_zap_markers_survive_redraw(self, selector):
        """The zap markers of the selected scan are still there after an update."""
        fig, ax = selector.ax1.figure, selector.ax1
        _send_click(fig, ax, 15, 0)
        with pytest.warns(TestWarning):
            _send_key(fig, ax, "z", 10, 0)
            _send_key(fig, ax, "z", 20, 0)
        with pytest.warns(PlotWarning):
            _send_key(fig, ax, "u", 15, 0)
        vertical = [line for line in ax.get_lines() if np.ptp(line.get_xdata()) == 0]
        assert len(vertical) == 2

    def test_matplotlib_keymaps_disabled_and_restored(self):
        """Matplotlib's own shortcuts on our keys can be switched off and back on."""
        import matplotlib as mpl

        saved = _disable_keymaps(DataSelector.keys)
        try:
            assert "r" not in mpl.rcParams["keymap.home"]
            assert "q" not in mpl.rcParams["keymap.quit"]
        finally:
            _restore_keymaps(saved)
        assert "r" in mpl.rcParams["keymap.home"]
        assert "q" in mpl.rcParams["keymap.quit"]

    def test_selector_kept_alive(self):
        """The selector is not garbage-collected while its figure is open."""
        import gc
        import weakref

        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax1, ax2 = fig.subplots(2, 1)
        ref = weakref.ref(
            DataSelector({"a": np.arange(30.0)}, {"a": np.zeros(30)}, ax1, ax2, test=True)
        )
        gc.collect()
        assert ref() is not None
        plt.close(fig)

    def test_select_data_does_not_hang(self):
        """With a non-interactive backend, select_data returns instead of waiting forever."""
        import threading

        import matplotlib as mpl

        if mpl.get_backend().lower() != "agg":
            pytest.skip("Needs the Agg backend")

        result = {}

        def run():
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result["info"] = select_data(np.arange(30.0), np.zeros(30))

        thread = threading.Thread(target=run, daemon=True)
        thread.start()
        thread.join(20)
        assert not thread.is_alive()
        assert result["info"] == {}
