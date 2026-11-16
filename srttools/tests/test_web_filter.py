"""Tests of the Bokeh filter page, calling its Python callbacks directly (no browser)."""

import glob
import os
import warnings

import numpy as np
import pytest

pytest.importorskip("bokeh")

from bokeh.document import Document  # noqa: E402

from srttools.imager import ScanSet  # noqa: E402
from srttools.web_filter import FilterApp  # noqa: E402

MEDICINA = os.path.join(os.path.dirname(__file__), "data", "medicina")


@pytest.fixture(scope="module")
def scanset(tmp_path_factory):
    """The Medicina Cas A map (51 Dec subscans), products in a temporary directory."""
    tmp = tmp_path_factory.mktemp("web")
    config = str(tmp / "casa.ini")
    with open(config, "w") as fobj:
        print(
            f"[local]\nworkdir = {tmp}\ndatadir = {MEDICINA}\n"
            f"productdir = {tmp / 'products'}\n\n[analysis]\n"
            "list_of_directories =\n    20221213-162927-22-22-cas_a/\npixel_size = 1.5\n",
            file=fobj,
        )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        scanset = ScanSet(config)
        scanset.calculate_images()
    return scanset


@pytest.fixture
def app(scanset):
    app = FilterApp(scanset)
    yield app
    app.session.reset()


def test_page_builds_with_map(app):
    """The page goes into a Bokeh document, and the map shows the image of the first channel."""
    doc = Document()
    app(doc)
    assert len(doc.roots) == 1
    img = app.image_source.data["image"][0]
    assert img.shape == app.scanset.images["Feed0_RCP"].shape
    assert np.all(np.isfinite(img))


def test_click_on_pixel_shows_subscans(app):
    """Clicking the brightest pixel shows the Dec subscans through it, none in RA."""
    img = app.scanset.images["Feed0_RCP"]
    y, x = np.unravel_index(np.nanargmax(img), img.shape)
    app.select_pixel(x + 0.3, y + 0.8)
    assert app.pixel == (x + 0.5, y + 0.5)
    lines = app.line_sources["dec"].data
    assert 1 <= len(lines["xs"]) <= 4
    assert len(app.line_sources["ra"].data["xs"]) == 0
    assert np.nanmax(np.concatenate(lines["ys"])) > 0.5 * np.nanmax(app.scanset["Feed0_RCP"])


def test_channel_and_map_kind(app):
    """Switching channel and map kind redraws the map from the matching image."""
    app.set_channel("Feed0_LCP")
    assert app.channel == "Feed0_LCP"
    app.set_map_kind("Outliers")
    img = app.image_source.data["image"][0]
    assert img.shape == app.scanset.images["Feed0_LCP"].shape
    assert np.all(np.isfinite(img))


def _click_peak(app):
    img = app.scanset.images[app.channel]
    y, x = np.unravel_index(np.nanargmax(img), img.shape)
    app.select_pixel(x, y)
    return app.line_sources["dec"].data["scan_id"][0]


def test_zap_apply_undo(app):
    """Zapping the peak of a subscan masks it, changes the map after Apply; Undo restores it."""
    sid = _click_peak(app)
    app.select_subscan(sid)
    s = app._selected_subscan()
    peak = s.coord[np.argmax(s.values)]
    app.add_range(peak - 0.05, peak + 0.05)
    before = np.array(app.scanset.images[app.channel])
    app.zap()
    assert len(app.masked_sources["dec"].data["x"]) > 0
    assert "1 unsaved" in app.edit_list.text
    app.apply()
    assert np.nansum(app.scanset.images[app.channel]) < np.nansum(before)
    app.undo()
    app.apply()
    assert np.allclose(app.scanset.images[app.channel], before, equal_nan=True)
    assert "No unsaved" in app.edit_list.text


def test_errors_shown_on_page(app):
    """A button pressed in the wrong state shows the reason on the page instead of failing."""
    app._report_errors(app.zap)()
    assert "Select a subscan first" in app.status.text


def test_save_from_page(app):
    """Save writes the edits and empties the list of unsaved edits."""
    sid = _click_peak(app)
    app.select_subscan(sid)
    app.flag()
    app.save()
    assert "No unsaved" in app.edit_list.text
    assert "Saved 1 edits" in app.status.text


def test_server_serves_independent_pages(scanset):
    """The server answers two page requests, each building its own widgets on a shared session."""
    import asyncio
    import threading

    from tornado.httpclient import AsyncHTTPClient

    from srttools.web_filter import _make_server

    codes = []

    def run():
        asyncio.set_event_loop(asyncio.new_event_loop())
        server, session = _make_server(scanset, port=0)
        server.start()

        async def get_twice():
            try:
                for _ in range(2):
                    response = await AsyncHTTPClient().fetch(f"http://localhost:{server.port}/")
                    codes.append(response.code)
            finally:
                server.io_loop.stop()

        server.io_loop.add_callback(get_twice)
        server.io_loop.start()
        server.stop()

    thread = threading.Thread(target=run)
    thread.start()
    thread.join(60)
    assert codes == [200, 200]


def test_sdtimage_interactive_calls_server(scanset, tmp_path, monkeypatch):
    """SDTimage --interactive hands the scan set to the web server, then makes the images."""
    import srttools.web_filter
    from srttools.imager import main_imager

    calls = []
    monkeypatch.setattr(
        srttools.web_filter, "serve", lambda scanset, port: calls.append((len(scanset), port))
    )
    outfile = str(tmp_path / "casa_dump.hdf5")
    scanset.write(outfile, overwrite=True)
    main_imager([outfile, "--interactive", "--port", "5123", "--noplot"])
    assert calls == [(len(scanset), 5123)]
    # The images are named after the config file
    assert glob.glob(scanset.meta["config_file"].replace(".ini", "*.fits"))


def test_click_elsewhere_after_zap(app):
    """Clicking a pixel the selected subscan does not cross deselects it; one it crosses keeps it."""
    sid = _click_peak(app)
    app.select_subscan(sid)
    s = app._selected_subscan()
    app.add_range(s.coord[0], s.coord[5])
    app.zap()
    # Same subscan, other pixel along it: still selected, and highlighted at its new position
    x0 = app.pixel[0]
    app.select_pixel(x0, 2)
    assert app.selected == sid
    source = app.line_sources["dec"]
    assert [source.data["scan_id"][i] for i in source.selected.indices] == [sid]
    # A pixel far from it in x: no selection, no marked intervals
    img = app.scanset.images[app.channel]
    app.select_pixel((x0 + img.shape[1] // 2) % img.shape[1], 2)
    assert app.selected is None and app.ranges == []
    assert source.selected.indices == []
    app._report_errors(app.zap)()
    assert "Select a subscan first" in app.status.text


def test_map_scaling_matches_old_display(app):
    """Both views use the old display's look: ds9-like log scale from the 20th percentile, gnuplot2."""
    from matplotlib import colormaps
    from matplotlib.colors import to_hex

    from srttools.utils import ds9_like_log_scale

    assert app.color_mapper.palette[-1] == to_hex(colormaps["gnuplot2"](1.0))
    ch = app.channel
    expo = np.asarray(app.scanset.images[f"{ch}-EXPO"])
    stat = "Outliers" if np.mean(expo[expo > 0]) > 6 else "Sdev"
    for kind, name in [("Image", ch), ("Outliers", f"{ch}-{stat}")]:
        app.set_map_kind(kind)
        img = np.array(app.scanset.images[name], dtype=float)
        if kind == "Outliers":
            bad = (img == 0) | ~np.isfinite(img)
            img[bad] = np.mean(img[~bad])
        expected = ds9_like_log_scale(img)
        assert np.allclose(app.image_source.data["image"][0], expected)
        assert np.isclose(app.color_mapper.low, np.percentile(expected, 20))
