"""Tests of the Bokeh filter page, calling its Python callbacks directly (no browser)."""

import glob
import os
import sys
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
    assert glob.glob(os.path.join(scanset.meta["productdir"], "casa*.fits"))


def test_sdtimage_interactive_without_bokeh(scanset, tmp_path, monkeypatch):
    """Without bokeh, SDTimage --interactive stops with an error that says how to install it."""
    from srttools.imager import main_imager

    monkeypatch.setitem(sys.modules, "srttools.web_filter", None)
    outfile = str(tmp_path / "casa_dump.hdf5")
    scanset.write(outfile, overwrite=True)
    with pytest.raises(ImportError, match="pip install bokeh"):
        main_imager([outfile, "--interactive", "--noplot"])


def test_sdtpreprocess_interactive_and_productdir(scanset, tmp_path, monkeypatch):
    """SDTpreprocess --interactive opens the scan list on the cleaned scans, which
    --productdir sends to the directory given on the command line."""
    import srttools.web_filter
    from srttools.imager import main_preprocess

    calls = []

    def fake_serve(scanset, port, view):
        calls.append((len(scanset.scan_list), port, view, scanset.meta["productdir"]))

    monkeypatch.setattr(srttools.web_filter, "serve", fake_serve)
    config = scanset.meta["config_file"]
    out = str(tmp_path / "out")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        main_preprocess(["-c", config, "--interactive", "--port", "5124", "--productdir", out])
    assert calls == [(len(scanset.scan_list), 5124, "scans", out)]
    assert glob.glob(os.path.join(out, "**", "*.hdf5"), recursive=True)
    with pytest.raises(ValueError, match="needs the config file"):
        main_preprocess(["some_scan.fits", "--interactive"])
    monkeypatch.setitem(sys.modules, "srttools.web_filter", None)
    with pytest.raises(ImportError, match="SDTpreprocess --interactive needs bokeh"):
        main_preprocess(["-c", config, "--interactive"])


def test_sdtpreprocess_interactive_single_files(scanset, tmp_path, monkeypatch):
    """SDTpreprocess --interactive also works on single files (e.g. to clean one subscan
    again), with the config file still required to know where the products go."""
    import srttools.web_filter
    from srttools.imager import main_preprocess

    calls = []

    def fake_serve(scanset, port, view):
        calls.append((sorted(scanset.scan_list), view))

    monkeypatch.setattr(srttools.web_filter, "serve", fake_serve)
    files = sorted(scanset.scan_list)[:2]
    out = str(tmp_path / "out")
    args = files + ["-c", scanset.meta["config_file"], "--interactive", "--productdir", out]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        main_preprocess(args)
    assert calls == [(files, "scans")]
    assert len(glob.glob(os.path.join(out, "**", "*.hdf5"), recursive=True)) == 2


def test_sdtimage_directories_from_command_line(scanset, tmp_path, monkeypatch):
    """SDTimage --productdir/--intermediatedir: with a config file, images go to the first and
    processed scans and dump to the second; with a dump alone, they reach the web filter."""
    import srttools.read_config
    import srttools.web_filter
    from srttools.imager import main_imager

    # Start from a fresh config: the overrides must not leak into the other tests
    monkeypatch.setattr(srttools.read_config, "SRT_tools_config", None)
    prod, mid = tmp_path / "prod", tmp_path / "mid"
    config = scanset.meta["config_file"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        main_imager(
            ["-c", config, "--noplot", "--productdir", str(prod)] + ["--intermediatedir", str(mid)]
        )
    assert (prod / "casa.fits").exists()
    assert (mid / "casa_dump.hdf5").exists()
    assert glob.glob(str(mid / "**" / "*.hdf5"), recursive=True)

    calls = []

    def fake_serve(scanset, port):
        calls.append((scanset.meta["productdir"], scanset.meta["intermediatedir"]))

    monkeypatch.setattr(srttools.web_filter, "serve", fake_serve)
    prod2 = tmp_path / "prod2"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        main_imager(
            [str(mid / "casa_dump.hdf5"), "--interactive", "--noplot", "--productdir", str(prod2)]
        )
    assert calls == [(str(prod2), str(prod2 / "intermediate"))]
    assert (prod2 / "casa.fits").exists()


@pytest.mark.parametrize(
    "main",
    [
        "srttools.imager:main_imager",
        "srttools.imager:main_preprocess",
        "srttools.calibration:main_cal",
        "srttools.calibration:main_lcurve",
        "srttools.crossscan:main_cross",
        "srttools.rfistat:main_rfistat",
    ],
)
def test_config_tools_have_directory_options(main, capsys):
    """Every tool that loads scans from a config file accepts --productdir and --intermediatedir."""
    import importlib

    module, func = main.split(":")
    with pytest.raises(SystemExit):
        getattr(importlib.import_module(module), func)(["--help"])
    out = capsys.readouterr().out
    assert "--productdir" in out and "--intermediatedir" in out


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
    """Both views keep the old look (ds9-like log scale from the 20th percentile, gnuplot2), over filled pixels."""
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
        filled = (expo > 0) & np.isfinite(img)
        if kind == "Outliers":
            filled &= img != 0
            img[~filled] = np.mean(img[filled])
        else:
            img[~filled] = np.min(img[filled])
        expected = ds9_like_log_scale(img)
        assert np.allclose(app.image_source.data["image"][0], expected)
        assert np.isclose(app.color_mapper.low, np.percentile(expected[filled], 20))


def test_empty_pixels_do_not_set_map_scale(app):
    """An offset of the filled pixels (e.g. no baseline subtraction) leaves the map unchanged."""
    ch = app.channel
    images = app.scanset.images
    expo = np.asarray(images[f"{ch}-EXPO"])
    expo[:2] = 0  # make sure some pixels are empty
    images[ch][:2] = 0
    app.update_map()
    before = np.array(app.image_source.data["image"][0])
    low = app.color_mapper.low
    original = np.array(images[ch])
    try:
        images[ch][expo > 0] += 500
        app.update_map()
        assert np.allclose(app.image_source.data["image"][0], before)
        assert np.isclose(app.color_mapper.low, low)
    finally:
        images[ch][:] = original
        app.scanset.calculate_images()


@pytest.fixture
def scans_app(scanset):
    closed = []
    app = FilterApp(scanset, view="scans", on_done=lambda: closed.append(True))
    app.closed = closed
    yield app
    app.session.reset()


def test_scans_view_lists_and_steps_through_scans(scans_app):
    """Without a map, the page lists all subscans, plots the first, and Next plots the following one."""
    app = scans_app
    doc = Document()
    app(doc)
    assert not hasattr(app, "map_figure")
    ids = app.scan_source.data["scan_id"]
    assert len(ids) == len(app.scanset.scan_list)
    assert app.selected == ids[0]
    assert len(app.line_sources["dec"].data["xs"]) == 1
    assert app.scan_figures["dec"].visible and not app.scan_figures["ra"].visible
    app.step_scan(1)
    assert app.selected == ids[1]
    assert list(app.scan_source.selected.indices) == [1]
    app.step_scan(-5)  # stops at the top of the list
    assert app.selected == ids[0]
    assert "Scan 1 of" in app.info.text


def test_scans_view_edits_update_the_list(scans_app):
    """Flagging and zapping the plotted scan show up in the list; Undo clears them again."""
    app = scans_app
    app.step_scan(2)
    app.flag()
    assert app.scan_source.data["state"][2] == "flagged"
    assert "Apply" not in app.status.text
    app.undo()
    s = app._selected_subscan()
    app.add_range(s.coord.min(), np.median(s.coord))
    app.zap()
    assert 0.3 < app.scan_source.data["masked"][2] < 0.7
    assert app.scan_source.data["state"][2] == ""


def test_close_does_not_drop_edits_silently(scans_app):
    """With unsaved edits and no map, Close first asks; it closes at the second press."""
    app = scans_app
    app.flag()
    app._report_errors(app.done)()
    assert "would be lost" in app.status.text and not app.closed
    app._report_errors(app.done)()
    assert app.closed


def test_save_to_another_directory_from_page(scans_app, tmp_path):
    """A bad "Save to" directory is refused at once; a good one outside the config gets the edits."""
    app = scans_app
    assert app.save_dir.value == app.session.default_intermediatedir
    assert app.save_note.text == ""
    blocker = tmp_path / "a_file"
    blocker.write_text("")
    app.save_dir.value = str(blocker / "sub")
    assert "cannot be written" in app.save_note.text
    app.flag()
    app._report_errors(app.save)()
    assert "cannot be written" in app.status.text and app.session.n_unsaved == 1
    app.save_dir.value = str(tmp_path / "elsewhere")
    assert "intermediatedir :" in app.save_note.text
    app.save()
    assert app.session.n_unsaved == 0
    assert glob.glob(str(tmp_path / "elsewhere" / "**" / "*.hdf5"), recursive=True)
