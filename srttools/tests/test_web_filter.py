"""Tests of the Bokeh filter page, calling its Python callbacks directly (no browser)."""

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
