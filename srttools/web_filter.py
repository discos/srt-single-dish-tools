"""Web page for the interactive filtering of a ScanSet, built with Bokeh.

The page shows the map on the left; clicking a pixel shows, on the right, the
subscans passing through it (Dec and RA scans separately). All edits go
through a :class:`srttools.filter_session.FilterSession`, so nothing is
written to disk until the user saves.

:class:`FilterApp` builds the page and holds the Python callbacks; it can be
served by a Bokeh server (``SDTimage --interactive``) or shown in a Jupyter
notebook. Bokeh is an optional dependency.
"""

import numpy as np
from bokeh.events import Tap
from bokeh.layouts import column, row
from bokeh.models import (
    ColumnDataSource,
    Div,
    LinearColorMapper,
    RadioButtonGroup,
    Select,
)
from bokeh.plotting import figure

from .filter_session import FilterSession
from .utils import ds9_like_log_scale

__all__ = ["FilterApp"]

LINE_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2"]
MAP_KINDS = ["Image", "Outliers"]


def _finite_image(img):
    """Copy of the image with NaNs and infinities replaced by the median of the rest."""
    img = np.array(img, dtype=float)
    bad = ~np.isfinite(img)
    if np.all(bad):
        return np.zeros_like(img)
    img[bad] = np.median(img[~bad])
    return img


class FilterApp:
    """The interactive filter page.

    Parameters
    ----------
    scanset : :class:`srttools.imager.ScanSet`
        The scan set to filter.

    Other Parameters
    ----------------
    session : :class:`srttools.filter_session.FilterSession`, optional
        An existing session; by default a new one is created.
    channel : str, optional
        The channel shown first (default: the first one).
    """

    def __init__(self, scanset, session=None, channel=None):
        self.scanset = scanset
        self.session = FilterSession(scanset) if session is None else session
        self.channel = channel if channel is not None else str(scanset.chan_columns[0])
        self.map_kind = MAP_KINDS[0]
        self.pixel = None
        self.subscans = []
        if scanset.images is None:
            scanset.calculate_images()

        self.image_source = ColumnDataSource(dict(image=[], x=[], y=[], dw=[], dh=[]))
        self.pixel_source = ColumnDataSource(dict(x=[], y=[]))
        self.line_sources = {
            d: ColumnDataSource(dict(xs=[], ys=[], scan_id=[], color=[])) for d in ("dec", "ra")
        }
        self.masked_sources = {d: ColumnDataSource(dict(x=[], y=[])) for d in ("dec", "ra")}

        self.map_figure = self._make_map_figure()
        self.scan_figures = {d: self._make_scan_figure(d) for d in ("dec", "ra")}

        self.channel_select = Select(
            title="Channel", value=self.channel, options=[str(c) for c in scanset.chan_columns]
        )
        self.channel_select.on_change("value", lambda attr, old, new: self.set_channel(new))
        self.map_kind_buttons = RadioButtonGroup(labels=MAP_KINDS, active=0)
        self.map_kind_buttons.on_change(
            "active", lambda attr, old, new: self.set_map_kind(MAP_KINDS[new])
        )
        self.info = Div(text="Click a pixel of the map to see the subscans through it.")

        self.update_map()
        self.layout = column(
            row(self.channel_select, self.map_kind_buttons),
            row(
                column(self.map_figure, self.info),
                column(self.scan_figures["dec"], self.scan_figures["ra"]),
            ),
        )

    # ----- Figures --------------------------------------------------------

    def _make_map_figure(self):
        fig = figure(
            width=420,
            height=420,
            match_aspect=True,
            tools="pan,wheel_zoom,box_zoom,reset,tap",
            title="Map (click a pixel)",
            x_axis_label="x (pixel)",
            y_axis_label="y (pixel)",
        )
        self.color_mapper = LinearColorMapper(palette="Inferno256")
        fig.image(
            image="image",
            x="x",
            y="y",
            dw="dw",
            dh="dh",
            source=self.image_source,
            color_mapper=self.color_mapper,
        )
        fig.rect(
            x="x",
            y="y",
            width=1,
            height=1,
            source=self.pixel_source,
            fill_alpha=0,
            line_color="cyan",
            line_width=2,
        )
        fig.on_event(Tap, lambda event: self.select_pixel(event.x, event.y))
        return fig

    def _make_scan_figure(self, direction):
        label = "Dec" if direction == "dec" else "RA"
        fig = figure(
            width=520,
            height=260,
            tools="pan,wheel_zoom,box_zoom,reset",
            title=f"{label} subscans through the pixel",
            x_axis_label=f"{label} (deg)",
            y_axis_label="Counts",
        )
        fig.multi_line(xs="xs", ys="ys", line_color="color", source=self.line_sources[direction])
        fig.scatter(
            x="x", y="y", size=3, color="gray", alpha=0.6, source=self.masked_sources[direction]
        )
        return fig

    # ----- Updates --------------------------------------------------------

    def _map_array(self):
        images = self.scanset.images
        ch = self.channel
        if self.map_kind == "Image":
            return ds9_like_log_scale(_finite_image(images[ch]))
        expo = np.asarray(images[f"{ch}-EXPO"])
        mean_expo = np.mean(expo[expo > 0]) if np.any(expo > 0) else 0
        stat = "Outliers" if mean_expo > 6 and f"{ch}-Outliers" in images else "Sdev"
        return _finite_image(images[f"{ch}-{stat}"])

    def update_map(self):
        """Redraw the map from the current images of the scan set."""
        img = self._map_array()
        ny, nx = img.shape
        self.image_source.data = dict(image=[img], x=[0], y=[0], dw=[nx], dh=[ny])
        self.color_mapper.low = float(np.percentile(img, 20 if self.map_kind == "Image" else 1))
        self.color_mapper.high = float(np.max(img))

    def update_subscans(self):
        """Redraw the subscans through the selected pixel."""
        if self.pixel is None:
            return
        self.subscans = self.session.subscans_through_pixel(*self.pixel, self.channel)
        for direction in ("dec", "ra"):
            lines = dict(xs=[], ys=[], scan_id=[], color=[])
            masked = dict(x=[], y=[])
            for i, s in enumerate(self.subscans):
                if s.direction != direction:
                    continue
                lines["xs"].append(s.coord)
                lines["ys"].append(np.where(s.good, s.values, np.nan))
                lines["scan_id"].append(s.scan_id)
                lines["color"].append(LINE_COLORS[i % len(LINE_COLORS)])
                masked["x"].extend(s.coord[~s.good])
                masked["y"].extend(s.values[~s.good])
            self.line_sources[direction].data = lines
            self.masked_sources[direction].data = masked
        x, y = self.pixel
        self.info.text = f"Pixel ({int(x)}, {int(y)}): {len(self.subscans)} subscans"

    def select_pixel(self, x, y):
        """Show the subscans through the pixel containing the point ``(x, y)``."""
        if x is None or y is None:
            return
        self.pixel = (np.floor(x) + 0.5, np.floor(y) + 0.5)
        self.pixel_source.data = dict(x=[self.pixel[0]], y=[self.pixel[1]])
        self.update_subscans()

    def set_channel(self, channel):
        """Show another channel."""
        self.channel = channel
        self.update_map()
        self.update_subscans()

    def set_map_kind(self, kind):
        """Show the image (``"Image"``) or the outlier map (``"Outliers"``)."""
        self.map_kind = kind
        self.update_map()

    def __call__(self, doc):
        """Add the page to a Bokeh document (used by the Bokeh server and in notebooks)."""
        doc.title = "SDTimage interactive filter"
        doc.add_root(self.layout)
