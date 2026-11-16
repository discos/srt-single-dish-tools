"""Web page for the interactive filtering of a ScanSet, built with Bokeh.

The page shows the map on the left; clicking a pixel shows, on the right, the
subscans passing through it (Dec and RA scans separately). All edits go
through a :class:`srttools.filter_session.FilterSession`, so nothing is
written to disk until the user saves.

:class:`FilterApp` builds the page and holds the Python callbacks; it can be
served by a Bokeh server (``SDTimage --interactive``) or shown in a Jupyter
notebook. Bokeh is an optional dependency.
"""

import os
import socket

import numpy as np
from bokeh.events import DocumentReady, SelectionGeometry, Tap
from bokeh.layouts import column, row
from bokeh.models import (
    BoxSelectTool,
    Button,
    ColumnDataSource,
    CustomJS,
    Div,
    LinearColorMapper,
    RadioButtonGroup,
    Select,
    TapTool,
)
from bokeh.plotting import figure

from .filter_session import FilterSession
from .utils import ds9_like_log_scale

__all__ = ["FilterApp", "serve", "ssh_hint"]

LINE_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2"]
MAP_KINDS = ["Image", "Outliers"]


def _gnuplot2_palette():
    """The gnuplot2 colour map of the old matplotlib display, as 256 hex colours."""
    from matplotlib import colormaps
    from matplotlib.colors import to_hex

    cmap = colormaps["gnuplot2"]
    return [to_hex(cmap(i / 255)) for i in range(256)]


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
    on_done : callable, optional
        Called by the "Close and continue" button (e.g. to stop the server).
    """

    def __init__(self, scanset, session=None, channel=None, on_done=None):
        self.on_done = on_done
        self.scanset = scanset
        self.session = FilterSession(scanset) if session is None else session
        self.channel = channel if channel is not None else str(scanset.chan_columns[0])
        self.map_kind = MAP_KINDS[0]
        self.pixel = None
        self.subscans = []
        self.selected = None  # scan_id of the selected subscan
        self.ranges = []  # (lo, hi) intervals marked on the subscan plots, degrees
        if scanset.images is None:
            scanset.calculate_images()

        self.image_source = ColumnDataSource(dict(image=[], x=[], y=[], dw=[], dh=[]))
        self.pixel_source = ColumnDataSource(dict(x=[], y=[]))
        self.line_sources = {
            d: ColumnDataSource(dict(xs=[], ys=[], scan_id=[], color=[])) for d in ("dec", "ra")
        }
        self.masked_sources = {d: ColumnDataSource(dict(x=[], y=[])) for d in ("dec", "ra")}
        self.range_sources = {d: ColumnDataSource(dict(left=[], right=[])) for d in ("dec", "ra")}

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
        self.status = Div(text="")
        self.edit_list = Div(text="")
        # Read by the browser to warn before closing the page with unsaved edits
        self.unsaved_flag = Div(text="0", visible=False)

        self.buttons = {}
        for name, label, callback in [
            ("zap", "Zap range", self.zap),
            ("flag", "Flag", self.flag),
            ("unflag", "Unflag", self.unflag),
            ("baseline", "Fit baseline", self.fit_baseline),
            ("align", "Align all", self.align),
            ("clear", "Clear ranges", self.clear_ranges),
            ("undo", "Undo", self.undo),
            ("apply", "Apply", self.apply),
            ("save", "Save to productdir", self.save),
            ("done", "Close and continue", self.done),
        ]:
            button = Button(
                label=label,
                button_type="primary" if name == "save" else "default",
                width=140 if name in ("save", "done") else 110,
            )
            button.on_click(self._report_errors(callback))
            self.buttons[name] = button

        self.update_map()
        self.update_edit_list()
        b = self.buttons
        self.layout = column(
            row(self.channel_select, self.map_kind_buttons),
            row(
                column(self.map_figure, self.info),
                column(self.scan_figures["dec"], self.scan_figures["ra"]),
            ),
            row(b["zap"], b["flag"], b["unflag"], b["baseline"], b["align"], b["clear"]),
            row(b["undo"], b["apply"], b["save"], b["done"], self.status),
            self.edit_list,
            self.unsaved_flag,
        )

    def _report_errors(self, callback):
        """Wrap a button callback: show errors on the page instead of losing them."""

        def wrapped(*args):
            try:
                callback()
            except Exception as e:
                self.status.text = f"<span style='color:#c00'>{e}</span>"

        return wrapped

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
        self.color_mapper = LinearColorMapper(palette=_gnuplot2_palette())
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
        strips = fig.vstrip(
            x0="left",
            x1="right",
            source=self.range_sources[direction],
            fill_color="red",
            fill_alpha=0.12,
            line_alpha=0,
        )
        lines = fig.multi_line(
            xs="xs",
            ys="ys",
            line_color="color",
            line_width=1.5,
            source=self.line_sources[direction],
            selection_line_width=3,
            nonselection_line_alpha=0.4,
        )
        fig.scatter(
            x="x", y="y", size=3, color="gray", alpha=0.6, source=self.masked_sources[direction]
        )
        # The range tool selects only the strips: selecting lines would change the subscan
        range_tool = BoxSelectTool(
            dimensions="width", renderers=[strips], description="Mark an interval"
        )
        fig.add_tools(TapTool(renderers=[lines]), range_tool)
        fig.toolbar.active_drag = range_tool
        source = self.line_sources[direction]
        source.selected.on_change("indices", lambda attr, old, new: self._on_line_tap(source, new))
        fig.on_event(
            SelectionGeometry,
            lambda event: self._on_range(event),
        )
        return fig

    def _on_line_tap(self, source, indices):
        if indices:
            self.select_subscan(source.data["scan_id"][indices[0]])

    def _on_range(self, event):
        geometry = event.geometry
        if geometry.get("type") == "rect" and event.final:
            self.add_range(geometry["x0"], geometry["x1"])

    # ----- Updates --------------------------------------------------------

    def _map_array(self):
        """The map to show, on the ds9-like log scale, and the mask of the filled pixels.

        Empty pixels (no samples, or not finite) do not enter the scaling:
        they are set to the minimum (Image) or to the mean (Outliers) of
        the others.
        """
        images = self.scanset.images
        ch = self.channel
        expo = np.asarray(images[f"{ch}-EXPO"])
        if self.map_kind == "Image":
            img = np.array(images[ch], dtype=float)
            filled = (expo > 0) & np.isfinite(img)
            fill_value = np.min
        else:
            mean_expo = np.mean(expo[expo > 0]) if np.any(expo > 0) else 0
            stat = "Outliers" if mean_expo > 6 and f"{ch}-Outliers" in images else "Sdev"
            img = np.array(images[f"{ch}-{stat}"], dtype=float)
            filled = (expo > 0) & (img != 0) & np.isfinite(img)
            fill_value = np.mean
        if not np.any(filled):
            return np.zeros_like(img), filled
        img[~filled] = fill_value(img[filled])
        return ds9_like_log_scale(img), filled

    def update_map(self):
        """Redraw the map from the current images of the scan set."""
        img, filled = self._map_array()
        ny, nx = img.shape
        self.image_source.data = dict(image=[img], x=[0], y=[0], dw=[nx], dh=[ny])
        # Limits as in the old matplotlib display, but from the filled pixels only
        shown = img[filled] if np.any(filled) else img
        self.color_mapper.low = float(np.percentile(shown, 20))
        self.color_mapper.high = float(np.max(shown))

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
        # Keep the selected subscan if it crosses this pixel too, otherwise drop it
        if self._selected_subscan() is None:
            self.selected = None
            self.ranges = []
        for source in self.line_sources.values():
            ids = source.data["scan_id"]
            source.selected.indices = [ids.index(self.selected)] if self.selected in ids else []
        x, y = self.pixel
        self.info.text = f"Pixel ({int(x)}, {int(y)}): {len(self.subscans)} subscans"
        self.update_ranges()

    def update_ranges(self):
        """Draw the marked ranges on the plot of the selected subscan's direction."""
        direction = self._selected_subscan().direction if self.selected is not None else None
        for d in ("dec", "ra"):
            ranges = self.ranges if d == direction else []
            self.range_sources[d].data = dict(
                left=[min(r) for r in ranges], right=[max(r) for r in ranges]
            )

    def update_edit_list(self):
        """Show the edits not yet saved."""
        session = self.session
        n_saved = len(session.edits) - session.n_unsaved
        lines = []
        for edit in session.edits[n_saved:]:
            what = edit["kind"]
            if what == "zap":
                what += " {:.4f}-{:.4f} deg".format(*sorted(edit["range"]))
            elif what == "baseline":
                what += " (q={:.3g}, m={:.3g})".format(*edit["pars"])
            elif what == "align":
                what += f" {len(edit['pars'])} subscans"
            target = "" if edit["scan_id"] is None else f"subscan {edit['scan_id']}, "
            chan = edit["channel"] if edit["channel"] is not None else "all channels"
            lines.append(f"<li>{target}{chan}: {what}</li>")
        self.unsaved_flag.text = str(session.n_unsaved)
        if lines:
            self.edit_list.text = (
                f"<b>{session.n_unsaved} unsaved edits</b> (in memory)<ul>{''.join(lines)}</ul>"
            )
        else:
            self.edit_list.text = "No unsaved edits."

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

    # ----- Editing --------------------------------------------------------

    def _selected_subscan(self):
        for s in self.subscans:
            if s.scan_id == self.selected:
                return s
        return None

    def _need_selection(self):
        if self._selected_subscan() is None:
            raise ValueError("Select a subscan first (click on its line)")
        return self.selected

    def _after_edit(self, message):
        self.update_subscans()
        self.update_edit_list()
        self.status.text = message + " Press Apply to update the map."

    def select_subscan(self, scan_id):
        """Select the subscan to edit (marked intervals are cleared if it changes)."""
        if self.selected == int(scan_id):
            return
        self.selected = int(scan_id)
        self.ranges = []
        self.update_ranges()
        self.status.text = f"Selected subscan {self.selected}"

    def add_range(self, lo, hi):
        """Mark an interval (degrees) along the selected subscan, for zap or baseline fit."""
        self.ranges.append((float(lo), float(hi)))
        self.update_ranges()

    def clear_ranges(self):
        """Remove the marked intervals."""
        self.ranges = []
        self.update_ranges()

    def zap(self):
        """Mask the samples of the selected subscan inside the marked intervals."""
        sid = self._need_selection()
        if not self.ranges:
            raise ValueError("Mark an interval first (box select along the x axis)")
        for lo, hi in self.ranges:
            self.session.zap(sid, self.channel, lo, hi)
        self.ranges = []
        self._after_edit(f"Zapped subscan {sid}.")

    def flag(self):
        """Mask the selected subscan in all channels."""
        sid = self._need_selection()
        self.session.flag(sid)
        self._after_edit(f"Flagged subscan {sid}.")

    def unflag(self):
        """Remove all masks of the selected subscan."""
        sid = self._need_selection()
        self.session.unflag(sid)
        self._after_edit(f"Unflagged subscan {sid}.")

    def fit_baseline(self):
        """Fit a line to the marked intervals of the selected subscan and subtract it."""
        sid = self._need_selection()
        if not self.ranges:
            raise ValueError("Mark the source-free intervals first (box select)")
        self.session.fit_baseline(sid, self.channel, self.ranges)
        self.ranges = []
        self._after_edit(f"Baseline subtracted from subscan {sid}.")

    def align(self):
        """Align the subscans through the pixel (same direction) to the selected one."""
        sid = self._need_selection()
        direction = self._selected_subscan().direction
        same = [s.scan_id for s in self.subscans if s.direction == direction]
        self.session.align(same, self.channel, reference=sid)
        self._after_edit(f"Aligned {len(same) - 1} subscans to subscan {sid}.")

    def undo(self):
        """Cancel the last edit."""
        self.session.undo()
        self._after_edit("Last edit cancelled.")

    def apply(self):
        """Recompute the images with the current edits."""
        self.scanset.calculate_images()
        self.update_map()
        self.status.text = "Map updated."

    def save(self):
        """Write the edits to the scan files in the product directory."""
        n = self.session.n_unsaved
        self.session.save()
        self.update_edit_list()
        self.status.text = f"Saved {n} edits to the product directory."

    def set_map_kind(self, kind):
        """Show the image (``"Image"``) or the outlier map (``"Outliers"``)."""
        self.map_kind = kind
        self.update_map()

    def done(self):
        """Close the page and let the calling program continue."""
        if self.on_done is None:
            self.status.text = "Nothing to close here: just continue."
            return
        self.status.text = "Closed: the program continues. You can close this page."
        self.on_done()

    def __call__(self, doc):
        """Add the page to a Bokeh document (used by the Bokeh server and in notebooks)."""
        doc.title = "SDTimage interactive filter"
        doc.add_root(self.layout)
        doc.js_on_event(
            DocumentReady,
            CustomJS(
                args=dict(flag=self.unsaved_flag),
                code="""
window.addEventListener("beforeunload", (event) => {
    if (flag.text !== "0") { event.preventDefault(); event.returnValue = ""; }
});
""",
            ),
        )


def ssh_hint(port, host=None, user=None):
    """Command to forward ``port`` of this machine to the same port of a local computer.

    Examples
    --------
    >>> ssh_hint(5006, host="calc1", user="obs")
    'ssh -N -L 5006:localhost:5006 obs@calc1'
    """
    import getpass

    host = socket.gethostname() if host is None else host
    user = getpass.getuser() if user is None else user
    return f"ssh -N -L {port}:localhost:{port} {user}@{host}"


def _make_server(scanset, port):
    """Bokeh server of the filter page, and the edit session shared by all its pages."""
    from bokeh.server.server import Server

    session = FilterSession(scanset)
    server = None

    def stop():
        # Short delay: the page first shows that the program continues
        server.io_loop.call_later(0.5, server.io_loop.stop)

    def make_page(doc):
        FilterApp(scanset, session=session, on_done=stop)(doc)

    # Only reachable from this machine (or through an ssh tunnel)
    server = Server({"/": make_page}, address="localhost", port=port, num_procs=1)
    return server, session


def serve(scanset, port=5006, open_browser=None):
    """Serve the filter page and wait until the user presses "Close and continue".

    Each browser page gets its own copy of the widgets, but all share one
    :class:`srttools.filter_session.FilterSession`: reloading the page keeps
    the edits. The edited scan set can be used by the caller afterwards,
    whether or not the edits were saved to disk.

    Parameters
    ----------
    scanset : :class:`srttools.imager.ScanSet`
    port : int
        Port of the web server (0: any free port).
    open_browser : bool, optional
        Open the page in a local browser. By default, only when not logged in
        through ssh.

    Returns
    -------
    session : :class:`srttools.filter_session.FilterSession`
    """
    server, session = _make_server(scanset, port)
    server.start()
    url = f"http://localhost:{server.port}/"
    remote = "SSH_CONNECTION" in os.environ
    print(f"\nInteractive filter: open {url} in a browser.")
    if remote:
        print(
            "You are logged in through ssh. On your own computer, run\n\n"
            f"    {ssh_hint(server.port)}\n\nthen open {url} there."
        )
    print('Press "Close and continue" on the page (or Ctrl-C here) when done.\n')
    if open_browser if open_browser is not None else not remote:
        server.io_loop.add_callback(server.show, "/")
    try:
        server.io_loop.start()
    except KeyboardInterrupt:
        pass
    finally:
        server.stop()
    if session.n_unsaved > 0:
        print(
            f"{session.n_unsaved} edits not saved to disk: they are used for the "
            "images of this run only."
        )
    return session
