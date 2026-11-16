"""Interactive filtering of a ScanSet, kept in memory.

A :class:`FilterSession` holds the edits made to the subscans of a
:class:`srttools.imager.ScanSet` (zapped intervals, flags, baseline
corrections). Each edit changes the ScanSet columns immediately, so that new
images can be computed at once, but nothing is written to disk until
:meth:`FilterSession.save` is called, and then only into a product
directory (the one of the configuration, or any other except the raw data
directory). Edits can be undone one by one.

This module has no graphics: front ends (e.g. the web page) call its methods.
"""

import copy
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .fit import align as fit_align
from .fit import linear_fit, linear_fun
from .io import get_channel_feed, root_name

__all__ = ["FilterSession", "Subscan"]


@dataclass
class Subscan:
    """The samples of one subscan in one channel, as seen by a front end."""

    scan_id: int
    name: str
    direction: str  # "ra" or "dec": the coordinate changing along the subscan
    coord: np.ndarray  # coordinate along the subscan, in degrees
    values: np.ndarray
    good: np.ndarray  # False where the sample is masked


def _inside(coord, ranges):
    """True where ``coord`` falls in any of the (lo, hi) ``ranges`` (borders in any order)."""
    inside = np.zeros(len(coord), dtype=bool)
    for r in ranges:
        lo, hi = sorted(r)
        inside |= (coord >= lo) & (coord <= hi)
    return inside


class FilterSession:
    """Edits to the subscans of a ScanSet, kept in memory until saved.

    Parameters
    ----------
    scanset : :class:`srttools.imager.ScanSet`
        The scan set to edit. Its channel columns and their ``-filt`` masks
        are modified in place by each edit.

    Attributes
    ----------
    edits : list of dict
        The edits, in the order they were made.
    """

    def __init__(self, scanset):
        self.scanset = scanset
        self.edits = []
        self._n_saved = 0
        self._columns = list(scanset.chan_columns) + [f"{c}-filt" for c in scanset.chan_columns]
        self._original = {col: np.array(scanset[col]) for col in self._columns}
        self._scan_masks = {}

    @property
    def n_unsaved(self):
        """Number of edits not yet written to disk."""
        return len(self.edits) - self._n_saved

    def _mask(self, scan_id):
        if scan_id not in self._scan_masks:
            self._scan_masks[scan_id] = np.asarray(self.scanset["Scan_id"] == scan_id)
        return self._scan_masks[scan_id]

    def _direction(self, scan_id):
        mask = self._mask(scan_id)
        if "direction" in self.scanset.colnames:
            return "ra" if self.scanset["direction"][mask][0] else "dec"
        ra, dec = self.scanset["ra"][mask][:, 0], self.scanset["dec"][mask][:, 0]
        return "ra" if np.ptp(ra) > np.ptp(dec) else "dec"

    def _coord(self, scan_id, channel):
        """Coordinate along the subscan, in degrees, for the feed of ``channel``."""
        feed = get_channel_feed(channel)
        col = self.scanset[self._direction(scan_id)]
        return np.degrees(np.asarray(col[self._mask(scan_id)][:, feed], dtype=float))

    def subscan(self, scan_id, channel):
        """Return a :class:`Subscan` with the current data of one subscan."""
        mask = self._mask(scan_id)
        return Subscan(
            scan_id=int(scan_id),
            name=self.scanset.scan_list[scan_id],
            direction=self._direction(scan_id),
            coord=self._coord(scan_id, channel),
            values=np.array(self.scanset[channel][mask], dtype=float),
            good=np.array(self.scanset[f"{channel}-filt"][mask], dtype=bool),
        )

    def subscans_through_pixel(self, x, y, channel):
        """List the subscans passing within one pixel of pixel ``(x, y)``."""
        feed = get_channel_feed(channel)
        near = (np.abs(self.scanset["x"][:, feed] - x) < 1) & (
            np.abs(self.scanset["y"][:, feed] - y) < 1
        )
        scan_ids = sorted(set(np.asarray(self.scanset["Scan_id"][near]).tolist()))
        return [self.subscan(sid, channel) for sid in scan_ids]

    def scan_summary(self, channel):
        """One row per subscan, for a list of scans in a front end.

        Returns
        -------
        summary : dict of lists
            ``scan_id``, ``name`` (file name without directory), ``direction``
            (``"ra"`` or ``"dec"``), ``masked`` (fraction of masked samples of
            ``channel``) and ``flagged`` (True if all samples of all channels
            are masked), in order of ``scan_id``.
        """
        ids, inverse = np.unique(np.asarray(self.scanset["Scan_id"]), return_inverse=True)
        n = np.bincount(inverse)
        masked = np.bincount(
            inverse, weights=~np.asarray(self.scanset[f"{channel}-filt"], dtype=bool)
        )
        any_good = np.zeros(len(ids))
        for c in self.scanset.chan_columns:
            any_good += np.bincount(
                inverse, weights=np.asarray(self.scanset[f"{c}-filt"], dtype=bool)
            )
        return dict(
            scan_id=[int(i) for i in ids],
            name=[os.path.basename(self.scanset.scan_list[i]) for i in ids],
            direction=[self._direction(i) for i in ids],
            masked=(masked / n).tolist(),
            flagged=(any_good == 0).tolist(),
        )

    # ----- Edits ----------------------------------------------------------

    def _add(self, edit):
        self.edits.append(edit)
        self._apply(edit)

    def zap(self, scan_id, channel, lo, hi):
        """Mask the samples of one channel of a subscan between ``lo`` and ``hi`` (degrees)."""
        self._add(dict(kind="zap", scan_id=int(scan_id), channel=channel, range=(lo, hi)))

    def flag(self, scan_id):
        """Mask the whole subscan, in all channels."""
        self._add(dict(kind="flag", scan_id=int(scan_id), channel=None))

    def unflag(self, scan_id):
        """Remove all masks of the subscan, in all channels (zaps included)."""
        self._add(dict(kind="unflag", scan_id=int(scan_id), channel=None))

    def fit_baseline(self, scan_id, channel, ranges):
        """Fit a line to the unmasked samples inside ``ranges`` and subtract it.

        Parameters
        ----------
        scan_id : int
        channel : str
        ranges : list of (lo, hi)
            Intervals of the coordinate along the subscan (degrees) with no
            source emission.
        """
        s = self.subscan(scan_id, channel)
        use = _inside(s.coord, ranges) & s.good
        if np.count_nonzero(use) < 2:
            raise ValueError("Select at least two unmasked samples for the baseline fit")
        q, m = linear_fit(s.coord[use], s.values[use], [np.median(s.values[use]), 0])
        self._add(dict(kind="baseline", scan_id=int(scan_id), channel=channel, pars=(q, m)))

    def align(self, scan_ids, channel, reference):
        """Subtract from each subscan the line that best aligns it to ``reference``.

        All the corrections are one edit (one undo).
        """
        others = [int(sid) for sid in scan_ids if sid != reference]
        xs, ys, kept = [], [], []
        for sid in [reference] + others:
            s = self.subscan(sid, channel)
            if not np.any(s.good):
                continue
            xs.append(s.coord[s.good])
            ys.append(s.values[s.good])
            kept.append(sid)
        if len(kept) < 2:
            return
        qs, ms = fit_align(xs, ys)
        pars = {sid: (q, m) for sid, q, m in zip(kept[1:], qs, ms)}
        self._add(dict(kind="align", scan_id=None, channel=channel, pars=pars))

    def undo(self):
        """Cancel the last edit (also when already saved; save again to update the files)."""
        if not self.edits:
            return
        self.edits.pop()
        self._n_saved = min(self._n_saved, len(self.edits))
        self._replay()

    def reset(self):
        """Cancel all edits."""
        self.edits = []
        self._n_saved = 0
        self._replay()

    # ----- Applying edits to the columns ----------------------------------

    def _replay(self):
        for col in self._columns:
            self.scanset[col][:] = self._original[col]
        for edit in self.edits:
            self._apply(edit)

    def _subtract_line(self, scan_id, channel, q, m):
        mask = self._mask(scan_id)
        coord = self._coord(scan_id, channel)
        self.scanset[channel][mask] = self.scanset[channel][mask] - linear_fun(coord, q, m)

    def _apply(self, edit):
        kind, sid, ch = edit["kind"], edit["scan_id"], edit["channel"]
        if kind == "zap":
            mask = self._mask(sid)
            zapped = _inside(self._coord(sid, ch), [edit["range"]])
            self.scanset[f"{ch}-filt"][mask] = self.scanset[f"{ch}-filt"][mask] & ~zapped
        elif kind in ("flag", "unflag"):
            mask = self._mask(sid)
            for c in self.scanset.chan_columns:
                self.scanset[f"{c}-filt"][mask] = kind == "unflag"
        elif kind == "baseline":
            self._subtract_line(sid, ch, *edit["pars"])
        elif kind == "align":
            for other, (q, m) in edit["pars"].items():
                self._subtract_line(other, ch, q, m)
        else:  # pragma: no cover
            raise ValueError(f"Unknown edit: {kind}")

    # ----- Saving ---------------------------------------------------------

    def _edited_scans(self):
        sids = set()
        for edit in self.edits:
            if edit["kind"] == "align":
                sids |= set(edit["pars"])
            else:
                sids.add(edit["scan_id"])
        return sorted(sids)

    @property
    def config_productdir(self):
        """Product directory of the configuration file, or None without a configuration."""
        from .read_config import read_config

        productdir = self.scanset.meta.get("productdir")
        config_file = self.scanset.meta.get("config_file")
        if productdir is None and config_file is not None:
            productdir = read_config(config_file).get("productdir")
        return None if productdir is None else os.path.abspath(productdir)

    @property
    def default_productdir(self):
        """Where :meth:`save` writes if no directory is given: never empty.

        The product directory of the configuration file or, without one,
        ``srttools_products`` in the current directory.
        """
        from .read_config import DEFAULT_PRODUCTDIR

        return self.config_productdir or os.path.abspath(DEFAULT_PRODUCTDIR)

    def _outfile(self, sname, productdir):
        """File where the scan ``sname`` is saved under ``productdir``."""
        from .scan import product_path_from_file_name

        rootdir, fn = product_path_from_file_name(
            sname,
            workdir=self.scanset.meta.get("workdir") or os.curdir,
            productdir=productdir,
            basename_only=bool(self.scanset.meta.get("basename_only", False)),
        )
        return root_name(os.path.join(rootdir, fn)) + ".hdf5"

    def check_productdir(self, productdir):
        """Return the absolute ``productdir``, or raise if edits cannot be saved there.

        Raises
        ------
        ValueError
            If the name is empty, if the directory cannot be created or
            written to, or if some scan would be saved in the directory of
            its raw data file: the raw data directory is never written to.
        """
        if productdir is None or not str(productdir).strip():
            raise ValueError("Give a directory where to save the edits")
        productdir = os.path.abspath(Path(str(productdir).strip()).expanduser())
        parent = productdir  # the deepest part of the path that already exists
        while not os.path.exists(parent) and os.path.dirname(parent) != parent:
            parent = os.path.dirname(parent)
        if not os.path.isdir(parent) or not os.access(parent, os.W_OK):
            raise ValueError(f"{productdir} cannot be written to (check {parent})")
        for sname in self.scanset.scan_list:
            rawdir = os.path.dirname(os.path.abspath(sname))
            if os.path.dirname(self._outfile(sname, productdir)) == rawdir:
                raise ValueError(
                    f"{productdir} would put the products in the raw data directory "
                    f"{rawdir}, which is never written to: choose another directory"
                )
        return productdir

    def config_hint(self, productdir):
        """Explain how a later run finds files saved in ``productdir``; empty if it already does."""
        productdir = os.path.abspath(Path(str(productdir).strip()).expanduser())
        if productdir == self.config_productdir:
            return ""
        return (
            "This is not the productdir of the config file: to use these files in "
            f"later runs, put `productdir : {productdir}` in its [local] section."
        )

    def save(self, productdir=None):
        """Write the edited subscans to their files in a product directory.

        Parameters
        ----------
        productdir : str, optional
            Where to save. Default: :attr:`default_productdir`, so there is
            always somewhere to save.

        Returns
        -------
        productdir : str
            The directory the files were written to.

        Raises
        ------
        ValueError
            If ``productdir`` would put the files next to the raw data.
        """
        from .scan import Scan

        productdir = self.check_productdir(productdir or self.default_productdir)
        config_file = self.scanset.meta.get("config_file")

        # Flags of the last flag/unflag edit of each scan
        flags = {}
        for edit in self.edits:
            if edit["kind"] in ("flag", "unflag"):
                flags[edit["scan_id"]] = edit["kind"] == "flag"

        for sid in self._edited_scans():
            sname = self.scanset.scan_list[sid]
            s = Scan(sname, config_file=config_file, nosave=True)
            outfile = self._outfile(sname, productdir)
            mask = self._mask(sid)
            if len(s) != np.count_nonzero(mask):
                raise ValueError(f"{sname} has changed since the scan set was created")
            for col in self._columns:
                s[col] = copy.deepcopy(np.asarray(self.scanset[col][mask]))
            if sid in flags:
                s.meta["FLAG"] = flags[sid]
            s.meta["productdir"] = productdir
            s.save(outfile)
        self._n_saved = len(self.edits)
        return productdir
