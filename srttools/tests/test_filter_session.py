"""Tests of the in-memory filter session (no graphics)."""

import glob
import os
import warnings

import numpy as np
import pytest

from astropy.table import Table
from srttools.filter_session import FilterSession
from srttools.imager import ScanSet
from srttools.scan import Scan

MEDICINA = os.path.join(os.path.dirname(__file__), "data", "medicina")
NPTS = 50


def _fake_scanset(slope=0.0, offsets=(0.0, 0.0, 0.0, 0.0)):
    """Two Dec scans (ids 0, 1) and two RA scans (ids 2, 3) crossing at pixel (5, 5).

    Each scan is a Gaussian source on top of an offset plus a linear slope (per degree).
    """
    deg = np.linspace(-1, 1, NPTS)
    rows = []
    for sid, (direction, pos) in enumerate([(0, 4.6), (0, 5.4), (1, 4.6), (1, 5.4)]):
        along = 5 + 5 * deg  # pixel coordinate along the scan
        if direction:  # RA scan: x varies
            x, y, ra, dec = along, np.full(NPTS, pos), deg, np.full(NPTS, (pos - 5) / 5)
        else:
            x, y, ra, dec = np.full(NPTS, pos), along, np.full(NPTS, (pos - 5) / 5), deg
        values = np.exp(-(deg**2) / 0.02) + offsets[sid] + slope * deg
        for i in range(NPTS):
            rows.append((ra[i], dec[i], x[i], y[i], values[i], 2 * values[i], sid, direction))
    names = ["ra", "dec", "x", "y", "Feed0_RCP", "Feed0_LCP", "Scan_id", "direction"]
    t = Table(rows=rows, names=names)
    for col in ["ra", "dec"]:
        t[col] = np.radians(np.array(t[col]))[:, np.newaxis]
    for col in ["x", "y"]:
        t[col] = np.array(t[col])[:, np.newaxis]
    t["direction"] = np.asarray(t["direction"], dtype=bool)
    for ch in ["Feed0_RCP", "Feed0_LCP"]:
        t[f"{ch}-filt"] = np.ones(len(t), dtype=bool)
    t.meta["scan_list"] = [f"scan{i}.fits" for i in range(4)]
    return ScanSet(t)


def test_subscans_through_pixel():
    """All four scans pass near pixel (5, 5); each comes with its direction and coordinate in deg."""
    session = FilterSession(_fake_scanset())
    subscans = session.subscans_through_pixel(5, 5, "Feed0_RCP")
    assert sorted(s.scan_id for s in subscans) == [0, 1, 2, 3]
    for s in subscans:
        assert s.direction == ("ra" if s.scan_id >= 2 else "dec")
        assert np.allclose(s.coord, np.linspace(-1, 1, NPTS))
        assert np.all(s.good)
    assert session.subscans_through_pixel(0, 0, "Feed0_RCP") == []


def test_zap_masks_one_channel_inside_interval():
    """A zap masks only the samples inside the interval, in one channel of one scan."""
    scanset = _fake_scanset()
    session = FilterSession(scanset)
    session.zap(1, "Feed0_RCP", 0.5, 0.2)  # order of the borders does not matter
    deg = np.linspace(-1, 1, NPTS)
    mask = np.asarray(scanset["Scan_id"] == 1)
    expected = ~((deg >= 0.2) & (deg <= 0.5))
    assert np.all(scanset["Feed0_RCP-filt"][mask] == expected)
    assert np.all(scanset["Feed0_LCP-filt"])
    assert np.all(scanset["Feed0_RCP-filt"][~mask])


def test_flag_and_unflag():
    """Flag masks all channels of a scan; unflag clears all its masks, including old zaps."""
    scanset = _fake_scanset()
    session = FilterSession(scanset)
    mask = np.asarray(scanset["Scan_id"] == 2)
    session.zap(2, "Feed0_RCP", -0.1, 0.1)
    session.flag(2)
    assert not np.any(scanset["Feed0_RCP-filt"][mask])
    assert not np.any(scanset["Feed0_LCP-filt"][mask])
    session.unflag(2)
    assert np.all(scanset["Feed0_RCP-filt"]) and np.all(scanset["Feed0_LCP-filt"])


def test_fit_baseline_removes_slope():
    """A linear baseline fitted outside the source removes offset and slope from one channel."""
    scanset = _fake_scanset(slope=0.5, offsets=(0.3, 0, 0, 0))
    session = FilterSession(scanset)
    lcp_before = np.array(scanset["Feed0_LCP"])
    session.fit_baseline(0, "Feed0_RCP", [(-1, -0.5), (0.5, 1)])
    mask = np.asarray(scanset["Scan_id"] == 0)
    deg = np.linspace(-1, 1, NPTS)
    assert np.allclose(scanset["Feed0_RCP"][mask], np.exp(-(deg**2) / 0.02), atol=1e-3)
    assert np.all(scanset["Feed0_LCP"] == lcp_before)


def test_align_to_reference():
    """Aligning removes the offsets of the other scans with respect to the reference."""
    scanset = _fake_scanset(offsets=(0.0, 0.4, -0.2, 0.1))
    session = FilterSession(scanset)
    session.align([0, 1], "Feed0_RCP", reference=0)
    m0 = np.asarray(scanset["Scan_id"] == 0)
    m1 = np.asarray(scanset["Scan_id"] == 1)
    assert np.allclose(scanset["Feed0_RCP"][m1], scanset["Feed0_RCP"][m0], atol=1e-3)


def test_undo_and_reset():
    """Undo restores the state before the last edit, reset the original one; edits are counted."""
    scanset = _fake_scanset(slope=0.5)
    original = scanset.copy()
    session = FilterSession(scanset)
    session.zap(0, "Feed0_RCP", -0.1, 0.1)
    after_zap = scanset.copy()
    session.fit_baseline(0, "Feed0_RCP", [(-1, -0.5), (0.5, 1)])
    assert session.n_unsaved == 2
    session.undo()
    assert session.n_unsaved == 1
    for col in ["Feed0_RCP", "Feed0_RCP-filt"]:
        assert np.all(scanset[col] == after_zap[col])
    session.reset()
    assert session.n_unsaved == 0
    for col in ["Feed0_RCP", "Feed0_RCP-filt"]:
        assert np.all(scanset[col] == original[col])


def test_save_needs_productdir():
    """Without a productdir, save refuses to write anything."""
    session = FilterSession(_fake_scanset())
    session.flag(0)
    with pytest.raises(ValueError, match="productdir"):
        session.save()


def test_save_writes_only_to_productdir(tmp_path):
    """Saved edits go to the scan files in productdir, and nothing is written next to the data."""
    config = str(tmp_path / "casa.ini")
    with open(config, "w") as fobj:
        print(
            f"[local]\nworkdir = {tmp_path}\ndatadir = {MEDICINA}\n"
            f"productdir = {tmp_path / 'products'}\n\n[analysis]\n"
            "list_of_directories =\n    20221213-162927-22-22-cas_a/\npixel_size = 1.5\n",
            file=fobj,
        )
    data_before = sorted(glob.glob(os.path.join(MEDICINA, "**"), recursive=True))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        scanset = ScanSet(config)
        session = FilterSession(scanset)
        sid = int(scanset["Scan_id"][0])
        coord = session.subscan(sid, "Feed0_RCP").coord
        session.zap(sid, "Feed0_RCP", coord[0], coord[10])
        session.save()
        assert session.n_unsaved == 0
        saved = Scan(scanset.scan_list[sid], config_file=config)
    assert not np.any(saved["Feed0_RCP-filt"][:10])
    assert np.all(saved["Feed0_LCP-filt"] == scanset["Feed0_LCP-filt"][scanset["Scan_id"] == sid])
    assert sorted(glob.glob(os.path.join(MEDICINA, "**"), recursive=True)) == data_before
