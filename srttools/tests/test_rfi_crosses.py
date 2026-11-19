import numpy as np
import pytest

from astropy.table import Table
from srttools.imager import ScanSet, find_outlier_subscans

NSUB = 40
NSAMP = 400
CHANS = ["Feed0_LCP", "Feed0_RCP"]


def _cross_arrays(seed=0, amplitude=0.0, burst_subscans=(), burst_peak=20.0):
    """Cross scans on 43 pixels, alternating arms, white noise of unit rms."""
    rng = np.random.default_rng(seed)
    along = np.linspace(1.5, 41.5, NSAMP)
    counts, coord, scan_id, arm = [], [], [], []
    t = np.arange(NSAMP)
    for i in range(NSUB):
        y = amplitude * np.exp(-((along - 21.5) ** 2) / (2 * 1.67**2))
        y = y + rng.normal(0, 1, NSAMP)
        if i in burst_subscans:
            y = y + burst_peak * np.exp(-((t - 150) ** 2) / (2 * 40**2))
        counts.append(y)
        coord.append(along)
        scan_id.append(np.full(NSAMP, i))
        arm.append(np.full(NSAMP, i % 2 == 0))
    return [np.concatenate(a) for a in (counts, coord, scan_id, arm)]


def test_find_outlier_subscans_flags_bursts():
    """Subscans with a long RFI burst are flagged, and only those."""
    counts, coord, scan_id, arm = _cross_arrays(burst_subscans=(3, 4, 17))
    bad = find_outlier_subscans(counts, coord, scan_id, arm, nsigma=5)
    assert sorted(bad) == [3, 4, 17]


@pytest.mark.parametrize("amplitude", [0.0, 100.0])
def test_find_outlier_subscans_clean(amplitude):
    """Clean data, with or without a bright source in every subscan, have no outliers."""
    counts, coord, scan_id, arm = _cross_arrays(amplitude=amplitude)
    assert len(find_outlier_subscans(counts, coord, scan_id, arm, nsigma=5)) == 0


def _scanset(**kwargs):
    counts, coord, scan_id, arm = _cross_arrays(**kwargs)
    t = Table()
    t["time"] = np.arange(counts.size) * 0.02
    t["x"] = np.where(arm, coord, 21.5)[:, np.newaxis]
    t["y"] = np.where(arm, 21.5, coord)[:, np.newaxis]
    t["direction"] = arm
    t["Scan_id"] = scan_id
    t["el"] = np.full((counts.size, 1), np.radians(45))
    for ch in CHANS:
        t[ch] = counts
        t[ch + "-filt"] = np.ones(counts.size, dtype=bool)
    t["x"].meta["frame"] = "icrs"
    t.meta["npix"] = np.array([43, 43])
    return ScanSet(t)


def test_calculate_avg_cross_rejects_rfi():
    """calculate_avg_cross flags RFI subscans with a warning, unless disabled."""
    scanset = _scanset(burst_subscans=(3,))
    with pytest.warns(UserWarning, match="RFI"):
        scanset.calculate_avg_cross(rfi_nsigma=5)
    for ch in CHANS:
        assert not np.any(scanset[ch + "-filt"][scanset["Scan_id"] == 3])
        assert np.all(scanset[ch + "-filt"][scanset["Scan_id"] != 3])

    scanset = _scanset(burst_subscans=(3,))
    scanset.calculate_avg_cross(rfi_nsigma=None)
    assert np.all(scanset[CHANS[0] + "-filt"])
