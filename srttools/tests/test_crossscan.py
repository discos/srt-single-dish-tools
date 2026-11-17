import glob

import numpy as np
import pytest

import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from srttools.crossscan import (
    _shape_warning,
    binned_profile,
    combine_amplitudes,
    consistency,
    fit_subscan,
    main_cross,
    measure_crosses,
    subscan_amplitudes,
    upper_limit,
)

BEAM_SIGMA = 1.67
NSAMP = 200


def _subscan(amplitude, rng, noise=0.0, center=0.0):
    """One subscan along +-20 pixels: Gaussian source + cubic baseline + white noise."""
    x = np.linspace(-20, 20, NSAMP)
    base = 3 + 0.1 * x - 0.01 * x**2 + 2e-4 * x**3
    y = base + amplitude * np.exp(-((x - center) ** 2) / (2 * BEAM_SIGMA**2))
    return x, y + rng.normal(0, noise, x.size) if noise > 0 else y


def _observation(amplitude, rng, nsub=20, noise=1.0):
    """Several subscans on two arms, concatenated as in a ScanSet."""
    xs, ys, ids, arms = [], [], [], []
    for i in range(nsub):
        x, y = _subscan(amplitude, rng, noise=noise)
        xs.append(x)
        ys.append(y)
        ids.append(np.full(x.size, i))
        arms.append(np.full(x.size, i % 2 == 0))
    return [np.concatenate(a) for a in (xs, ys, ids, arms)]


@pytest.mark.parametrize("joint", [False, True])
def test_fit_subscan_noiseless(joint):
    """The forced amplitude of a noiseless source on a curved baseline is exact."""
    x, y = _subscan(5.0, None)
    amp, resid = fit_subscan(x, y, BEAM_SIGMA, joint=joint)
    assert np.isclose(amp, 5.0, rtol=1e-3)
    assert np.allclose(resid, 5.0 * np.exp(-(x**2) / (2 * BEAM_SIGMA**2)), atol=5e-3)


def test_fit_subscan_default_is_joint():
    """The default is the joint fit of baseline and source, chosen by the Monte Carlo."""
    rng = np.random.default_rng(3)
    x, y = _subscan(5.0, rng, noise=1.0)
    assert fit_subscan(x, y, BEAM_SIGMA)[0] == fit_subscan(x, y, BEAM_SIGMA, joint=True)[0]


def test_fit_subscan_too_few_points():
    """A subscan with too few good samples gives NaN instead of raising."""
    amp, resid = fit_subscan(np.arange(5.0), np.zeros(5), BEAM_SIGMA)
    assert np.isnan(amp)
    assert np.all(np.isnan(resid))


def test_combine_amplitudes():
    """Mean and error from the scatter between subscans; the median is also available."""
    amps = np.array([1.0, 2.0, 3.0, 4.0, 100.0])
    flux, err = combine_amplitudes(amps, statistic="mean")
    assert np.isclose(flux, amps.mean())
    assert np.isclose(err, np.std(amps, ddof=1) / np.sqrt(5))
    flux, _ = combine_amplitudes(amps, statistic="median")
    assert flux == 3.0
    flux, _ = combine_amplitudes(np.array([1.0, np.nan, 3.0]))
    assert flux == 2.0


def test_upper_limit():
    """Upper limit is max(flux, 0) + n sigma, with the calibration error in quadrature."""
    assert upper_limit(-1.0, 1.0, nsigma=3) == 3.0
    assert upper_limit(2.0, 1.0, nsigma=3) == 5.0
    assert np.isclose(upper_limit(10.0, 1.0, nsigma=3, cal_rel_err=0.1), 10 + 3 * np.sqrt(2))


def test_consistency():
    """The difference between two estimates is given in units of its error."""
    assert np.isclose(consistency((1.0, 0.3), (2.0, 0.4)), -2.0)


@pytest.mark.parametrize("joint", [False, True])
def test_amplitude_error_coverage(joint):
    """With white noise, the error from the subscan scatter matches the real scatter."""
    rng = np.random.default_rng(1234)
    pulls = []
    for _ in range(100):
        x, y, ids, arm = _observation(1.0, rng)
        _, amps = subscan_amplitudes(x, y, ids, BEAM_SIGMA, joint=joint)
        flux, err = combine_amplitudes(amps)
        pulls.append((flux - 1.0) / err)
    assert abs(np.mean(pulls)) < 0.3
    assert 0.8 < np.std(pulls) < 1.2


def test_binned_profile():
    """The binned residual profile shows the source, with errors from the subscan scatter."""
    rng = np.random.default_rng(0)
    x, y, ids, arm = _observation(10.0, rng, noise=0.5)
    centers, prof, err = binned_profile(x, y, ids, BEAM_SIGMA, bin_width=1.0)
    i0 = np.argmin(np.abs(centers))
    assert np.isclose(prof[i0], 10.0 * np.exp(-(centers[i0] ** 2) / (2 * BEAM_SIGMA**2)), atol=1)
    assert np.all(err[np.isfinite(err)] > 0)
    assert np.abs(prof[np.abs(centers) > 10]).max() < 1


RA0, DEC0 = np.radians(248.26), np.radians(23.83)
FWHM_DEG = 0.0463
CHANS = ["Feed0_LCP", "Feed0_RCP"]


def _cross_table(amplitude, seed=0, nsub=40, src_offset_deg=0.0, src_width=1.0):
    """Cross scans as stored in a ScanSet: ra/dec in radians, one feed, two channels."""
    rng = np.random.default_rng(seed)
    sigma = src_width * FWHM_DEG / 2.3548
    off = np.linspace(-0.24, 0.24, NSAMP)
    cols = {k: [] for k in ["ra", "dec", "direction", "Scan_id", *CHANS]}
    for i in range(nsub):
        isra = i % 2 == 0
        dra = off / np.cos(DEC0) if isra else np.zeros(NSAMP)
        ddec = np.zeros(NSAMP) if isra else off
        cols["ra"].append(RA0 + np.radians(dra))
        cols["dec"].append(DEC0 + np.radians(ddec))
        cols["direction"].append(np.full(NSAMP, isra))
        cols["Scan_id"].append(np.full(NSAMP, i))
        src = amplitude * np.exp(-((off - src_offset_deg) ** 2) / (2 * sigma**2))
        for ch in CHANS:
            cols[ch].append(src + 0.5 * off + rng.normal(0, 1, NSAMP))
    t = Table()
    for k, v in cols.items():
        v = np.concatenate(v)
        t[k] = v[:, np.newaxis] if k in ["ra", "dec"] else v
    for ch in CHANS:
        t[ch + "-filt"] = np.ones(len(t), dtype=bool)
    return t


def test_measure_crosses():
    """Fluxes per channel, arm and in total are right; flagged samples are ignored."""
    t = _cross_table(3.0)
    bad = t["Scan_id"] == 5
    t["Feed0_LCP"][bad] += 1000
    t["Feed0_LCP-filt"][bad] = False
    res = measure_crosses(t, RA0, DEC0, FWHM_DEG, channels=CHANS, rfi_nsigma=None)
    summary = res["summary"]
    assert set(summary["label"]) == {
        "Feed0_LCP RA", "Feed0_LCP Dec", "Feed0_RCP RA", "Feed0_RCP Dec",
        "Feed0_LCP", "Feed0_RCP", "RA", "Dec", "TOTAL",
    }  # fmt: skip
    total = summary[summary["label"] == "TOTAL"][0]
    assert abs(total["flux"] - 3.0) < 4 * total["err"]
    assert total["nsub"] == 40
    assert summary[summary["label"] == "Feed0_LCP"][0]["nsub"] == 39
    assert total["detected"]
    assert abs(res["consistency"]["RA-Dec"]) < 4
    assert abs(res["consistency"]["Feed0_LCP-Feed0_RCP"]) < 4
    assert abs(res["consistency"]["joint-window"]) < 4
    centers, prof, err = res["profiles"]["TOTAL"]
    assert np.isclose(prof[np.argmin(np.abs(centers))], 3.0, atol=4 * total["err"])


def test_measure_crosses_extended_source():
    """For a source twice as wide as the beam, the joint and window fits disagree."""
    t = _cross_table(50.0, src_width=2.0)
    res = measure_crosses(t, RA0, DEC0, FWHM_DEG, channels=CHANS)
    assert res["consistency"]["joint-window"] < -5


def test_shape_warning():
    """The joint-window warning explains the likely cause, depending on the sign."""
    assert _shape_warning(1.0, 3.0) is None
    assert "underestimated" in _shape_warning(-4.0, 3.0)
    assert "--beam-fwhm" in _shape_warning(4.0, 3.0)


def test_measure_crosses_offset_position():
    """The flux is measured at the requested position, not at the scan center."""
    t = _cross_table(10.0, src_offset_deg=0.08)
    res = measure_crosses(t, RA0, DEC0 + np.radians(0.08), FWHM_DEG, channels=CHANS)
    dec_arm = res["summary"][res["summary"]["label"] == "Dec"][0]
    assert abs(dec_arm["flux"] - 10.0) < 4 * dec_arm["err"]


def _fake_dump(fname, **meta):
    """Write a dump file of a simulated cross scan, as SDTimage would."""
    t = _cross_table(5.0)
    t["x"] = np.zeros((len(t), 1))
    t["y"] = np.zeros((len(t), 1))
    t.meta.update(RA=RA0 * u.rad, Dec=DEC0 * u.rad, SOURCE="Fake", **meta)
    for ch in CHANS:
        t[ch].meta.update(frequency=6200 * u.MHz, bandwidth=1500 * u.MHz)
    t.write(fname, serialize_meta=True)
    return fname


@pytest.mark.usefixtures("no_pyplot_figures")
def test_sdtcross_command_line(tmp_path):
    """SDTcross reads a dump file and writes the FITS, ECSV and plot files."""
    fname = _fake_dump(str(tmp_path / "fake_dump.hdf5"))
    with pytest.warns(UserWarning, match="No calibration"):
        main_cross([fname, "--beam-fwhm", str(FWHM_DEG * 60), "--rfi-nsigma", "0"])
    root = str(tmp_path / "fake")
    header = fits.getheader(root + "_cross.fits")
    assert header["SOURCE"] == "Fake"
    assert header["DETECTED"]
    assert abs(header["CONS_SHP"]) < 4
    assert abs(header["FLUX"] - 5.0) < 4 * header["FLUXERR"]
    assert len(fits.getdata(root + "_cross.fits", "SUBSCANS")) == 40
    summary = Table.read(root + "_cross.ecsv")
    assert "TOTAL" in summary["label"]
    for suffix in ["_cross.png", "_cross_total.png"]:
        assert (tmp_path / ("fake" + suffix)).exists()


def test_sdtcross_writes_to_productdir_of_dump(tmp_path):
    """Without -o, the results go to the productdir stored in the dump, not next to it."""
    (tmp_path / "intermediate").mkdir()
    fname = _fake_dump(
        str(tmp_path / "intermediate" / "fake_dump.hdf5"), productdir=str(tmp_path / "prod")
    )
    with pytest.warns(UserWarning, match="No calibration"):
        main_cross([fname, "--beam-fwhm", str(FWHM_DEG * 60), "--rfi-nsigma", "0", "--noplot"])
    assert (tmp_path / "prod" / "fake_cross.fits").exists()
    assert not glob.glob(str(tmp_path / "intermediate" / "*_cross*"))
