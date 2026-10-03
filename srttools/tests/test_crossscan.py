import numpy as np
import pytest

from srttools.crossscan import (
    binned_profile,
    combine_amplitudes,
    consistency,
    fit_subscan,
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
