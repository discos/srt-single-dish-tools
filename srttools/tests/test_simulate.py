import os
import shutil
import subprocess as sp

import numpy as np
import pytest

from astropy.io import fits
from srttools.io import mkdir_p
from srttools.simulate import (
    _default_map_shape,
    _offset_arm_shape,
    _red_noise,
    main_simulate,
    save_scan,
    sim_crossscans,
    sim_position_switching,
    simulate_map,
    simulate_scan,
    simulate_sun,
)


class TestSimulate:
    @classmethod
    def setup_class(cls):
        cls.outdir = os.path.join("sim")
        cls.emptydir = os.path.join("sim", "empty")
        cls.pswdir = os.path.join("sim", "psw")
        for d in [cls.emptydir, cls.pswdir]:
            mkdir_p(d)

    def test_script_is_installed(self):
        sp.check_call("SDTfake -h".split(" "))

    def test_sim_scan(self):
        """Test the simulation of a single scan."""
        times, position, shape = simulate_scan()
        save_scan(
            times,
            position,
            np.zeros_like(position),
            {"Ch0": shape, "Ch1": shape},
            os.path.join(self.outdir, "output.fits"),
        )

    def test_sim_scan_other_columns(self):
        """Test the simulation of a single scan."""
        times, position, shape = simulate_scan()
        save_scan(
            times,
            position,
            np.zeros_like(position),
            {"Ch0": shape, "Ch1": shape},
            os.path.join(self.outdir, "output.fits"),
            scan_type="Any",
            other_columns={"Pippo": np.zeros_like(shape)},
        )

    def test_position_switching_legacy(self):
        """Test the simulation of an empty map."""
        sim_position_switching(self.pswdir, legacy_cal_format=True, nbin=128)

    def test_sim_crossscans(self):
        sim_crossscans(1, self.pswdir)

    def test_sim_crossscans_noise_options(self):
        """Smoke test: cross scans with seeded 1/f noise and a pointing offset."""
        sim_crossscans(
            1,
            self.pswdir,
            noise_amplitude=0.5,
            rng=1234,
            fknee=1.0,
            pol_correlation=0.5,
            pointing_offset=(0.01, -0.01),
        )

    def test_sim_crossscans_source_centered_in_both_directions(self, tmp_path):
        """The source is at the nominal position also in the crosses scanned backwards."""
        sim_crossscans(2, str(tmp_path), noise_amplitude=0.0, rng=1)
        for fname, col, ref in [("Ra", "raj2000", 185.0), ("Dec", "decj2000", 75.0)]:
            for i in range(2):
                data = fits.getdata(tmp_path / f"{i}_{fname}.fits", "DATA TABLE")
                pos = np.degrees(data[col]) - ref
                if col == "raj2000":
                    pos *= np.cos(np.radians(75.0))
                step = np.median(np.abs(np.diff(pos)))
                y = data["Ch0"]
                good = y > 0.1 * y.max()
                centroid = np.sum(pos[good] * y[good]) / np.sum(y[good])
                assert np.abs(centroid) < 0.2 * step, (fname, i, centroid / step)

    def test_sim_map_empty(self):
        """Test the simulation of an empty map."""
        out_ra, _ = simulate_map(width_ra=2, width_dec=2.0, outdir=self.emptydir)
        probe = os.path.join(out_ra, "Ra0.fits")
        assert os.path.exists(probe)
        with fits.open(probe) as hdul:
            assert hdul[0].header["Declination Offset"] != 0.0
        os.unlink(probe)

    def test_sim_sun_map_small(self):
        """Test the simulation of an empty map."""
        out_ra, _ = simulate_sun(length_ra=2, length_dec=2.0, outdir=self.emptydir)
        probe = os.path.join(out_ra, "Ra0.fits")
        assert os.path.exists(probe)
        os.unlink(probe)

    def test_sim_map_small(self):
        """Test the simulation of an empty map."""
        out_ra, _ = simulate_map(length_ra=2, length_dec=2.0, outdir=self.emptydir)
        probe = os.path.join(out_ra, "Ra0.fits")
        assert os.path.exists(probe)
        os.unlink(probe)

    @pytest.mark.usefixtures("no_pyplot_figures")
    def test_sim_map_empty_messy(self):
        """Test the simulation of an empty map."""
        simulate_map(
            width_ra=2,
            width_dec=2.0,
            outdir=self.emptydir,
            baseline="messy",
            nbin=100,
            debug=True,
        )

    def test_sim_map_empty_slope(self):
        """Test the simulation of an empty map."""
        simulate_map(width_ra=2, width_dec=2.0, outdir=self.emptydir, baseline="slope")

    def test_raises_wrong_map_shape(self):
        with pytest.raises(ValueError):
            _default_map_shape(np.zeros((3, 4)), np.ones((3, 6)))

    def test_use_script(self):
        main_simulate("--no-cal -s 0.005 -g 10 10 1 1 -o sima -b 1".split(" "))
        shutil.rmtree("sima")
        main_simulate("--no-cal -g 10 10 1 1 -o simb -b messy".split(" "))
        shutil.rmtree("simb")
        main_simulate("--no-cal -g 10 10 1 1 -o simc -b slope".split(" "))
        shutil.rmtree("simc")
        main_simulate("--integration-time .02 -g 10 10 1 1 -o simd".split(" "))
        shutil.rmtree("simd")
        main_simulate("--no-cal --scan-speed 3. -g 10 10 1 1 " "-o sime".split(" "))
        shutil.rmtree("sime")

    def test_use_wrong_baseline(self):
        with pytest.raises(ValueError):
            main_simulate("-b qwerty -g 10 10 1 1 -o sime".split(" "))

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.outdir)


def test_simulate_scan_seed_reproducible():
    """The same seed gives identical scans, including baseline and 1/f noise."""
    kw = dict(baseline="messy", fknee=1.0, noise_amplitude=1.0)
    _, _, a = simulate_scan(rng=42, **kw)
    _, _, b = simulate_scan(rng=42, **kw)
    _, _, c = simulate_scan(rng=43, **kw)
    assert np.array_equal(a, b)
    assert not np.array_equal(a, c)


def test_simulate_scan_white_noise_level():
    """With a zero source and flat baseline, the scatter equals noise_amplitude."""
    _, _, scan = simulate_scan(
        shape=lambda x: np.zeros_like(x), noise_amplitude=3.0, nsamples=20000, rng=0
    )
    assert np.isclose(np.std(scan), 3.0, rtol=0.03)


def test_red_noise_spectrum():
    """1/f noise power at f << fknee exceeds that at f >> fknee by about fknee/f."""
    dt, n, fknee = 0.01, 4096, 1.0
    noise = _red_noise(n, dt, fknee, noise_amplitude=1.0, ntrials=200, rng=1)[:, 0, :]
    power = np.mean(np.abs(np.fft.rfft(noise, axis=1)) ** 2, axis=0) / n
    freqs = np.fft.rfftfreq(n, dt)
    low = power[(freqs > 0.05) & (freqs < 0.1)].mean()
    high = power[(freqs > 10) & (freqs < 20)].mean()
    # Expected ratio of (fknee/f)^1 over the two bands
    expected = np.mean(1 / freqs[(freqs > 0.05) & (freqs < 0.1)]) / np.mean(
        1 / freqs[(freqs > 10) & (freqs < 20)]
    )
    assert np.isclose(low / high, expected, rtol=0.2)


@pytest.mark.parametrize("corr", [0.0, 0.6, 1.0])
def test_red_noise_pol_correlation(corr):
    """The correlation of the 1/f noise between channels equals pol_correlation."""
    noise = _red_noise(512, 0.04, 1.0, 1.0, nchan=2, correlation=corr, ntrials=400, rng=2)
    measured = np.mean([np.corrcoef(trial)[0, 1] for trial in noise])
    assert np.isclose(measured, corr, atol=0.05)


def test_offset_arm_shape():
    """An offset shifts the profile along the arm and attenuates it across the arm."""
    sigma = 2.5 / 60

    def func(x):
        return 100 * np.exp(-(x**2) / (2 * sigma**2))

    shifted = _offset_arm_shape(func, along=0.01, across=sigma, beam_sigma=sigma)
    x = np.linspace(-0.1, 0.1, 2001)
    y = shifted(x)
    assert np.isclose(x[np.argmax(y)], 0.01, atol=1e-4)
    assert np.isclose(y.max(), 100 * np.exp(-0.5), rtol=1e-4)
