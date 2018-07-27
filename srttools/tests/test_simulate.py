import numpy as np
import shutil
import os
import pytest
import subprocess as sp
from astropy.io import fits
from srttools.io import mkdir_p
from srttools.simulate import simulate_map, simulate_scan, \
    sim_position_switching
from srttools.simulate import _default_map_shape, save_scan, main_simulate


class TestSimulate(object):
    @classmethod
    def setup_class(cls):
        cls.outdir = os.path.join('sim')
        cls.emptydir = os.path.join('sim', 'empty')
        cls.pswdir = os.path.join('sim', 'psw')
        for d in [cls.emptydir, cls.pswdir]:
            mkdir_p(d)

    def test_script_is_installed(self):
        sp.check_call('SDTfake -h'.split(' '))

    def test_sim_scan(self):
        """Test the simulation of a single scan."""
        times, position, shape = simulate_scan()
        save_scan(times, position, np.zeros_like(position),
                  {'Ch0': shape, 'Ch1': shape},
                  os.path.join(self.outdir, 'output.fits'))

    def test_sim_scan_other_columns(self):
        """Test the simulation of a single scan."""
        times, position, shape = simulate_scan()
        save_scan(times, position, np.zeros_like(position),
                  {'Ch0': shape, 'Ch1': shape},
                  os.path.join(self.outdir, 'output.fits'),
                  scan_type="Any",
                  other_columns={'Pippo': np.zeros_like(shape)})

    def test_position_switching(self):
        """Test the simulation of an empty map."""
        sim_position_switching(self.pswdir)

    def test_sim_map_empty(self):
        """Test the simulation of an empty map."""
        out_ra, _ = \
            simulate_map(width_ra=2, width_dec=2., outdir=self.emptydir)
        probe = os.path.join(out_ra, 'Ra0.fits')
        assert os.path.exists(probe)
        with fits.open(probe) as hdul:
            assert hdul[0].header['Declination Offset'] != 0.

    def test_sim_map_empty_messy(self):
        """Test the simulation of an empty map."""
        simulate_map(width_ra=2, width_dec=2., outdir=self.emptydir,
                     baseline='messy', nbin=100)

    def test_sim_map_empty_slope(self):
        """Test the simulation of an empty map."""
        simulate_map(width_ra=2, width_dec=2., outdir=self.emptydir,
                     baseline='slope')

    def test_raises_wrong_map_shape(self):
        with pytest.raises(ValueError):
            _default_map_shape(np.zeros((3, 4)), np.ones((3, 6)))

    def test_use_script(self):
        main_simulate('--no-cal -s 0.005 -g 10 10 1 1 -o sima -b 1'.split(' '))
        shutil.rmtree('sima')
        main_simulate('--no-cal -g 10 10 1 1 -o simb -b messy'.split(' '))
        shutil.rmtree('simb')
        main_simulate('--no-cal -g 10 10 1 1 -o simc -b slope'.split(' '))
        shutil.rmtree('simc')
        main_simulate('--integration-time .02 -g 10 10 1 1 -o simd'.split(' '))
        shutil.rmtree('simd')
        main_simulate('--no-cal --scan-speed 3. -g 10 10 1 1 '
                      '-o sime'.split(' '))
        shutil.rmtree('sime')

    def test_use_wrong_baseline(self):
        with pytest.raises(ValueError):
            main_simulate(
                '-b qwerty -g 10 10 1 1 -o sime'.split(' '))

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.outdir)
