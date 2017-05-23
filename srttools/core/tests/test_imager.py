# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function)
from ..read_config import read_config, SRT_tools_config
import numpy as np
import numpy.random as ra
import matplotlib.pyplot as plt
from astropy.table import Table
from ..imager import ScanSet
from ..simulate import simulate_map, save_scan
from ..global_fit import display_intermediate
from ..calibration import CalibratorTable
from ..io import mkdir_p
import os
import glob
import subprocess as sp

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(x):
        return x



np.random.seed(1241347)


def _2d_gauss(x, y, sigma=3 / 60.):
    """A Gaussian beam"""
    return np.exp(-(x ** 2 + y ** 2) / (2 * sigma**2))


def gauss_src_func(x, y):
    return 50 * _2d_gauss(x, y, sigma=3 / 60)


def calibrator_scan_func(x):
    return 100 * _2d_gauss(x, 0, sigma=3 / 60)


def sim_config_file(filename):
    """Create a sample config file, to be modified by hand."""
    string = """
[local]
workdir : .
datadir : .

[analysis]
projection : ARC
interpolation : spline
prefix : test_
list_of_directories :
    gauss_ra
    gauss_dec
calibrator_directories :
    calibration

noise_threshold : 5

pixel_size : 0.8
        """
    with open(filename, 'w') as fobj:
        print(string, file=fobj)
    return string


def sim_map(obsdir_ra, obsdir_dec):
    simulate_map(count_map=gauss_src_func,
                 length_ra=30.,
                 length_dec=30.,
                 outdir=(obsdir_ra, obsdir_dec), mean_ra=180,
                 mean_dec=70, speed=2.,
                 spacing=0.5, srcname='Dummy')


def sim_calibrators(ncross, caldir):
    src_ra = 185
    src_dec = 75
    timedelta = 0
    speed = 2.  # arcmin/s
    dt = 0.04
    dtheta = speed * dt
    scan_values = np.arange(-2, 2, dtheta/60)
    zero_values = np.zeros_like(scan_values)

    for i in tqdm(range(ncross)):
        ras = src_ra + scan_values / np.cos(np.radians(src_dec))
        decs = src_dec + zero_values
        times = np.arange(scan_values.size) * dt + timedelta

        scan = calibrator_scan_func(scan_values) + \
               ra.normal(0, 0.2, scan_values.size)
        save_scan(times, ras, decs, {'Ch0': scan, 'Ch1': scan},
                  filename=os.path.join(caldir, '{}_Ra.fits'.format(i)),
                  src_ra=src_ra, src_dec=src_dec, srcname='DummyCal')
        timedelta = times[-1] + 1

        ras = src_ra + zero_values
        decs = src_dec + scan_values
        times = np.arange(scan_values.size) * dt + timedelta

        scan = calibrator_scan_func(scan_values) + \
               ra.normal(0, 0.2, scan_values.size)
        save_scan(times, ras, decs, {'Ch0': scan, 'Ch1': scan},
                  filename=os.path.join(caldir, '{}_Dec.fits'.format(i)),
                  src_ra=src_ra, src_dec=src_dec, srcname='DummyCal')
        timedelta = times[-1] + 1


class TestScanSet(object):
    @classmethod
    def setup_class(klass):
        import os
        global DEBUG_MODE
        DEBUG_MODE = True

        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, 'data')
        klass.obsdir_ra = os.path.join(klass.datadir, 'sim', 'gauss_ra')
        klass.obsdir_dec = os.path.join(klass.datadir, 'sim', 'gauss_dec')
        klass.config_file = \
            os.path.abspath(os.path.join(klass.datadir, 'sim',
                                         'test_config_sim.ini'))
        klass.caldir = os.path.join(klass.datadir, 'sim', 'calibration')
        mkdir_p(klass.obsdir_ra)
        mkdir_p(klass.obsdir_dec)
        mkdir_p(klass.caldir)
        # First off, simulate a beamed observation  -------

        print('Setting up simulated data.')
        sim_config_file(klass.config_file)
        print('Fake calibrators: DummyCal, 1 Jy.')
        sim_calibrators(5, klass.caldir)
        print('Fake map: Point-like (but Gaussian beam shape), 0.5 Jy.')
        sim_map(klass.obsdir_ra, klass.obsdir_dec)

        caltable = CalibratorTable()
        caltable.from_scans(glob.glob(os.path.join(klass.caldir,
                                                   '*.fits')))

        caltable.update()
        klass.calfile = os.path.join(klass.datadir, 'calibrators.hdf5')
        caltable.write(klass.calfile,
                       path="config", overwrite=True)


        klass.config = read_config(klass.config_file)
        klass.scanset = ScanSet(klass.config_file, norefilt=False)

        klass.scanset.write('test.hdf5', overwrite=True)

        plt.ioff()

    def test_0_prepare(self):
        pass

    def test_1_meta_saved_and_loaded_correctly(self):
        scanset = ScanSet('test.hdf5',
                          config_file=self.config_file)
        for k in scanset.meta.keys():
            assert np.all(scanset.meta[k] == self.scanset.meta[k])
        assert sorted(scanset.meta.keys()) == sorted(self.scanset.meta.keys())
        assert scanset.scan_list == self.scanset.scan_list

    def test_2_rough_image(self):
        '''Test image production.'''

        # scanset = ScanSet(Table.read('test.hdf5', path='scanset'),
        #                   config_file=self.config_file)

        scanset = ScanSet('test.hdf5',
                          config_file=self.config_file)

        images = scanset.calculate_images()

        img = images['Ch0']

        fig = plt.figure('img')
        plt.imshow(img, origin='lower')
        plt.colorbar()
        plt.savefig('img.png')
        plt.close(fig)

    def test_3_rough_image_altaz(self):
        '''Test image production.'''

        scanset = ScanSet(Table.read('test.hdf5', path='scanset'),
                          config_file=self.config_file)

        images = scanset.calculate_images(altaz=True)

        img = images['Ch0']

        fig = plt.figure('img_altaz')
        plt.imshow(img, origin='lower')
        plt.colorbar()
        plt.savefig('img_altaz.png')
        scanset.save_ds9_images(save_sdev=True, altaz=True)
        plt.close(fig)

    def test_4_image_stdev(self):
        '''Test image production.'''

        scanset = ScanSet(Table.read('test.hdf5', path='scanset'),
                          config_file=self.config_file)

        images = scanset.calculate_images()

        img = images['Ch0-Sdev']

        fig = plt.figure('log(img-Sdev)')
        plt.imshow(np.log10(img), origin='lower')
        plt.colorbar()
        plt.ioff()
        plt.savefig('img_sdev.png')
        plt.close(fig)

    def test_5_image_scrunch(self):
        '''Test image production.'''

        scanset = ScanSet(Table.read('test.hdf5', path='scanset'),
                          config_file=self.config_file)

        images = scanset.calculate_images(scrunch=True)

        img = images['Ch0']

        fig = plt.figure('img - scrunched')
        plt.imshow(img, origin='lower')
        plt.colorbar()
        img = images['Ch0-Sdev']
        plt.savefig('img_scrunch.png')
        plt.close(fig)

        fig = plt.figure('img - scrunched - sdev')
        plt.imshow(img, origin='lower')
        plt.colorbar()
        plt.ioff()
        plt.savefig('img_scrunch_sdev.png')
        plt.close(fig)

    def test_6_calibrate_image(self):
        scanset = ScanSet(Table.read('test.hdf5', path='scanset'),
                          config_file=self.config_file)
        images = scanset.calculate_images()

        images = scanset.calculate_images(calibration=self.calfile)

        assert np.allclose(np.sum(images['Ch0']), 0.5)

    def test_7_ds9_image(self):
        '''Test image production.'''

        scanset = ScanSet('test.hdf5',
                          config_file=self.config_file)

        scanset.save_ds9_images(save_sdev=True)

    def test_8_global_fit_image(self):
        '''Test image production.'''

        scanset = ScanSet(Table.read('test.hdf5', path='scanset'),
                          config_file=self.config_file)
        images = scanset.calculate_images()
        nx, ny = images['Ch0'].shape
        excluded = [[nx//2, ny//2, nx//4]]
        scanset.fit_full_images(excluded=excluded, chans='Ch0')
        os.path.exists("out_iter_Ch0_002.txt")
        scanset.fit_full_images(excluded=excluded, chans='Ch1')
        os.path.exists("out_iter_Ch1_000.txt")
        display_intermediate(scanset, chan="Ch0", excluded=excluded,
                             parfile="out_iter_Ch0_002.txt")
        os.path.exists("out_iter_Ch1_002.png")

    def test_9_find_scan_through_pixel(self):
        scanset = ScanSet('test.hdf5',
                          config_file=self.config_file)

        scanset.calculate_images()
        scanset.find_scans_through_pixel(125, 125, test=True)

    @classmethod
    def teardown_class(klass):
        """Clean up the mess."""
        os.unlink('img.png')
        os.unlink('img_altaz.png')
        os.unlink('img_scrunch.png')
        os.unlink('delta_altaz.png')
        os.unlink('altaz_with_src.png')
        os.unlink('img_sdev.png')
        # os.unlink('test.hdf5')
        for d in klass.config['list_of_directories']:
            hfiles = \
                glob.glob(os.path.join(klass.config['datadir'], d, '*.hdf5'))
            for h in hfiles:
                os.unlink(h)
        out_iter_files = glob.glob('out_iter_*.txt')
        for o in out_iter_files:
            os.unlink(o)
        out_fits_files = glob.glob(os.path.join(klass.config['datadir'],
                                                'test_config*.fits'))
        out_hdf5_files = glob.glob(os.path.join(klass.config['datadir'], 'sim',
                                                '*/', '*.hdf5'))

        for o in out_fits_files + out_hdf5_files:
            os.unlink(o)
