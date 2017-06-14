# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function)
from ..read_config import read_config
import numpy as np
import matplotlib.pyplot as plt
from ..imager import ScanSet, main_imager
from ..simulate import simulate_map
from ..global_fit import display_intermediate
from ..calibration import CalibratorTable
from ..io import mkdir_p
import os
import glob
import astropy.units as u

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
                 length_ra=50.,
                 length_dec=50.,
                 outdir=(obsdir_ra, obsdir_dec), mean_ra=180,
                 mean_dec=70, speed=2.,
                 spacing=0.5, srcname='Dummy')


class TestScanSet(object):
    @classmethod
    def setup_class(klass):
        import os

        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, 'data')
        klass.obsdir_ra = os.path.join(klass.datadir, 'sim', 'gauss_ra')
        klass.obsdir_dec = os.path.join(klass.datadir, 'sim', 'gauss_dec')
        klass.config_file = \
            os.path.abspath(os.path.join(klass.datadir, 'sim',
                                         'test_config_sim.ini'))
        klass.caldir = os.path.join(klass.datadir, 'sim', 'calibration')
        # First off, simulate a beamed observation  -------

        if not os.path.exists(klass.config_file):
            print('Setting up simulated data.')
            sim_config_file(klass.config_file)

        if (not os.path.exists(klass.obsdir_ra)) or \
                (not os.path.exists(klass.obsdir_dec)):
            mkdir_p(klass.obsdir_ra)
            mkdir_p(klass.obsdir_dec)
            print('Fake map: Point-like (but Gaussian beam shape), 0.5 Jy.')
            sim_map(klass.obsdir_ra, klass.obsdir_dec)

        caltable = CalibratorTable()
        caltable.from_scans(glob.glob(os.path.join(klass.caldir,
                                                   '*.fits')), debug=True)

        caltable.update()
        klass.calfile = os.path.join(klass.datadir, 'calibrators.hdf5')
        caltable.write(klass.calfile, overwrite=True)

        klass.config = read_config(klass.config_file)

        klass.scanset = ScanSet(klass.config_file, norefilt=False,
                                debug=True)
        klass.scanset.write('test.hdf5', overwrite=True)

        plt.ioff()

    def test_0_prepare(self):
        pass

    def test_interactive_quit(self):
        scanset = ScanSet('test.hdf5')
        imgsel = scanset.interactive_display('Ch0', test=True)
        fake_event = type('event', (), {})()
        fake_event.key = 'q'
        fake_event.xdata, fake_event.ydata = (130, 30)

        retval = imgsel.on_key(fake_event)
        assert retval == (130, 30, 'q')

    def test_interactive_scans(self):
        scanset = ScanSet('test.hdf5')
        imgsel = scanset.interactive_display('Ch0', test=True)
        fake_event = type('event', (), {})()
        fake_event.key = 'a'
        fake_event.xdata, fake_event.ydata = (62, 62)

        imgsel.on_key(fake_event)

    def test_use_command_line(self):
        main_imager(('test.hdf5 -u Jy/beam ' +
                     '--calibrate {}'.format(self.calfile) +
                     ' -o bubu.hdf5 --debug').split(' '))

    def test_1_meta_saved_and_loaded_correctly(self):
        scanset = ScanSet('test.hdf5')
        for k in scanset.meta.keys():
            assert np.all(scanset.meta[k] == self.scanset.meta[k])
        assert sorted(scanset.meta.keys()) == sorted(self.scanset.meta.keys())
        assert scanset.scan_list == self.scanset.scan_list
        assert scanset.meta['calibrator_directories'] == \
            list(set(scanset.meta['calibrator_directories']))
        assert scanset.meta['list_of_directories'] == \
            list(set(scanset.meta['list_of_directories']))

    def test_2_rough_image(self):
        '''Test image production.'''

        scanset = ScanSet('test.hdf5')

        images = scanset.calculate_images()

        img = images['Ch0']

        fig = plt.figure('img')
        plt.imshow(img, origin='lower')
        plt.colorbar()
        plt.savefig('img.png')
        plt.close(fig)

    def test_3_rough_image_altaz(self):
        '''Test image production.'''

        scanset = ScanSet('test.hdf5')

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

        scanset = ScanSet('test.hdf5')

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

        scanset = ScanSet('test.hdf5')

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

    def test_6a_calibrate_image_pixel(self):
        scanset = ScanSet('test.hdf5')

        scanset.calculate_images()
        images = scanset.calculate_images(calibration=self.calfile,
                                          map_unit="Jy/pixel")

        img = images['Ch0']
        center = img.shape[0] // 2, img.shape[1] // 2
        shortest_side = np.min(img.shape)
        X, Y = np.meshgrid(np.arange(img.shape[1]), np.arange(img.shape[0]))
        good = (X-center[1])**2 + (Y-center[0])**2 <= (shortest_side//4)**2
        assert np.allclose(np.sum(images['Ch0'][good]), 0.5, 0.05)

    def test_6b_calibrate_image_beam(self):
        scanset = ScanSet('test.hdf5')

        scanset.calculate_images()
        images = scanset.calculate_images(calibration=self.calfile,
                                          map_unit="Jy/beam")

        assert np.allclose(np.max(images['Ch0']), 0.5, atol=0.05)

    def test_6c_calibrate_image_sr(self):
        scanset = ScanSet('test.hdf5')

        scanset.calculate_images()
        images = scanset.calculate_images(calibration=self.calfile,
                                          map_unit="Jy/sr")

        images_pix = scanset.calculate_images(calibration=self.calfile,
                                              map_unit="Jy/pixel")

        pixel_area = scanset.meta['pixel_size'] ** 2
        assert np.allclose(images['Ch0'],
                           images_pix['Ch0'] / pixel_area.to(u.sr).value,
                           rtol=0.05)

    def test_7_ds9_image(self):
        '''Test image production.'''

        scanset = ScanSet('test.hdf5')

        scanset.save_ds9_images(save_sdev=True)

    def test_8_global_fit_image(self):
        '''Test image production.'''

        scanset = ScanSet('test.hdf5')
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

    def test_9a_find_scan_through_pixel(self):
        scanset = ScanSet('test.hdf5')

        scanset.calculate_images()
        _, _, _, _, _, _, _, coord = \
            scanset.find_scans_through_pixel(62, 62, test=True)

        dec_scan = os.path.join(self.obsdir_dec, 'Dec99.fits')
        assert dec_scan in coord
        assert coord[dec_scan] == 'dec'
        ra_scan = os.path.join(self.obsdir_ra, 'Ra100.fits')
        assert ra_scan in coord
        assert coord[ra_scan] == 'ra'

    def test_9b_find_scan_through_pixel(self):
        scanset = ScanSet('test.hdf5')

        scanset.calculate_images()
        _, _, _, _, _, _, _, coord = \
            scanset.find_scans_through_pixel(62, 0, test=True)

        dec_scan = os.path.join(self.obsdir_dec, 'Dec99.fits')
        assert dec_scan in coord
        assert coord[dec_scan] == 'dec'
        ra_scan = os.path.join(self.obsdir_ra, 'Ra0.fits')
        assert ra_scan in coord
        assert coord[ra_scan] == 'ra'

    def test_9c_find_scan_through_invalid_pixel(self):
        scanset = ScanSet('test.hdf5')

        scanset.calculate_images()
        _, _, _, _, _, _, _, coord = \
            scanset.find_scans_through_pixel(62, -2, test=True)
        assert coord == {}
        _, _, _, _, _, _, _, coord = \
            scanset.find_scans_through_pixel(62, 64, test=True)
        assert coord == {}

    @classmethod
    def teardown_class(klass):
        """Clean up the mess."""
        os.unlink('img.png')
        os.unlink('img_altaz.png')
        os.unlink('img_scrunch.png')
        os.unlink('delta_altaz.png')
        os.unlink('altaz_with_src.png')
        os.unlink('img_sdev.png')
        os.unlink('test.hdf5')
        os.unlink('test_scan_list.txt')
        os.unlink('img_scrunch_sdev.png')
        os.unlink('bubu.hdf5')
        for d in klass.config['list_of_directories']:
            hfiles = \
                glob.glob(os.path.join(klass.config['datadir'], d, '*.hdf5'))
            for h in hfiles:
                os.unlink(h)
        out_iter_files = glob.glob('out_iter_*')
        for o in out_iter_files:
            os.unlink(o)
        out_fits_files = glob.glob(os.path.join(klass.config['datadir'],
                                                'test_config*.fits'))
        out_hdf5_files = glob.glob(os.path.join(klass.config['datadir'], 'sim',
                                                '*/', '*.hdf5'))

        for o in out_fits_files + out_hdf5_files:
            os.unlink(o)
