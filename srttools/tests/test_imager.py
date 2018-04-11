# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function)
import numpy as np
import hashlib

try:
    import contextlib2 as contextlib
    FileNotFoundError = IOError
except ImportError:
    import contextlib

try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

try:
    import pyregion
    HAS_PYREGION = True
except ImportError:
    HAS_PYREGION = False

from srttools import ScanSet
from srttools import Scan
from srttools import CalibratorTable
from srttools.calibration import HAS_STATSM
from srttools.read_config import read_config
from srttools.imager import main_imager, main_preprocess, \
    _excluded_regions_from_args
from srttools.simulate import simulate_map
from srttools.global_fit import display_intermediate
from srttools.io import mkdir_p
from srttools.interactive_filter import intervals
import copy
import os
import glob
import astropy.units as u
import shutil
import pytest
import logging
import subprocess as sp
import astropy

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(x):
        return x


@pytest.fixture()
def logger():
    logger = logging.getLogger('Some.Logger')
    logger.setLevel(logging.INFO)

    return logger


np.random.seed(1241347)


def _md5(file):
    with open(file, 'rb') as fobj:
        string = fobj.read()

    return hashlib.md5(string).hexdigest()


def _2d_gauss(x, y, sigma=2.5 / 60.):
    """A Gaussian beam"""
    return np.exp(-(x ** 2 + y ** 2) / (2 * sigma**2))


def gauss_src_func(x, y):
    return 25 * _2d_gauss(x, y, sigma=2.5 / 60)


def sim_config_file(filename, add_garbage=False, prefix=None):
    """Create a sample config file, to be modified by hand."""
    string0 = """
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
    defective
"""
    string1 = """
calibrator_directories :
    calibration
"""
    string2 = """

skydip_directories :
    gauss_skydip

noise_threshold : 5

pixel_size : 0.8

[debugging]

debug_file_format : eps

"""
    if prefix is None:
        prefix = os.getcwd()
    import tempfile
    string = string0
    with open(filename, 'w') as fobj:
        print(string0, file=fobj)
        if add_garbage:
            for _ in range(100):
                garbage = '    ' + \
                    tempfile.NamedTemporaryFile(prefix=prefix).name[1:]
                print(garbage, file=fobj)
                string += garbage + '\n'
        print(string1, file=fobj)
        string += string1
        if add_garbage:
            for _ in range(100):
                garbage = '    ' + \
                    tempfile.NamedTemporaryFile(prefix=prefix).name[1:]
                print(garbage, file=fobj)
                string += garbage + '\n'
        print(string2, file=fobj)
        string += string2
    return string


def sim_map(obsdir_ra, obsdir_dec):
    simulate_map(count_map=gauss_src_func,
                 length_ra=30.,
                 length_dec=30.,
                 outdir=(obsdir_ra, obsdir_dec), mean_ra=180,
                 mean_dec=45, speed=1.5,
                 spacing=0.5, srcname='Dummy', channel_ratio=0.8,
                 baseline="flat")


class TestScanSet(object):
    @classmethod
    def setup_class(klass):
        import os

        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, 'data')
        klass.sim_dir = os.path.join(klass.datadir, 'sim')

        klass.obsdir_ra = os.path.join(klass.datadir, 'sim', 'gauss_ra')
        klass.obsdir_dec = os.path.join(klass.datadir, 'sim', 'gauss_dec')
        klass.config_file = \
            os.path.abspath(os.path.join(klass.sim_dir, 'test_config_sim.ini'))
        klass.caldir = os.path.join(klass.datadir, 'sim', 'calibration')
        klass.simulated_flux = 0.25
        # First off, simulate a beamed observation  -------

        sim_config_file(klass.config_file, add_garbage=True,
                        prefix="./")

        if (not os.path.exists(klass.obsdir_ra)) or \
                (not os.path.exists(klass.obsdir_dec)):
            mkdir_p(klass.obsdir_ra)
            mkdir_p(klass.obsdir_dec)
            print('Fake map: Point-like (but Gaussian beam shape), '
                  '{} Jy.'.format(klass.simulated_flux))
            sim_map(klass.obsdir_ra, klass.obsdir_dec)

        defective_dir = os.path.join(klass.sim_dir, 'defective')
        if not os.path.exists(defective_dir):
            shutil.copytree(os.path.join(klass.datadir, 'calibrators'),
                            defective_dir)

        # Copy skydip scan
        skydip_dir = os.path.join(klass.datadir, 'gauss_skydip')
        new_skydip_dir = os.path.join(klass.sim_dir, 'gauss_skydip')
        if os.path.exists(skydip_dir) and not os.path.exists(new_skydip_dir):
            shutil.copytree(skydip_dir, new_skydip_dir)
        caltable = CalibratorTable()
        caltable.from_scans(glob.glob(os.path.join(klass.caldir,
                                                   '*.fits')), debug=True)

        caltable.update()
        klass.calfile = os.path.join(klass.datadir, 'calibrators.hdf5')
        caltable.write(klass.calfile, overwrite=True)

        klass.config = read_config(klass.config_file)
        klass.raonly = os.path.abspath(os.path.join(klass.datadir,
                                                    'test_raonly.ini'))
        klass.deconly = os.path.abspath(os.path.join(klass.datadir,
                                                     'test_deconly.ini'))

        if HAS_PYREGION:
            excluded_xy, excluded_radec = \
                _excluded_regions_from_args([os.path.join(klass.datadir,
                                                          "center.reg")])
        else:
            excluded_xy, excluded_radec = None, None

        klass.scanset = ScanSet(klass.config_file, nosub=False,
                                norefilt=False,
                                debug=True, avoid_regions=excluded_radec)
        klass.scanset.write('test.hdf5', overwrite=True)

        klass.stdinfo = {}
        klass.stdinfo['FLAG'] = None
        klass.stdinfo['zap'] = intervals()
        klass.stdinfo['base'] = intervals()
        klass.stdinfo['fitpars'] = np.array([0, 0])

        def scan_no(scan_str):
            basename = os.path.splitext(os.path.basename(scan_str))[0]
            return int(basename.replace('Dec', '').replace('Ra', ''))

        klass.dec_scans = \
            dict([(scan_no(s), s)
                 for s in klass.scanset.scan_list if 'Dec' in s])
        klass.ra_scans = \
            dict([(scan_no(s), s)
                 for s in klass.scanset.scan_list if 'Ra' in s])
        klass.n_ra_scans = max(list(klass.ra_scans.keys()))
        klass.n_dec_scans = max(list(klass.dec_scans.keys()))

        if HAS_MPL:
            plt.ioff()

    def test_prepare(self):
        pass

    def test_script_is_installed(self):
        sp.check_call('SDTimage -h'.split(' '))

    def test_preprocess_single_files(self):
        files = glob.glob(os.path.join(self.obsdir_ra, '*.fits'))

        main_preprocess(files[:2] + ['--debug', '-c', self.config_file])
        for file in files[:2]:
            # I used debug_file_format : eps in the config
            if HAS_MPL:
                assert os.path.exists(file.replace('.fits', '_0.eps'))

    def test_script_is_installed_prep(self):
        sp.check_call('SDTpreprocess -h'.split(' '))

    def test_preprocess_no_config(self):
        with pytest.raises(ValueError) as excinfo:
            main_preprocess([])
        assert "Please specify the config file!" in str(excinfo)

    def test_preprocess_invalid(self):
        with pytest.warns(UserWarning) as record:
            main_preprocess([self.config_file])
        assert "is not in a known format" in record[0].message.args[0]

    def test_preprocess_config(self):
        main_preprocess(['-c', self.config_file])

    def test_imager_no_config(self):
        with pytest.raises(ValueError) as excinfo:
            main_imager([])
        assert "Please specify the config file!" in str(excinfo)

    def test_load_table_and_config(self):
        from astropy.table import Table
        table = Table.read('test.hdf5', path='scanset')
        scanset = ScanSet(table, config_file=self.config_file)
        for k in self.config.keys():
            assert scanset.meta[k] == self.config[k]

    def test_raonly(self):
        scanset = ScanSet(self.raonly)
        assert np.all(scanset['direction'])

    def test_deconly(self):
        scanset = ScanSet(self.deconly)
        assert not np.any(scanset['direction'])

    def test_multiple_tables(self):
        # scanset_all = ScanSet('test.hdf5')
        scanset = ScanSet([self.raonly, self.deconly])
        assert len(scanset.scan_list) == 122

    def test_wrong_file_name_raises(self):
        scanset = ScanSet('test.hdf5')
        with pytest.raises(astropy.io.registry.IORegistryError):
            scanset.write('asdlkfjsd.fjsdkf')

    @pytest.mark.skipif('not HAS_MPL')
    def test_interactive_quit(self):
        scanset = ScanSet('test.hdf5')
        imgsel = scanset.interactive_display('Feed0_RCP', test=True)
        fake_event = type('event', (), {})()
        fake_event.key = 'q'
        fake_event.xdata, fake_event.ydata = (130, 30)

        retval = imgsel.on_key(fake_event)
        assert retval == (130, 30, 'q')

    @pytest.mark.skipif('HAS_MPL')
    def test_interactive_quit_raises(self):
        scanset = ScanSet('test.hdf5')
        with pytest.raises(ImportError) as excinfo:
            imgsel = scanset.interactive_display('Feed0_RCP', test=True)
            assert "matplotlib is not installed" in str(excinfo)

    @pytest.mark.skipif('not HAS_MPL')
    def test_interactive_scans_all_calibrated_channels(self, capsys):
        scanset = ScanSet('test.hdf5')
        scanset.calibrate_images(calibration=self.calfile)
        images = scanset.images
        ysize, xsize = images['Feed0_RCP'].shape

        imgsel = scanset.interactive_display(test=True)
        fake_event = type('event', (), {})()
        fake_event.key = 'a'
        fake_event.xdata, fake_event.ydata = (xsize//2, ysize-1)

        imgsel.on_key(fake_event)
        fake_event.key = 'h'
        fake_event.xdata, fake_event.ydata = (xsize//2, ysize-1)
        out, err = capsys.readouterr()
        assert "a    open a window to filter all" in out

        imgsel.on_key(fake_event)
        fake_event.key = 'v'
        fake_event.xdata, fake_event.ydata = (xsize//2, ysize-1)
        imgsel.on_key(fake_event)

    def test_use_command_line(self):
        main_imager(('test.hdf5 -u Jy/beam ' +
                     '--calibrate {}'.format(self.calfile) +
                     ' -o bubu.hdf5 --debug --scrunch-channels').split(' '))

    def test_use_command_line_config(self):
        main_imager(['-c', self.config_file])

    def test_meta_saved_and_loaded_correctly(self):
        scanset = ScanSet('test.hdf5')
        for k in scanset.meta.keys():
            assert np.all(scanset.meta[k] == self.scanset.meta[k])
        assert sorted(scanset.meta.keys()) == sorted(self.scanset.meta.keys())
        assert scanset.scan_list == self.scanset.scan_list
        assert sorted(scanset.meta['calibrator_directories']) == \
            sorted(list(set(scanset.meta['calibrator_directories'])))
        assert sorted(scanset.meta['list_of_directories']) == \
            sorted(list(set(scanset.meta['list_of_directories'])))

    def test_barycenter_times(self):
        '''Test image production.'''

        scanset = ScanSet('test.hdf5')

        scanset.barycenter_times()

        assert 'barytime' in scanset.colnames
        assert np.all(np.abs(scanset['barytime'] -
                             scanset['time']) < 9 * 60 / 86400)

    def test_apply_bad_user_filt(self):
        scanset = ScanSet('test.hdf5')
        with pytest.raises(ValueError):
            scanset.apply_user_filter()

    def test_apply_user_filt(self):
        '''Test apply user filts.'''

        def user_fun(scanset):
            return scanset['barytime'] - np.floor(scanset['barytime'])

        scanset = ScanSet('test.hdf5')
        scanset.barycenter_times()
        phase_in_sec = scanset.apply_user_filter(user_fun, 'Phase_in_sec')

        assert np.min(phase_in_sec) >= 0
        assert np.max(phase_in_sec) <= 1
        assert np.all(phase_in_sec == scanset['Phase_in_sec'])

    def test_rough_image_nooffsets_nofilt(self):
        '''Test image production.'''

        scanset = ScanSet('test.hdf5')
        scanset.remove_column('Feed0_RCP-filt')
        images = scanset.calculate_images(no_offsets=True)

        img = images['Feed0_RCP']
        if HAS_MPL:
            fig = plt.figure('img')
            plt.imshow(img, origin='lower')
            plt.colorbar()
            plt.savefig('img_nooff_nofilt.png')
            plt.close(fig)

    def test_hor_images(self):
        '''Test image production.'''

        scanset = ScanSet('test.hdf5')
        scanset.remove_column('Feed0_RCP-filt')
        images = scanset.calculate_images(no_offsets=True, direction=0)

        img = images['Feed0_RCP']

    def test_ver_images(self):
        '''Test image production.'''

        scanset = ScanSet('test.hdf5')
        scanset.remove_column('Feed0_RCP-filt')
        images = scanset.calculate_images(no_offsets=True, direction=1)

        img = images['Feed0_RCP']

    def test_rough_image(self):
        '''Test image production.'''

        scanset = ScanSet('test.hdf5')

        images = scanset.calculate_images()

        img = images['Feed0_RCP']

        if HAS_MPL:
            fig = plt.figure('img')
            plt.imshow(img, origin='lower')
            plt.colorbar()
            plt.savefig('img.png')
            plt.close(fig)

    def test_rough_image_altaz(self):
        '''Test image production.'''
        scanset = ScanSet('test.hdf5')

        images = scanset.calculate_images(altaz=True)

        img = images['Feed0_RCP']

        if HAS_MPL:
            fig = plt.figure('img_altaz')
            plt.imshow(img, origin='lower')
            plt.colorbar()
            plt.savefig('img_altaz.png')
            scanset.save_ds9_images(save_sdev=True, altaz=True)
            plt.close(fig)

    def test_image_stdev(self):
        '''Test image production.'''

        scanset = ScanSet('test.hdf5')

        images = scanset.calculate_images()

        img = images['Feed0_RCP-Sdev']
        if HAS_MPL:
            fig = plt.figure('log(img-Sdev)')
            plt.imshow(np.log10(img), origin='lower')
            plt.colorbar()
            plt.ioff()
            plt.savefig('img_sdev.png')
            plt.close(fig)

    def test_image_scrunch(self):
        '''Test image production.'''

        scanset = ScanSet('test.hdf5')

        scanset.calculate_images()
        images = scanset.scrunch_images()

        img = images['TOTAL']

        if HAS_MPL:
            fig = plt.figure('img - scrunched')
            plt.imshow(img, origin='lower')
            plt.colorbar()
            img = images['TOTAL-Sdev']
            plt.savefig('img_scrunch.png')
            plt.close(fig)

            fig = plt.figure('img - scrunched - sdev')
            plt.imshow(img, origin='lower')
            plt.colorbar()
            plt.ioff()
            plt.savefig('img_scrunch_sdev.png')
            plt.close(fig)

    def test_calc_and_calibrate_image_pixel(self):
        scanset = ScanSet('test.hdf5')

        scanset.calibrate_images(calibration=self.calfile,
                                 map_unit="Jy/pixel")
        images = scanset.images
        img = images['Feed0_RCP']
        center = img.shape[0] // 2, img.shape[1] // 2
        shortest_side = np.min(img.shape)
        X, Y = np.meshgrid(np.arange(img.shape[1]), np.arange(img.shape[0]))
        good = (X-center[1])**2 + (Y-center[0])**2 <= (shortest_side//4)**2
        rtol = 0.2

        print(np.sum(images['Feed0_RCP'][good]), self.simulated_flux)
        assert np.isclose(np.sum(images['Feed0_RCP'][good]),
                          self.simulated_flux, rtol=rtol)

    def test_destripe(self):
        '''Test image production.'''

        scanset = ScanSet('test.hdf5')

        scanset.destripe_images(calibration=self.calfile,
                                map_unit="Jy/pixel", npix_tol=10)
        images = scanset.images
        img = images['Feed0_RCP']
        center = img.shape[0] // 2, img.shape[1] // 2
        shortest_side = np.min(img.shape)
        X, Y = np.meshgrid(np.arange(img.shape[1]), np.arange(img.shape[0]))
        good = (X-center[1])**2 + (Y-center[0])**2 <= (shortest_side//4)**2
        print(np.sum(images['Feed0_RCP'][good]), self.simulated_flux)
        assert np.isclose(np.sum(images['Feed0_RCP'][good]),
                          self.simulated_flux, rtol=0.1)

    def test_destripe_scrunch(self):
        '''Test image production.'''

        scanset = ScanSet('test.hdf5')

        scanset.destripe_images(calibration=self.calfile,
                                map_unit="Jy/pixel",
                                calibrate_scans=True)
        scanset.scrunch_images()
        images = scanset.images
        img = images['TOTAL']
        center = img.shape[0] // 2, img.shape[1] // 2
        shortest_side = np.min(img.shape)
        X, Y = np.meshgrid(np.arange(img.shape[1]), np.arange(img.shape[0]))
        good = (X-center[1])**2 + (Y-center[0])**2 <= (shortest_side//4)**2

        print(np.sum(images['TOTAL'][good]), self.simulated_flux)
        assert np.isclose(np.sum(images['TOTAL'][good]),
                          self.simulated_flux, rtol=0.1)

    def test_calibrate_image_pixel(self):
        scanset = ScanSet('test.hdf5')

        images = scanset.calculate_images(calibration=self.calfile,
                                          map_unit="Jy/pixel")

        img = images['Feed0_RCP']
        center = img.shape[0] // 2, img.shape[1] // 2
        shortest_side = np.min(img.shape)
        X, Y = np.meshgrid(np.arange(img.shape[1]), np.arange(img.shape[0]))
        good = (X-center[1])**2 + (Y-center[0])**2 <= (shortest_side//4)**2
        rtol = 0.1 if HAS_STATSM else 0.15

        print(np.sum(images['Feed0_RCP'][good]), self.simulated_flux)
        assert np.isclose(np.sum(images['Feed0_RCP'][good]),
                          self.simulated_flux, rtol=rtol)

    def test_calibrate_image_beam(self):
        scanset = ScanSet('test.hdf5')

        scanset.calculate_images()
        images = scanset.calculate_images(calibration=self.calfile,
                                          map_unit="Jy/beam")

        print(np.max(images['Feed0_RCP']), self.simulated_flux)
        assert np.allclose(np.max(images['Feed0_RCP']), self.simulated_flux,
                           atol=0.1)

    def test_calibrate_image_junk_unit_fails(self):
        scanset = ScanSet('test.hdf5')

        scanset.calculate_images()
        with pytest.raises(ValueError) as excinfo:
            images = scanset.calculate_images(calibration=self.calfile,
                                              map_unit="junk")
            assert "Unit for calibration not recognized" in str(excinfo)

    def test_calibrate_image_sr(self):
        scanset = ScanSet('test.hdf5')

        scanset.calculate_images()
        images = scanset.calculate_images(calibration=self.calfile,
                                          map_unit="Jy/sr")

        images_pix = scanset.calculate_images(calibration=self.calfile,
                                              map_unit="Jy/pixel")

        pixel_area = scanset.meta['pixel_size'] ** 2
        print(images['Feed0_RCP'],
              images_pix['Feed0_RCP'] / pixel_area.to(u.sr).value)
        assert np.allclose(images['Feed0_RCP'],
                           images_pix['Feed0_RCP'] / pixel_area.to(u.sr).value,
                           rtol=0.05)

    def test_calibrate_scanset_pixel(self):
        scanset = ScanSet('test.hdf5')
        images_standard = scanset.calculate_images(calibration=self.calfile,
                                                   map_unit="Jy/pixel")
        images = scanset.calculate_images(calibration=self.calfile,
                                          map_unit="Jy/pixel",
                                          calibrate_scans=True)
        print(images['Feed0_RCP'], images_standard['Feed0_RCP'])
        assert np.allclose(images['Feed0_RCP'], images_standard['Feed0_RCP'],
                           rtol=0.05)

    def test_calibrate_scanset_beam(self):
        scanset = ScanSet('test.hdf5')
        images_standard = scanset.calculate_images(calibration=self.calfile,
                                                   map_unit="Jy/beam")
        images = scanset.calculate_images(calibration=self.calfile,
                                          map_unit="Jy/beam",
                                          calibrate_scans=True)

        print(images['Feed0_RCP'], images_standard['Feed0_RCP'])
        assert np.allclose(images['Feed0_RCP'], images_standard['Feed0_RCP'],
                           atol=1e-3)

    def test_calibrate_scanset_sr(self):
        scanset = ScanSet('test.hdf5')
        images_standard = scanset.calculate_images(calibration=self.calfile,
                                                   map_unit="Jy/sr")
        images = scanset.calculate_images(calibration=self.calfile,
                                          map_unit="Jy/sr",
                                          calibrate_scans=True)

        good = images['Feed0_RCP'] > 1

        print(images['Feed0_RCP'][good], images_standard['Feed0_RCP'][good])
        assert np.allclose(images['Feed0_RCP'][good],
                           images_standard['Feed0_RCP'][good], rtol=0.05)

    def test_ds9_image(self):
        '''Test image production.'''

        scanset = ScanSet('test.hdf5')

        scanset.save_ds9_images(save_sdev=True)

    def test_ds9_image_not_save_sdev(self):
        '''Test image production.'''

        scanset = ScanSet('test.hdf5')

        scanset.save_ds9_images(save_sdev=False)

    def test_ds9_image_destripe(self):
        '''Test image production.'''

        scanset = ScanSet('test.hdf5')

        scanset.save_ds9_images(destripe=True)

    def test_global_fit_image(self):
        '''Test image production.'''
        scanset = ScanSet('test.hdf5')
        # It works with no parameters, before calculating images,
        # with no_offsets
        scanset.fit_full_images(no_offsets=True)
        # It works after calculating images
        images = scanset.calculate_images()
        nx, ny = images['Feed0_RCP'].shape
        excluded = [[nx // 2, ny // 2, nx // 4]]
        scanset.fit_full_images(excluded=excluded, chans='Feed0_RCP')
        os.path.exists("out_iter_Feed0_RCP_002.txt")
        scanset.fit_full_images(excluded=excluded, chans='Feed0_LCP')
        os.path.exists("out_iter_Feed0_LCP_000.txt")

        if not HAS_MPL:
            with pytest.raises(ImportError) as excinfo:
                display_intermediate(scanset, chan="Feed0_RCP",
                                     excluded=excluded,
                                     parfile="out_iter_Feed0_RCP_002.txt")

            assert "display_intermediate: matplotlib" in str(excinfo)
        else:
            display_intermediate(scanset, chan="Feed0_RCP",
                                 excluded=excluded,
                                 parfile="out_iter_Feed0_RCP_002.txt")
            os.path.exists("out_iter_Feed0_LCP_002.png")

        # It works after calculating images

    def test_find_scan_through_pixel0(self):
        scanset = ScanSet('test.hdf5')

        images = scanset.calculate_images()
        ysize, xsize = images['Feed0_RCP'].shape
        _, _, _, _, _, _, _, coord = \
            scanset.find_scans_through_pixel(xsize//2, ysize-1, test=True)

        dec_scan = os.path.join(self.obsdir_dec,
                                self.dec_scans[self.n_dec_scans // 2])
        assert dec_scan in coord
        assert coord[dec_scan] == 'dec'
        ra_scan = os.path.join(self.obsdir_ra,
                               self.ra_scans[self.n_ra_scans - 1])
        assert ra_scan in coord
        assert coord[ra_scan] == 'ra'

    def test_find_scan_through_pixel1(self):
        scanset = ScanSet('test.hdf5')
        for i in scanset.chan_columns:
            scanset.remove_column(i + '-filt')

        images = scanset.calculate_images()
        ysize, xsize = images['Feed0_RCP'].shape
        _, _, _, _, _, _, _, coord = \
            scanset.find_scans_through_pixel(xsize//2, 0, test=True)

        dec_scan = os.path.join(self.obsdir_dec,
                                self.dec_scans[self.n_dec_scans // 2])
        assert dec_scan in coord
        assert coord[dec_scan] == 'dec'
        ra_scan = os.path.join(self.obsdir_ra, 'Ra0.fits')
        assert ra_scan in coord
        assert coord[ra_scan] == 'ra'

    def test_find_scan_through_invalid_pixel(self):
        scanset = ScanSet('test.hdf5')

        images = scanset.calculate_images()
        ysize, xsize = images['Feed0_RCP'].shape
        _, _, _, _, _, _, _, coord = \
            scanset.find_scans_through_pixel(xsize//2, -2, test=True)
        assert coord == {}
        _, _, _, _, _, _, _, coord = \
            scanset.find_scans_through_pixel(xsize//2, ysize + 2, test=True)
        assert coord == {}

    def test_find_scan_through_pixel_bad_scan(self, logger, caplog):
        scanset = ScanSet('test.hdf5')
        images = scanset.calculate_images()
        ysize, xsize = images['Feed0_RCP'].shape
        x, y = xsize // 2, 0
        good_entries = np.logical_and(
                np.abs(scanset['x'][:, 0] - x) < 1,
                np.abs(scanset['y'][:, 0] - y) < 1)

        sids = list(set(scanset['Scan_id'][good_entries]))
        scanset.scan_list[sids[0]] = 'skd'
        scanset.find_scans_through_pixel(x, y, test=True)
        assert "Errors while opening scan skd" in caplog.text

    def test_update_scan_invalid(self):
        scanset = ScanSet('test.hdf5')

        images = scanset.calculate_images()
        ysize, xsize = images['Feed0_RCP'].shape
        ra_xs, ra_ys, dec_xs, dec_ys, scan_ids, ra_masks, dec_masks, coord = \
            scanset.find_scans_through_pixel(xsize//2, 0, test=True)

        sname = 'xkd'
        coord['xkd'] = None
        scan_ids['xkd'] = None

        info = {sname: copy.copy(self.stdinfo)}
        info[sname]['FLAG'] = True
        scanset.update_scan(sname, scan_ids[sname], coord[sname],
                            info[sname]['zap'],
                            info[sname]['fitpars'], info[sname]['FLAG'],
                            test=True)

    def test_update_scan_flag(self):
        scanset = ScanSet('test.hdf5')

        images = scanset.calculate_images()
        ysize, xsize = images['Feed0_RCP'].shape
        ra_xs, ra_ys, dec_xs, dec_ys, scan_ids, ra_masks, dec_masks, coord = \
            scanset.find_scans_through_pixel(xsize//3, 0, test=True)

        sname = list(scan_ids.keys())[0]

        info = {sname: copy.copy(self.stdinfo)}
        info[sname]['FLAG'] = True
        sid = scan_ids[sname]
        mask = scanset['Scan_id'] == sid
        before = scanset['Feed0_RCP-filt'][mask]
        scanset.update_scan(sname, scan_ids[sname], coord[sname],
                            info[sname]['zap'],
                            info[sname]['fitpars'], info[sname]['FLAG'],
                            test=True)
        after = scanset['Feed0_RCP-filt'][mask]
        assert np.all(before != after)
        s = Scan(sname)
        assert np.all(after == s['Feed0_RCP-filt'])
        assert s.meta['FLAG'] is True
        os.unlink(sname.replace('fits', 'hdf5'))

    def test_update_scan_fit(self):
        scanset = ScanSet('test.hdf5')

        images = scanset.calculate_images()
        ysize, xsize = images['Feed0_RCP'].shape
        ra_xs, ra_ys, dec_xs, dec_ys, scan_ids, ra_masks, dec_masks, coord = \
            scanset.find_scans_through_pixel(xsize//2, 0, test=True)

        sname = list(scan_ids.keys())[0]

        info = {sname: copy.copy(self.stdinfo)}
        info[sname]['fitpars'] = np.array([0.1, 0.3])
        sid = scan_ids[sname]
        mask = scanset['Scan_id'] == sid
        before = scanset['Feed0_RCP'][mask]
        scanset.update_scan(sname, scan_ids[sname], coord[sname],
                            info[sname]['zap'],
                            info[sname]['fitpars'], info[sname]['FLAG'],
                            test=True)
        after = scanset['Feed0_RCP'][mask]
        assert np.all(before != after)
        s = Scan(sname)
        assert np.all(after == s['Feed0_RCP'])
        os.unlink(sname.replace('fits', 'hdf5'))

    def test_update_scan_zap(self):
        scanset = ScanSet('test.hdf5')

        images = scanset.calculate_images()
        ysize, xsize = images['Feed0_RCP'].shape
        ra_xs, ra_ys, dec_xs, dec_ys, scan_ids, ra_masks, dec_masks, coord = \
            scanset.find_scans_through_pixel(xsize//2, 0, test=True)

        sname = sorted(list(dec_xs.keys()))[0]
        s = Scan(sname)

        info = {sname: copy.copy(self.stdinfo)}
        info[sname]['zap'].xs = [np.float(s['dec'][0]), np.float(s['dec'][10])]
        sid = scan_ids[sname]
        mask = scanset['Scan_id'] == sid
        before = scanset['Feed0_RCP-filt'][mask]
        scanset.update_scan(sname, scan_ids[sname], coord[sname],
                            info[sname]['zap'],
                            info[sname]['fitpars'], info[sname]['FLAG'],
                            test=True)
        after = scanset['Feed0_RCP-filt'][mask]

        assert np.all(before[:10] != after[:10])
        s = Scan(sname)

        assert np.all(np.array(after, dtype=bool) ==
                      np.array(s['Feed0_RCP-filt'], dtype=bool))
        os.unlink(sname.replace('fits', 'hdf5'))

    def test_imager_global_fit_invalid(self):
        '''Test image production.'''
        with pytest.raises(ValueError) as excinfo:
            main_imager('test.hdf5 -g -e 10 10 2 1'.split(' '))
            assert "Exclusion region has to be specified as " in str(excinfo)

    def test_imager_global_fit_valid(self):
        '''Test image production.'''
        # Get information on images
        scanset = ScanSet('test.hdf5')
        scanset.fit_full_images(no_offsets=True)
        # It works after calculating images
        images = scanset.calculate_images()
        nx, ny = images['Feed0_RCP'].shape
        excluded = [[nx//2, ny//2, nx//4]]

        main_imager('test.hdf5 -g '
                    '-e {} {} {}'.format(*(excluded[0])).split(' '))

    @pytest.mark.skipif('not HAS_PYREGION')
    def test_global_fit_image_using_ds9_region(self):
        scanset = ScanSet('test.hdf5')
        # It works after calculating images
        images = scanset.calculate_images()
        nx, ny = images['Feed0_RCP'].shape

        regstr = 'image;circle({},{},{})'.format(nx // 2, ny // 2, nx // 4)
        with open('region.reg', 'w') as fobj:
            print(regstr, file=fobj)

        main_imager('test.hdf5 -g --exclude region.reg'.split())
        os.unlink('region.reg')

    @pytest.mark.skipif('not HAS_PYREGION')
    def test_baseline_using_ds9_region(self):
        regstr = 'fk5;circle(180,45,960")'
        with open('region.reg', 'w') as fobj:
            print(regstr, file=fobj)

        main_imager(('-c {} --refilt '.format(self.config_file) +
                     '--sub --exclude region.reg').split())
        os.unlink('region.reg')

    @pytest.mark.skipif('not HAS_PYREGION')
    def test_baseline_preprocess_using_ds9_region(self):
        regstr = 'fk5;circle(180,45,960")'
        with open('region.reg', 'w') as fobj:
            print(regstr, file=fobj)

        main_preprocess(('-c {} '.format(self.config_file) +
                         '--sub --exclude region.reg').split())
        os.unlink('region.reg')

    @pytest.mark.skipif('not HAS_PYREGION')
    def test_global_fit_image_using_ds9_region_garbage_warns(self, logger,
                                                             caplog):
        regstr = 'physical;circle(30,30,1)'
        with open('region.reg', 'w') as fobj:
            print(regstr, file=fobj)
        main_imager('test.hdf5 -g --exclude region.reg'.split())
        assert 'Only regions in fk5' in caplog.text
        os.unlink('region.reg')

    @pytest.mark.skipif('not HAS_PYREGION')
    def test_global_fit_image_using_ds9_region_noncircular_warns(self, logger,
                                                                 caplog):
        regstr = 'image;line(100,100,200,200)'
        with open('region.reg', 'w') as fobj:
            print(regstr, file=fobj)
        main_imager('test.hdf5 -g --exclude region.reg'.split())
        assert 'Only circular regions' in caplog.text
        os.unlink('region.reg')

    def test_imager_sample_config(self):
        if os.path.exists('sample_config_file.ini'):
            os.unlink('sample_config_file.ini')
        with pytest.raises(SystemExit):
            main_imager(['--sample-config'])
        assert os.path.exists('sample_config_file.ini')

    @classmethod
    def teardown_class(klass):
        """Clean up the mess."""
        with contextlib.suppress(FileNotFoundError):
            if HAS_MPL:
                os.unlink('img.png')
                os.unlink('img_altaz.png')
                os.unlink('img_scrunch.png')
                os.unlink('delta_altaz.png')
                os.unlink('altaz_with_src.png')
                os.unlink('img_sdev.png')
                os.unlink('img_scrunch_sdev.png')
            os.unlink('test.hdf5')
            os.unlink('test_scan_list.txt')
            os.unlink('bubu.hdf5')
            for d in klass.config['list_of_directories']:
                hfiles = \
                    glob.glob(os.path.join(klass.config['datadir'], d,
                                           '*.hdf5'))
                for h in hfiles:
                    os.unlink(h)
            out_iter_files = glob.glob('out_iter_*')
            for o in out_iter_files:
                os.unlink(o)
            out_fits_files = glob.glob(os.path.join(klass.config['datadir'],
                                                    'test_config*.fits'))
            out_hdf5_files = glob.glob(os.path.join(klass.config['datadir'],
                                                    'sim',
                                                    '*/', '*.hdf5'))

            for o in out_fits_files + out_hdf5_files:
                os.unlink(o)
