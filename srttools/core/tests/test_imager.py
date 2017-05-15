# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function)
from ..read_config import read_config, SRT_tools_config
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from ..imager import ScanSet
import os
import glob


class TestScanSet(object):
    @classmethod
    def setup_class(klass):
        import os
        global DEBUG_MODE
        DEBUG_MODE = True

        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, 'data')

        klass.config_file = \
            os.path.abspath(os.path.join(klass.datadir, 'test_config.ini'))

        klass.config = read_config(klass.config_file)

    def test_1_scanset(self):
        '''Test that sets of data are read.'''
        plt.ioff()

        scanset = ScanSet(self.config_file, norefilt=False)

        scanset.write('test.hdf5', overwrite=True)

    def test_2_rough_image(self):
        '''Test image production.'''

        scanset = ScanSet(Table.read('test.hdf5', path='scanset'),
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
        plt.savefig('img_altax.png')
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

    def test_7_ds9_image(self):
        '''Test image production.'''

        scanset = ScanSet(Table.read('test.hdf5', path='scanset'),
                          config_file=self.config_file)

        scanset.save_ds9_images(save_sdev=True)

    def test_8_global_fit_image(self):
        '''Test image production.'''

        scanset = ScanSet(Table.read('test.hdf5', path='scanset'),
                          config_file=self.config_file)
        excluded = [[125, 125, 30]]
        scanset.fit_full_images(excluded=excluded, chans='Ch0')

    @classmethod
    def teardown_class(klass):
        """Clean up the mess."""
        os.unlink('img.png')
        os.unlink('img_altax.png')
        os.unlink('img_scrunch.png')
        os.unlink('img_scrunch_sdev.png')
        os.unlink('img_sdev.png')
        os.unlink('test.hdf5')
        for d in klass.config['list_of_directories']:
            hfiles = glob.glob(os.path.join(klass.config['datadir'], d, '*.hdf5'))
            print(hfiles)
            for h in hfiles:
                os.unlink(h)
        out_iter_files = glob.glob('out_iter_*.txt')
        for o in out_iter_files:
            os.unlink(o)
