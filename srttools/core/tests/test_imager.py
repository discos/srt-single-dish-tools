# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function)
from ..read_config import read_config, SRT_tools_config
import numpy as np
import matplotlib.pyplot as plt
import unittest
from astropy.table import Table
from ..imager import ScanSet
import os
import glob


class Test2_ScanSet(unittest.TestCase):
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

    def step_1_scanset(self):
        '''Test that sets of data are read.'''
        plt.ioff()

        scanset = ScanSet(self.config_file, norefilt=False)

        scanset.write('test.hdf5', overwrite=True)

    def step_2_rough_image(self):
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

    def step_3_rough_image_altaz(self):
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

    def step_4_image_stdev(self):
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

    def step_5_image_scrunch(self):
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
    #
    # def step_6_interactive_image(self):
    #     '''Test image production.'''
    #
    #     scanset = ScanSet(Table.read('test.hdf5', path='scanset'),
    #                       config_file=self.config)
    #
    #     scanset.calculate_images()
    #     scanset.interactive_display()

    def step_7_ds9_image(self):
        '''Test image production.'''

        scanset = ScanSet(Table.read('test.hdf5', path='scanset'),
                          config_file=self.config_file)

        scanset.save_ds9_images(save_sdev=True)

    def step_8_global_fit_image(self):
        '''Test image production.'''

        scanset = ScanSet(Table.read('test.hdf5', path='scanset'),
                          config_file=self.config_file)
        excluded = [[125, 125, 30]]
        scanset.fit_full_images(excluded=excluded, chans='Ch0')


    def step_999_cleanup(self):
        """Clean up the mess."""
        os.unlink('img.png')
        os.unlink('img_altax.png')
        os.unlink('img_scrunch.png')
        os.unlink('img_scrunch_sdev.png')
        os.unlink('img_sdev.png')
        os.unlink('test.hdf5')
        for d in self.config['list_of_directories']:
            hfiles = glob.glob(os.path.join(self.config['datadir'], d, '*.hdf5'))
            print(hfiles)
            for h in hfiles:
                os.unlink(h)

    def test_all(self):
        self.step_1_scanset()

        self.step_2_rough_image()
        self.step_3_rough_image_altaz()
        self.step_4_image_stdev()
        self.step_5_image_scrunch()
        # self.step_6_interactive_image()
        self.step_7_ds9_image()
        self.step_8_global_fit_image()
        self.step_999_cleanup()

#
# class Test3_MultiFeed(unittest.TestCase):
#     @classmethod
#     def setup_class(klass):
#         import os
#         global DEBUG_MODE
#         print('Setting up class')
#         DEBUG_MODE = True
#
#         klass.curdir = os.path.dirname(__file__)
#         klass.datadir = os.path.join(klass.curdir, '..', '..', 'TEST_DATASET')
#
#         klass.config = \
#             os.path.abspath(os.path.join(klass.curdir, '..', '..',
#                                          'TEST_DATASET',
#                                          'test_config.ini'))
#
#         read_config(klass.config)
#
#     def step1_scanset(self):
#         '''Test that sets of data are read also with multifeed.'''
#         plt.ioff()
#         scanset = ScanSet(self.config, norefilt=True)
#
#         scanset.write('test_multifeed.hdf5', overwrite=True)
#
#     def step2_img(self):
#         scanset = ScanSet(Table.read('test_multifeed.hdf5', path='scanset'),
#                           config_file=self.config)
#
#         images = scanset.calculate_images()
#
#         scanset.save_ds9_images('multifeed.fits')
#
#         img = images['Ch0']
#
#         plt.figure('img 0')
#         plt.imshow(img, origin='lower')
#         plt.colorbar()
#
#         img = images['Ch0-Sdev']
#
#         plt.figure('log(img 0 -Sdev)')
#         plt.imshow(np.log10(img), origin='lower')
#         plt.colorbar()
#
#     def step3_scrunch(self):
#         scanset = ScanSet(Table.read('test_multifeed.hdf5', path='scanset'),
#                           config_file=self.config)
#         images = scanset.calculate_images(scrunch=True)
#         scanset.save_ds9_images('scrunch.fits', scrunch=True)
#
#         img = images['Ch0']
#
#         plt.figure('img - scrunched')
#         plt.imshow(img, origin='lower')
#         plt.colorbar()
#         img = images['Ch0-Sdev']
#
#         plt.figure('log(img - scrunched - sdev)')
#         plt.imshow(np.log10(img), origin='lower')
#         plt.colorbar()
#         plt.ioff()
#         plt.show()
#
#     def step4_nocorrect(self):
#         scanset = ScanSet(Table.read('test_multifeed.hdf5', path='scanset'),
#                           config_file=self.config)
#         scanset.save_ds9_images('multifeed_nooffsets.fits', no_offsets=True)
#
#     def step5_altaz(self):
#         scanset = ScanSet(Table.read('test_multifeed.hdf5', path='scanset'),
#                           config_file=self.config)
#         scanset.save_ds9_images('multifeed_nooffsets_altaz.fits',
#                                 no_offsets=True,
#                                 altaz=True)
#
#     def step6_interactive_image(self):
#         '''Test image production.'''
#
#         scanset = ScanSet(Table.read('test_multifeed.hdf5', path='scanset'),
#                           config_file=self.config)
#
#         scanset.interactive_display()
#
#     def test_all_multifeed(self):
#         '''Test that sets of data are read also with multifeed.'''
#         self.step1_scanset()
#         self.step2_img()
#         self.step3_scrunch()
#         self.step4_nocorrect()
#         self.step6_interactive_image()
