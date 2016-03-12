# -*- coding: utf-8 -*-

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
from ..io import DEBUG_MODE
from ..read_config import read_config
import numpy as np
import matplotlib.pyplot as plt
import unittest
from astropy.table import Table
from ..scan import Scan, ScanSet


# class Test1_Scan(unittest.TestCase):
#     @classmethod
#     def setup_class(klass):
#         import os
#         global DEBUG_MODE
#         DEBUG_MODE = True
#
#         klass.curdir = os.path.dirname(__file__)
#         klass.datadir = os.path.join(klass.curdir, '..', '..', 'TEST_DATASET')
#
#         klass.fname = \
#             os.path.abspath(
#                 os.path.join(klass.datadir, '20140603-103246-scicom-3C157',
#                              '20140603-103246-scicom-3C157_003_003.fits'))
#
#         klass.config_file = \
#             os.path.abspath(os.path.join(klass.curdir, '..', '..',
#                                          'TEST_DATASET',
#                                          'test_config.ini'))
#
#         read_config(klass.config_file)
#
#     def step_1_scan(self):
#         '''Test that data are read.'''
#
#         scan = Scan(self.fname)
#
#         scan.write('scan.hdf5', overwrite=True)
#
#     def step_2_read_scan(self):
#         scan = Scan('scan.hdf5')
#         plt.ion()
#         for col in scan.chan_columns():
#             plt.plot(scan['time'], scan[col])
#         plt.ioff()
#         plt.show()
#
#         return scan
#
#     def test_all(self):
#         self.step_1_scan()
#         self.step_2_read_scan()
#

class Test2_ScanSet(unittest.TestCase):
    @classmethod
    def setup_class(klass):
        import os
        global DEBUG_MODE
        DEBUG_MODE = True

        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, '..', '..', 'TEST_DATASET')

        klass.fname = \
            os.path.abspath(
                os.path.join(klass.datadir, '20140603-103246-scicom-3C157',
                             '20140603-103246-scicom-3C157_003_003.fits'))

        klass.config = \
            os.path.abspath(os.path.join(klass.curdir, '..', '..',
                                         'TEST_DATASET',
                                         'test_config.ini'))

        read_config(klass.config)

    def step_1_scanset(self):
        '''Test that sets of data are read.'''
        plt.ioff()

        scanset = ScanSet(self.config, norefilt=False)
        print(scanset)

        scanset.write('test.hdf5', overwrite=True)

    def step_2_rough_image(self):
        '''Test image production.'''

        scanset = ScanSet(Table.read('test.hdf5', path='scanset'),
                          config_file=self.config)

        images = scanset.calculate_images()

        img = images['Ch0']

        plt.figure('img')
        plt.imshow(img, origin='lower')
        plt.colorbar()
        plt.show()

    def step_3_rough_image_altaz(self):
        '''Test image production.'''

        scanset = ScanSet(Table.read('test.hdf5', path='scanset'),
                          config_file=self.config)

        images = scanset.calculate_images(altaz=True)

        img = images['Ch0']

        plt.figure('img_altaz')
        plt.imshow(img, origin='lower')
        plt.colorbar()
        plt.show()

    def step_4_image_stdev(self):
        '''Test image production.'''

        scanset = ScanSet(Table.read('test.hdf5', path='scanset'),
                          config_file=self.config)

        images = scanset.calculate_images()

        img = images['Ch0-Sdev']

        plt.figure('log(img-Sdev)')
        plt.imshow(np.log10(img), origin='lower')
        plt.colorbar()
        plt.ioff()
        plt.show()

    def step_5_image_scrunch(self):
        '''Test image production.'''

        scanset = ScanSet(Table.read('test.hdf5', path='scanset'),
                          config_file=self.config)

        images = scanset.calculate_images(scrunch=True)

        img = images['Ch0']

        plt.figure('img - scrunched')
        plt.imshow(img, origin='lower')
        plt.colorbar()
        img = images['Ch0-Sdev']

        plt.figure('img - scrunched - sdev')
        plt.imshow(img, origin='lower')
        plt.colorbar()
        plt.ioff()
        plt.show()

    def step_6_interactive_image(self):
        '''Test image production.'''

        scanset = ScanSet(Table.read('test.hdf5', path='scanset'),
                          config_file=self.config)

        scanset.calculate_images()
        scanset.interactive_display()

    def step_7_ds9_image(self):
        '''Test image production.'''

        scanset = ScanSet(Table.read('test.hdf5', path='scanset'),
                          config_file=self.config)

        scanset.save_ds9_images(save_sdev=True)

    def test_all(self):
        self.step_1_scanset()
        plt.ion()
        self.step_2_rough_image()
        # self.step_3_rough_image_altaz()
        # self.step_4_image_stdev()
        self.step_5_image_scrunch()
        plt.ioff()
        self.step_6_interactive_image()
        self.step_7_ds9_image()

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
