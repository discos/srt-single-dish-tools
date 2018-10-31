from __future__ import (absolute_import, division,
                        print_function)
import os
import subprocess as sp
import threading
import time
import pytest
from ..read_config import read_config
from ..monitor import main_monitor, create_dummy_config, HAS_WATCHDOG
from ..scan import product_path_from_file_name


class TestMonitor(object):
    @classmethod
    def setup_class(klass):
        import os

        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, 'data')
        klass.specdir = os.path.join(klass.datadir, 'spectrum')
        klass.config_file = \
            os.path.abspath(os.path.join(klass.datadir, 'test_config.ini'))

        config = read_config(klass.config_file)

        dummy = os.path.join(klass.datadir, "srt_data_dummy.fits")
        klass.proddir, _ = \
            product_path_from_file_name(dummy, workdir=config['workdir'],
                                        productdir=config['productdir'])

        klass.config = config

        klass.file_empty_init = \
            os.path.abspath(os.path.join(klass.datadir, 'spectrum',
                                         "srt_data.fits"))
        klass.file_empty = os.path.abspath(dummy)
        klass.file_empty_hdf5 = \
            os.path.abspath(os.path.join(klass.datadir,
                                         "srt_data_dummy.hdf5"))
        klass.file_empty_pdf0 = \
            os.path.abspath(os.path.join(klass.datadir,
                                         "srt_data_dummy_0.jpg"))
        klass.file_empty_pdf1 = \
            os.path.abspath(os.path.join(klass.datadir,
                                         "srt_data_dummy_1.jpg"))
        klass.file_empty_hdf5_alt = \
            os.path.abspath(os.path.join(klass.proddir,
                                         "srt_data_dummy.hdf5"))
        klass.file_empty_pdf0_alt = \
            os.path.abspath(os.path.join(klass.proddir,
                                         "srt_data_dummy_0.jpg"))
        klass.file_empty_pdf1_alt = \
            os.path.abspath(os.path.join(klass.proddir,
                                         "srt_data_dummy_1.jpg"))
        klass.file_index = 'index.html'
        klass.dummy_config = 'monitor_config.ini'

        if os.path.exists(klass.file_empty):
            os.unlink(klass.file_empty)
        if os.path.exists(klass.file_empty_pdf0):
            os.unlink(klass.file_empty_pdf0)
        if os.path.exists(klass.file_empty_pdf1):
            os.unlink(klass.file_empty_pdf1)

    def teardown_method(self):
        files = [
            self.file_empty,
            self.file_empty_hdf5,
            self.file_empty_pdf0,
            self.file_empty_pdf1,
            self.file_empty_hdf5_alt,
            self.file_empty_pdf0_alt,
            self.file_empty_pdf1_alt
        ]

        for fname in files:
            if os.path.exists(fname):
                os.unlink(fname)

    @pytest.mark.skipif('not HAS_WATCHDOG')
    def test_monitor_installed(self):
        sp.check_call('SDTmonitor -h'.split())

    @pytest.mark.skipif('not HAS_WATCHDOG')
    def test_all(self):
        def process():
            main_monitor([self.datadir, '--test'])

        w = threading.Thread(name='worker', target=process)
        w.start()
        time.sleep(1)

        sp.check_call('cp {} {}'.format(self.file_empty_init,
                                        self.file_empty).split())

        time.sleep(8)
        w.join()

        for fname in [self.file_empty_pdf0, self.file_empty_pdf1,
                      'latest_0.jpg', 'latest_1.jpg']:
            assert os.path.exists(fname)
            os.unlink(fname)

    @pytest.mark.skipif('not HAS_WATCHDOG')
    def test_all_new_with_config(self):
        fname = self.config_file

        def process():
            main_monitor([self.datadir, '--test', '-c', fname])

        w = threading.Thread(name='worker', target=process)
        w.start()
        time.sleep(1)

        sp.check_call('cp {} {}'.format(self.file_empty_init,
                                        self.file_empty).split())

        time.sleep(8)
        w.join()

        for fname in [self.file_empty_pdf0_alt, self.file_empty_pdf1_alt,
                      'latest_0.jpg', 'latest_1.jpg']:
            assert os.path.exists(fname)
            os.unlink(fname)

    @pytest.mark.skipif('not HAS_WATCHDOG')
    def test_polling(self):
        def process():
            main_monitor([self.datadir, '--test', '--polling'])

        w = threading.Thread(name='worker', target=process)
        w.start()
        time.sleep(1)

        sp.check_call('cp {} {}'.format(self.file_empty_init,
                                        self.file_empty).split())

        time.sleep(8)
        w.join()

        for fname in [self.file_empty_pdf0, self.file_empty_pdf1,
                      'latest_0.jpg', 'latest_1.jpg']:
            assert os.path.exists(fname)
            os.unlink(fname)

    @pytest.mark.skipif('not HAS_WATCHDOG')
    def test_nosave(self):
        def process():
            main_monitor([self.datadir, '--test', '--nosave'])

        w = threading.Thread(name='worker', target=process)
        w.start()
        time.sleep(1)

        sp.check_call('cp {} {}'.format(self.file_empty_init,
                                        self.file_empty).split())

        time.sleep(8)
        w.join()

        for fname in [self.file_empty_pdf0, self.file_empty_pdf1,
                      self.file_empty_hdf5]:
            assert not os.path.exists(fname)

        for fname in ['latest_0.jpg', 'latest_1.jpg']:
            assert os.path.exists(fname)
            os.unlink(fname)


    @pytest.mark.skipif('not HAS_WATCHDOG')
    def test_delete_old_images(self):
        def process():
            main_monitor([self.datadir, '--test'])

        files = [self.file_empty_pdf0, self.file_empty_pdf1]
        files += ['latest_{}.jpg'.format(i) for i in range(8)]

        for fname in files[2:]:
            sp.check_call('touch {}'.format(fname).split())

        w = threading.Thread(name='worker', target=process)
        w.start()
        time.sleep(1)

        sp.check_call('cp {} {}'.format(self.file_empty_init,
                                        self.file_empty).split())

        time.sleep(8)
        w.join()

        for fname in files[4:]:
            assert not os.path.exists(fname)
        for fname in files[:4]:
            assert os.path.exists(fname)
            os.unlink(fname)

    @classmethod
    def teardown_class(klass):
        if os.path.exists(klass.file_index):
            os.unlink(klass.file_index)
        if os.path.exists(klass.dummy_config):
            os.unlink(klass.dummy_config)
