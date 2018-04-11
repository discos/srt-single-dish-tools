import os
import subprocess as sp
import threading
import time
import pytest
from ..monitor import main_monitor, create_dummy_config, HAS_WATCHDOG


class TestMonitor(object):
    @classmethod
    def setup_class(klass):
        import os

        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, 'data')
        klass.specdir = os.path.join(klass.datadir, 'spectrum')

        klass.file_empty_init = \
            os.path.abspath(os.path.join(klass.datadir, 'spectrum',
                                         "srt_data.fits"))
        klass.file_empty = \
            os.path.abspath(os.path.join(klass.datadir,
                                         "srt_data_dummy.fits"))
        klass.file_empty_hdf5 = \
            os.path.abspath(os.path.join(klass.datadir,
                                         "srt_data_dummy.hdf5"))
        klass.file_empty_pdf0 = \
            os.path.abspath(os.path.join(klass.datadir,
                                         "srt_data_dummy_0.jpg"))
        klass.file_empty_pdf1 = \
            os.path.abspath(os.path.join(klass.datadir,
                                         "srt_data_dummy_1.jpg"))
        if os.path.exists(klass.file_empty):
            os.unlink(klass.file_empty)
        if os.path.exists(klass.file_empty_pdf0):
            os.unlink(klass.file_empty_pdf0)
        if os.path.exists(klass.file_empty_pdf1):
            os.unlink(klass.file_empty_pdf1)

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
        fname = create_dummy_config()

        def process():
            main_monitor([self.datadir, '--test', '-c', fname])

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

    @classmethod
    def teardown_class(klass):
        if os.path.exists(klass.file_empty):
            os.unlink(klass.file_empty)
        if os.path.exists(klass.file_empty_hdf5):
            os.unlink(klass.file_empty_hdf5)
