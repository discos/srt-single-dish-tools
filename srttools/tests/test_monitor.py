
import os
import subprocess as sp
import threading
import time
import pytest
import urllib
import shutil

from ..read_config import read_config
from ..monitor import main_monitor, create_dummy_config, HAS_WATCHDOG
from ..scan import product_path_from_file_name
from ..utils import look_for_files_or_bust


STANDARD_TIMEOUT = 10


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
        klass.file_empty_init_single_feed = \
            os.path.abspath(os.path.join(klass.datadir, 'spectrum',
                                         "new_sardara.fits5"))

        klass.file_empty = os.path.abspath(dummy)
        klass.file_empty_single_feed = os.path.abspath(dummy) + '5'

        klass.file_empty_hdf5 = \
            os.path.abspath(os.path.join(klass.datadir,
                                         "srt_data_dummy.hdf5"))
        klass.file_empty_pdf0 = \
            os.path.abspath(os.path.join(klass.datadir,
                                         "srt_data_dummy_0.png"))
        klass.file_empty_pdf1 = \
            os.path.abspath(os.path.join(klass.datadir,
                                         "srt_data_dummy_1.png"))
        klass.file_empty_hdf5_alt = \
            os.path.abspath(os.path.join(klass.proddir,
                                         "srt_data_dummy.hdf5"))
        klass.file_empty_pdf0_alt = \
            os.path.abspath(os.path.join(klass.proddir,
                                         "srt_data_dummy_0.png"))
        klass.file_empty_pdf1_alt = \
            os.path.abspath(os.path.join(klass.proddir,
                                         "srt_data_dummy_1.png"))
        klass.file_empty_pdf10 = \
            os.path.abspath(os.path.join(klass.proddir,
                                         "srt_data_dummy5_0.png"))
        klass.file_empty_hdf5_SF = \
            os.path.abspath(os.path.join(klass.proddir,
                                         "srt_data_dummy5.hdf5"))
        klass.file_index = 'index.html'
        klass.dummy_config = 'monitor_config.ini'

        if os.path.exists(klass.file_empty):
            os.unlink(klass.file_empty)
        if os.path.exists(klass.file_empty_single_feed):
            os.unlink(klass.file_empty_single_feed)
        if os.path.exists(klass.file_empty_pdf0):
            os.unlink(klass.file_empty_pdf0)
        if os.path.exists(klass.file_empty_pdf1):
            os.unlink(klass.file_empty_pdf1)

    def teardown_method(self):
        files = [
            self.file_empty,
            self.file_empty_single_feed,
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

        shutil.copy(self.file_empty_init, self.file_empty)

        files = [self.file_empty_pdf0, self.file_empty_pdf1,
                 'latest_0.png', 'latest_1.png']
        look_for_files_or_bust(files, STANDARD_TIMEOUT)

        w.join()

        for fname in files:
            os.unlink(fname)

    @pytest.mark.skipif('not HAS_WATCHDOG')
    def test_all_new_with_config(self):
        fname = self.config_file

        def process():
            main_monitor([self.datadir, '--test', '-c', fname])

        w = threading.Thread(name='worker', target=process)
        w.start()
        time.sleep(1)

        shutil.copy(self.file_empty_init, self.file_empty)

        files = [self.file_empty_pdf0_alt, self.file_empty_pdf1_alt,
                      'latest_0.png', 'latest_1.png']
        look_for_files_or_bust(files, STANDARD_TIMEOUT)

        w.join()

        for fname in files:
            os.unlink(fname)

    @pytest.mark.skipif('not HAS_WATCHDOG')
    def test_verbose(self):
        fname = self.config_file

        def process():
            main_monitor([self.datadir, '--test', '-v', '-c', fname])

        w = threading.Thread(name='worker', target=process)
        w.start()
        time.sleep(1)

        shutil.copy(self.file_empty_init, self.file_empty)

        files = [self.file_empty_pdf0_alt, self.file_empty_pdf1_alt,
                      'latest_0.png', 'latest_1.png']
        look_for_files_or_bust(files, STANDARD_TIMEOUT)

        w.join()

        for fname in files:
            os.unlink(fname)

    @pytest.mark.skipif('not HAS_WATCHDOG')
    def test_a_single_feed(self):
        fname = self.config_file

        def process():
            main_monitor([self.datadir, '--test', '-c', fname])

        w = threading.Thread(name='worker', target=process)
        w.start()
        time.sleep(1)

        shutil.copy(self.file_empty_init_single_feed,
                    self.file_empty_single_feed)

        files = [self.file_empty_pdf10, self.file_empty_hdf5_SF,
                 'latest_10.png']
        look_for_files_or_bust(files, STANDARD_TIMEOUT)

        w.join()

        for fname in files:
            os.unlink(fname)

    @pytest.mark.skipif('not HAS_WATCHDOG')
    def test_polling(self):
        def process():
            main_monitor([self.datadir, '--test', '--polling'])

        w = threading.Thread(name='worker', target=process)
        w.start()
        time.sleep(1)

        shutil.copy(self.file_empty_init, self.file_empty)

        files = [self.file_empty_pdf0, self.file_empty_pdf1,
                 'latest_0.png', 'latest_1.png']
        look_for_files_or_bust(files, STANDARD_TIMEOUT)

        w.join()

        for fname in files:
            os.unlink(fname)

    @pytest.mark.skipif('not HAS_WATCHDOG')
    def test_nosave(self):
        def process():
            main_monitor([self.datadir, '--test', '--nosave'])

        w = threading.Thread(name='worker', target=process)
        w.start()
        time.sleep(1)

        shutil.copy(self.file_empty_init, self.file_empty)

        files = ['latest_0.png', 'latest_1.png']
        look_for_files_or_bust(files, STANDARD_TIMEOUT)

        w.join()

        for fname in files:
            os.unlink(fname)

        for fname in [self.file_empty_pdf0, self.file_empty_pdf1,
                      self.file_empty_hdf5]:
            assert not os.path.exists(fname)

    @pytest.mark.skipif('not HAS_WATCHDOG')
    def test_workers(self):
        fname = self.config_file

        files = [
            self.file_empty_pdf10,
            self.file_empty_pdf10.replace('dummy5', 'dummy4'),
            self.file_empty_hdf5_SF,
            self.file_empty_hdf5_SF.replace('dummy5', 'dummy4'),
            'latest_8.png',
            'latest_10.png'
        ]

        def process():
            main_monitor([self.datadir, '--test', '-c', fname,
                          '--timeout', '30',
                          '--workers', '2'])

        w = threading.Thread(name='worker', target=process)
        w.start()
        time.sleep(1)

        shutil.copy(
            self.file_empty_init_single_feed,
            self.file_empty_single_feed
        )
        shutil.copy(
            self.file_empty_init_single_feed,
            self.file_empty_single_feed.replace('fits5', 'fits4')
        )

        look_for_files_or_bust(files, 30)

        w.join()

        for i, fname in enumerate(files):
            os.unlink(fname)

        os.unlink(self.file_empty_single_feed.replace('fits5', 'fits4'))

    @pytest.mark.skipif('not HAS_WATCHDOG')
    def test_delete_old_images(self):
        def process():
            main_monitor([self.datadir, '--test'])

        files = [self.file_empty_pdf0, self.file_empty_pdf1]
        files += ['latest_{}.png'.format(i) for i in range(8)]

        for fname in files[2:]:
            sp.check_call('touch {}'.format(fname).split())

        w = threading.Thread(name='worker', target=process)
        w.start()
        time.sleep(1)

        shutil.copy(self.file_empty_init, self.file_empty)

        time.sleep(STANDARD_TIMEOUT)
        w.join()

        for fname in files[4:]:
            assert not os.path.exists(fname)
        for fname in files[:4]:
            assert os.path.exists(fname)
            os.unlink(fname)

    @pytest.mark.skipif('not HAS_WATCHDOG')
    def test_http_server_found(self):
        def process():
            main_monitor([self.datadir, '--test',
                          '--http-server-port', '10000', '--timeout', '30'])

        w = threading.Thread(name='worker', target=process)
        w.start()
        time.sleep(1)

        shutil.copy(self.file_empty_init, self.file_empty)

        time.sleep(7)
        urls = [
            '',
            'index.html',
            'latest_0.png',
            'latest_1.png'#,
            # 'latest_0.png?%d' % int(time.time())
        ]
        # Check that the files are provided by the HTTP server
        for url in urls:
            url = 'http://127.0.0.1:10000/' + url
            print(url)
            # import socket
            # try:
            r = urllib.request.urlopen(url, timeout=5)
            assert r.code == 200
            # except socket.timeout as e:
            #     w.join()
            #     raise

        print("Done")

        w.join()

        for fname in [self.file_empty_pdf0, self.file_empty_pdf1,
                      'latest_0.png', 'latest_1.png']:
            assert os.path.exists(fname)
            os.unlink(fname)

    @pytest.mark.skipif('not HAS_WATCHDOG')
    def test_http_server_not_found(self):
        def process():
            main_monitor([self.datadir, '--test', '--http-server-port', '10000'])

        w = threading.Thread(name='worker', target=process)
        w.start()
        time.sleep(1)

        shutil.copy(self.file_empty_init, self.file_empty)

        time.sleep(7)

        # Ask for non-existent files
        for i in [20, 21]:
            url = 'http://127.0.0.1:10000/latest_{}.png'.format(i)
            with pytest.raises(urllib.error.HTTPError):
                r = urllib.request.urlopen(url, timeout=10)

        w.join()

        for fname in [self.file_empty_pdf0, self.file_empty_pdf1,
                      'latest_0.png', 'latest_1.png']:
            assert os.path.exists(fname)
            os.unlink(fname)

    @pytest.mark.skipif('not HAS_WATCHDOG')
    def test_http_server_forbidden(self):
        def process():
            main_monitor([self.datadir, '--test', '--http-server-port', '10000'])

        w = threading.Thread(name='worker', target=process)
        w.start()
        time.sleep(1)

        shutil.copy(self.file_empty_init, self.file_empty)

        time.sleep(7)

        # Create a dummy 'latest_0.txt' file.
        # The file should not be accessible since its name does not match
        # 'index.html' nor 'index.html' nor 'latest_X.EXT' pattern, with EXT
        # equal to the configured file extension, in this case is 'png'
        forbidden_file = 'latest_0.txt'
        open(forbidden_file, 'w').close()
        url = 'http://127.0.0.1:10000/' + forbidden_file
        with pytest.raises(urllib.error.HTTPError):
            r = urllib.request.urlopen(url, timeout=10)

        w.join()

        for fname in [self.file_empty_pdf0, self.file_empty_pdf1,
                      'latest_0.png', 'latest_1.png', forbidden_file]:
            assert os.path.exists(fname)
            os.unlink(fname)

    @classmethod
    def teardown_class(klass):
        if os.path.exists(klass.file_index):
            os.unlink(klass.file_index)
        if os.path.exists(klass.dummy_config):
            os.unlink(klass.dummy_config)
