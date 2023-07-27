import os
import subprocess as sp
import threading
import time
import pytest
import urllib
import shutil
import socket
import asyncio
import json
import base64
import numpy as np

from srttools.read_config import read_config

try:
    from srttools.monitor import MAX_FEEDS
    from srttools.monitor import Monitor
    from tornado.websocket import websocket_connect

    HAS_DEPENDENCIES = True
except ImportError:
    HAS_DEPENDENCIES = False

try:
    import pytest_asyncio

    HAS_PYTEST_ASYNCIO = True
except ImportError:
    HAS_PYTEST_ASYNCIO = False


from srttools.scan import product_path_from_file_name
from srttools.utils import look_for_files_or_bust

STANDARD_TIMEOUT = 10


def get_free_tcp_port():
    """This method provides an available TCP port to test
    the WebSocket server without getting a SocketError"""
    tcp = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    tcp.bind(("", 0))
    addr, port = tcp.getsockname()
    tcp.close()
    return port


def compare_images(received_string, image_file):
    image_string = base64.b64encode(open(image_file, "rb").read())
    image_string = image_string.decode("utf-8")
    assert received_string == image_string


class WebSocketClient(object):
    def __init__(self, url):
        self._url = url
        self._ws = None

    async def __aenter__(self):
        self._ws = await websocket_connect(self._url)
        return self

    async def __aexit__(self, exc_type, exc, tb):
        if self._ws:
            self._ws.close()

    async def get_messages(self, n, timeout=None):
        messages = []
        t0 = time.time()
        while True:
            elapsed = time.time() - t0
            if len(messages) == n or (timeout and elapsed >= timeout):
                break
            try:
                remaining = timeout - elapsed if timeout else 0
                messages.append(await asyncio.wait_for(self._ws.read_message(), timeout=remaining))
            except asyncio.TimeoutError:
                pass
        return messages


class TestMonitor(object):
    @classmethod
    def setup_class(klass):
        import os

        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, "data")
        klass.specdir = os.path.join(klass.datadir, "spectrum")
        klass.config_file = os.path.abspath(os.path.join(klass.datadir, "test_config.ini"))

        config = read_config(klass.config_file)

        dummy = os.path.join(klass.datadir, "srt_data_dummy.fits")
        klass.proddir, _ = product_path_from_file_name(
            dummy, workdir=config["workdir"], productdir=config["productdir"]
        )

        klass.config = config

        klass.file_empty_init = os.path.abspath(
            os.path.join(klass.datadir, "spectrum", "srt_data.fits")
        )
        klass.file_empty_init_single_feed = os.path.abspath(
            os.path.join(klass.datadir, "spectrum", "new_sardara.fits5")
        )

        klass.file_empty = os.path.abspath(dummy)
        klass.file_empty_single_feed = os.path.abspath(dummy) + "5"

        klass.file_empty_hdf5 = os.path.abspath(os.path.join(klass.datadir, "srt_data_dummy.hdf5"))
        klass.file_empty_pdf0 = os.path.abspath(os.path.join(klass.datadir, "srt_data_dummy_0.png"))
        klass.file_empty_pdf1 = os.path.abspath(os.path.join(klass.datadir, "srt_data_dummy_1.png"))
        klass.file_empty_hdf5_alt = os.path.abspath(
            os.path.join(klass.proddir, "srt_data_dummy.hdf5")
        )
        klass.file_empty_pdf0_alt = os.path.abspath(
            os.path.join(klass.proddir, "srt_data_dummy_0.png")
        )
        klass.file_empty_pdf1_alt = os.path.abspath(
            os.path.join(klass.proddir, "srt_data_dummy_1.png")
        )
        klass.file_empty_pdf10 = os.path.abspath(
            os.path.join(klass.proddir, "srt_data_dummy5_0.png")
        )
        klass.file_empty_hdf5_SF = os.path.abspath(
            os.path.join(klass.proddir, "srt_data_dummy5.hdf5")
        )
        klass.file_index = "index.html"
        klass.dummy_config = "monitor_config.ini"

        if os.path.exists(klass.file_empty):
            os.unlink(klass.file_empty)
        if os.path.exists(klass.file_empty_single_feed):
            os.unlink(klass.file_empty_single_feed)
        if os.path.exists(klass.file_empty_pdf0):
            os.unlink(klass.file_empty_pdf0)
        if os.path.exists(klass.file_empty_pdf1):
            os.unlink(klass.file_empty_pdf1)

    def setup_method(self):
        self.monitor = None

    def teardown_method(self):
        if self.monitor:
            # This ensures to stop the monitor even when an exception is raised in a test
            self.monitor.stop()
        files = [
            self.file_empty,
            self.file_empty_single_feed,
            self.file_empty_hdf5,
            self.file_empty_pdf0,
            self.file_empty_pdf1,
            self.file_empty_hdf5_alt,
            self.file_empty_pdf0_alt,
            self.file_empty_pdf1_alt,
        ]

        for fname in files:
            if os.path.exists(fname):
                os.unlink(fname)

    @pytest.mark.xfail(strict=False)
    @pytest.mark.skipif("not HAS_DEPENDENCIES")
    def test_monitor_installed(self):
        sp.check_call("SDTmonitor -h".split())

    @pytest.mark.xfail(strict=False)
    @pytest.mark.skipif("not HAS_DEPENDENCIES")
    def test_all(self):
        port = get_free_tcp_port()
        self.monitor = Monitor([self.datadir], port=port)
        self.monitor.start()

        time.sleep(1)

        shutil.copy(self.file_empty_init, self.file_empty)

        files = ["latest_0.png", "latest_1.png"]
        look_for_files_or_bust(files, STANDARD_TIMEOUT)

        for fname in files:
            os.unlink(fname)

    @pytest.mark.xfail(strict=False)
    @pytest.mark.skipif("not HAS_DEPENDENCIES")
    def test_all_new_with_config(self):
        port = get_free_tcp_port()
        self.monitor = Monitor([self.datadir], config_file=self.config_file, port=port)
        self.monitor.start()

        time.sleep(1)

        shutil.copy(self.file_empty_init, self.file_empty)

        files = ["latest_0.png", "latest_1.png"]
        look_for_files_or_bust(files, STANDARD_TIMEOUT)

        for fname in files:
            os.unlink(fname)

    @pytest.mark.xfail(strict=False)
    @pytest.mark.skipif("not HAS_DEPENDENCIES")
    def test_verbose(self):
        port = get_free_tcp_port()
        self.monitor = Monitor(
            [self.datadir],
            verbosity=1,
            config_file=self.config_file,
            port=port,
        )
        self.monitor.start()

        time.sleep(1)

        shutil.copy(self.file_empty_init, self.file_empty)

        files = ["latest_0.png", "latest_1.png"]
        look_for_files_or_bust(files, STANDARD_TIMEOUT)

        for fname in files:
            os.unlink(fname)

    @pytest.mark.xfail(strict=False)
    @pytest.mark.skipif("not HAS_DEPENDENCIES")
    def test_a_single_feed(self):
        port = get_free_tcp_port()
        self.monitor = Monitor([self.datadir], config_file=self.config_file, port=port)
        self.monitor.start()

        time.sleep(1)

        shutil.copy(self.file_empty_init_single_feed, self.file_empty_single_feed)

        files = ["latest_10.png"]
        look_for_files_or_bust(files, STANDARD_TIMEOUT)

        for fname in files:
            os.unlink(fname)

    @pytest.mark.xfail(strict=False)
    @pytest.mark.skipif("not HAS_DEPENDENCIES")
    def test_polling(self):
        port = get_free_tcp_port()
        self.monitor = Monitor([self.datadir], polling=True, port=port)
        self.monitor.start()

        time.sleep(1)

        shutil.copy(self.file_empty_init, self.file_empty)

        files = ["latest_0.png", "latest_1.png"]
        look_for_files_or_bust(files, STANDARD_TIMEOUT)

        for fname in files:
            os.unlink(fname)

    @pytest.mark.xfail(strict=False)
    @pytest.mark.skipif("not HAS_DEPENDENCIES")
    def test_workers(self):
        files = ["latest_8.png", "latest_10.png"]

        port = get_free_tcp_port()
        self.monitor = Monitor([self.datadir], config_file=self.config_file, workers=2, port=port)
        self.monitor.start()

        time.sleep(1)

        shutil.copy(self.file_empty_init_single_feed, self.file_empty_single_feed)
        shutil.copy(
            self.file_empty_init_single_feed,
            self.file_empty_single_feed.replace("fits5", "fits4"),
        )

        look_for_files_or_bust(files, STANDARD_TIMEOUT)

        for i, fname in enumerate(files):
            os.unlink(fname)

        os.unlink(self.file_empty_single_feed.replace("fits5", "fits4"))

    @pytest.mark.xfail(strict=False)
    @pytest.mark.skipif("not HAS_DEPENDENCIES")
    def test_delete_old_images(self):
        files = ["latest_{}.png".format(i) for i in range(8)]

        for fname in files[2:]:
            sp.check_call("touch {}".format(fname).split())

        port = get_free_tcp_port()
        self.monitor = Monitor([self.datadir], port=port)
        self.monitor.start()

        time.sleep(1)

        shutil.copy(self.file_empty_init, self.file_empty)

        look_for_files_or_bust(files[:2], STANDARD_TIMEOUT)

        for fname in files[2:]:
            assert not os.path.exists(fname)

        for fname in files[:2]:
            os.unlink(fname)

    @pytest.mark.xfail(strict=False)
    @pytest.mark.skipif("not HAS_DEPENDENCIES")
    def test_http_server(self):
        port = get_free_tcp_port()
        self.monitor = Monitor([self.datadir], port=port)
        self.monitor.start()

        time.sleep(1)

        shutil.copy(self.file_empty_init, self.file_empty)

        # Check that index.html is provided by the HTTP server
        url = "http://127.0.0.1:{}/".format(port)
        r = urllib.request.urlopen(url, timeout=5)
        assert r.code == 200

        files = ["latest_0.png", "latest_1.png"]
        look_for_files_or_bust(files, STANDARD_TIMEOUT)

        for fname in files:
            assert os.path.exists(fname)
            os.unlink(fname)

    @pytest.mark.xfail(strict=False)
    @pytest.mark.asyncio
    @pytest.mark.skipif("not HAS_DEPENDENCIES or not HAS_PYTEST_ASYNCIO")
    async def test_websocket_server(self):
        port = get_free_tcp_port()
        self.monitor = Monitor([self.datadir], port=port)
        self.monitor.start()

        time.sleep(1)

        ws_url = "ws://localhost:{}/images".format(port)
        async with WebSocketClient(ws_url) as ws:
            # Retrieve the starting images, they should be 14 blanks,
            # we ask for 15 and make sure they are actually 14
            l = await ws.get_messages(15, timeout=1)
            assert len(l) == 14
            for image_string in l:
                image = json.loads(image_string)
                assert image["image"] == ""

            # Now trigger the process of a file
            shutil.copy(self.file_empty_init, self.file_empty)

            files = ["latest_0.png", "latest_1.png"]
            look_for_files_or_bust(files, STANDARD_TIMEOUT)

            # Ask the new images
            l = await ws.get_messages(2, timeout=1)
            assert len(l) == 2
            for image_string in l:
                image = json.loads(image_string)
                assert image["index"] in [0, 1]
                compare_images(image["image"], "latest_{}.png".format(image["index"]))

            for fname in files:
                assert os.path.exists(fname)
                os.unlink(fname)

    @pytest.mark.xfail(strict=False)
    @pytest.mark.skipif("HAS_DEPENDENCIES")
    def test_dependencies_missing(self):
        with pytest.warns(UserWarning) as record:
            from srttools.monitor import Monitor
        at_least_one_warning = False
        for string in ["watchdog", "tornado"]:
            at_least_one_warning = at_least_one_warning or np.any(
                [f"install {string}" in r.message.args[0] for r in record]
            )
        assert at_least_one_warning

    @classmethod
    def teardown_class(klass):
        if os.path.exists(klass.file_index):
            os.unlink(klass.file_index)
        if os.path.exists(klass.dummy_config):
            os.unlink(klass.dummy_config)
