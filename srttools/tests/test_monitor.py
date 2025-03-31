import multiprocessing as mp
import queue
import base64
import json
import os
import shutil
import socket
import subprocess as sp
import threading
import time
import urllib
import glob
import numpy as np
import pytest
import copy

from srttools.read_config import read_config
from srttools.monitor import main_monitor  # import the CLI regardless of dependencies

try:
    from tornado.websocket import websocket_connect
    from tornado.ioloop import IOLoop
    from srttools.monitor import Monitor
    from srttools.monitor.common import stop_event, MAX_PROCS

    HAS_DEPENDENCIES = True
except ImportError:
    HAS_DEPENDENCIES = False


from srttools.scan import product_path_from_file_name
from srttools.utils import look_for_files_or_bust


STANDARD_TIMEOUT = 20


@pytest.fixture(scope="session", autouse=True)
def always_spawn():
    mp.set_start_method("spawn", force=True)


def get_free_tcp_port():
    """This method provides an available TCP port to test
    the WebSocket server without getting a SocketError"""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as tcp:
        tcp.bind(("localhost", 0))
        _, port = tcp.getsockname()
        return port


def wait_for_message(caplog, message, clear=False, timeout=STANDARD_TIMEOUT):
    expiry = time.time() + timeout
    if clear:
        caplog.clear()
    while time.time() < expiry:
        if message in caplog.text:
            return True
        time.sleep(0.01)
    return False


def compare_images(received_string, image_file):
    image_string = base64.b64encode(open(image_file, "rb").read())
    image_string = image_string.decode("utf-8")
    assert received_string == image_string


class WebSocketClient:
    def __init__(self, url):
        self.url = url
        self.messages = queue.Queue()
        self.connection = None
        IOLoop.current().add_callback(self.connect)

    def __del__(self):
        self.close()

    def close(self):
        self.connection.close()

    def connect(self):
        websocket_connect(self.url, callback=self.on_connect)

    def on_connect(self, future):
        try:
            self.connection = future.result()
            self.connection.read_message(self.on_message)
        except Exception as e:
            print(e)

    def on_message(self, message):
        if message is None:
            return
        else:
            self.messages.put(message.result())
        self.connection.read_message(self.on_message)

    def get_messages(self, n, timeout=None):
        t0 = time.time()
        messages = []
        while True:
            elapsed = time.time() - t0
            try:
                messages.append(self.messages.get(timeout=0.01))
            except queue.Empty:
                continue
            if len(messages) == n or (timeout and elapsed >= timeout):
                break
        return messages


@pytest.mark.skipif("not HAS_DEPENDENCIES")
class TestMonitor:
    @classmethod
    def setup_class(klass):
        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, "data")
        klass.specdir = os.path.join(klass.datadir, "spectrum")
        klass.config_file = os.path.abspath(os.path.join(klass.datadir, "test_config.ini"))

        klass.config = read_config(klass.config_file)

        dummy = os.path.join(klass.datadir, "srt_data_dummy.fits")
        klass.proddir, _ = product_path_from_file_name(
            dummy, workdir=klass.config["workdir"], productdir=klass.config["productdir"]
        )

        klass.file_empty_init = os.path.abspath(
            os.path.join(klass.datadir, "spectrum", "srt_data.fits")
        )
        klass.file_empty_init_single_feed = os.path.abspath(
            os.path.join(klass.datadir, "spectrum", "new_sardara.fits5")
        )

        klass.file_empty = os.path.abspath(dummy)
        klass.file_empty_single_feed = os.path.abspath(dummy) + "5"

        klass.file_empty_hdf5 = os.path.abspath(os.path.join(klass.datadir, "srt_data_dummy.hdf5"))
        klass.file_empty_pdf0 = os.path.abspath(
            os.path.join(klass.datadir, "srt_data_dummy_000.png")
        )
        klass.file_empty_pdf1 = os.path.abspath(
            os.path.join(klass.datadir, "srt_data_dummy_001.png")
        )
        klass.file_empty_hdf5_alt = os.path.abspath(
            os.path.join(klass.proddir, "srt_data_dummy.hdf5")
        )
        klass.file_empty_pdf0_alt = os.path.abspath(
            os.path.join(klass.proddir, "srt_data_dummy_000.png")
        )
        klass.file_empty_pdf1_alt = os.path.abspath(
            os.path.join(klass.proddir, "srt_data_dummy_001.png")
        )
        klass.file_empty_pdf10 = os.path.abspath(
            os.path.join(klass.proddir, "srt_data_dummy5_000.png")
        )
        klass.file_empty_hdf5_SF = os.path.abspath(
            os.path.join(klass.proddir, "srt_data_dummy5.hdf5")
        )
        klass.file_index = "index.html"
        klass.dummy_config = "monitor_config.ini"

        klass.removefiles = [
            klass.file_empty,
            klass.file_empty_single_feed,
            klass.file_empty_hdf5,
            klass.file_empty_pdf0,
            klass.file_empty_pdf1,
            klass.file_empty_hdf5_alt,
            klass.file_empty_pdf0_alt,
            klass.file_empty_pdf1_alt,
            klass.file_empty_pdf10,
            klass.file_empty_hdf5_SF,
        ]

    @classmethod
    def teardown_class(klass):
        if os.path.exists(klass.file_index):
            os.unlink(klass.file_index)
        if os.path.exists(klass.dummy_config):
            os.unlink(klass.dummy_config)

    def setup_method(self):
        self.monitor = None

        removeimages = glob.glob("latest_*.png")
        for fname in self.removefiles + removeimages:
            if os.path.exists(fname):
                os.unlink(fname)

    def teardown_method(self):
        if self.monitor:
            # This ensures to stop the monitor even when an exception is raised in a test
            self.monitor.stop()

        removeimages = glob.glob("latest_*.png")
        for fname in self.removefiles + removeimages:
            if os.path.exists(fname):
                os.unlink(fname)

    def test_monitor_installed(self):
        sp.check_call("SDTmonitor -h".split())

    def test_stop_without_start(self):
        port = get_free_tcp_port()
        self.monitor = Monitor([self.datadir], localhost=True, port=port)
        self.monitor.stop()

    def test_start_twice(self):
        port = get_free_tcp_port()
        self.monitor = Monitor([self.datadir], localhost=True, port=port)
        self.monitor.start()
        time.sleep(1)
        self.monitor.start()

    def test_all(self, caplog):
        port = get_free_tcp_port()
        self.monitor = Monitor([self.datadir], localhost=True, port=port)
        self.monitor.start()

        time.sleep(1)

        shutil.copy(self.file_empty_init, self.file_empty)
        wait_for_message(caplog, f"Loading file {self.file_empty}") or pytest.fail(
            "Processing not started!"
        )
        files = ["latest_000.png", "latest_001.png"]
        wait_for_message(caplog, f"Completed file {self.file_empty}") or pytest.fail(
            "Processing not completed!"
        )
        look_for_files_or_bust(files, STANDARD_TIMEOUT)

    def test_all_new_with_config(self, caplog):
        port = get_free_tcp_port()
        self.monitor = Monitor(
            [self.datadir], config_file=self.config_file, localhost=True, port=port
        )
        self.monitor.start()

        time.sleep(1)

        shutil.copy(self.file_empty_init, self.file_empty)
        wait_for_message(caplog, f"Loading file {self.file_empty}") or pytest.fail(
            "Processing not started!"
        )
        files = ["latest_000.png", "latest_001.png"]
        wait_for_message(caplog, f"Completed file {self.file_empty}") or pytest.fail(
            "Processing not completed!"
        )
        look_for_files_or_bust(files, STANDARD_TIMEOUT)

    def test_verbose(self, caplog):
        port = get_free_tcp_port()
        self.monitor = Monitor(
            [self.datadir],
            verbosity=1,
            config_file=self.config_file,
            localhost=True,
            port=port,
        )
        self.monitor.start()

        time.sleep(1)

        shutil.copy(self.file_empty_init, self.file_empty)
        wait_for_message(caplog, f"Loading file {self.file_empty}") or pytest.fail(
            "Processing not started!"
        )
        files = ["latest_000.png", "latest_001.png"]
        wait_for_message(caplog, f"Completed file {self.file_empty}") or pytest.fail(
            "Processing not completed!"
        )
        look_for_files_or_bust(files, STANDARD_TIMEOUT)

    def test_fake_file(self, caplog):
        # We just create a fake file, it should be aborted
        port = get_free_tcp_port()
        self.monitor = Monitor(
            [self.datadir], config_file=self.config_file, localhost=True, port=port
        )
        self.monitor.start()

        time.sleep(1)
        shutil.copy(self.config_file, self.file_empty)
        wait_for_message(caplog, f"Loading file {self.file_empty}") or pytest.fail(
            "Processing not started!"
        )
        wait_for_message(caplog, f"Aborted file {self.file_empty}") or pytest.fail(
            "Processing not aborted!"
        )

        time.sleep(2)

    def test_a_single_feed(self, caplog):
        port = get_free_tcp_port()
        self.monitor = Monitor(
            [self.datadir], config_file=self.config_file, localhost=True, port=port
        )
        self.monitor.start()

        time.sleep(1)

        shutil.copy(self.file_empty_init_single_feed, self.file_empty_single_feed)
        wait_for_message(caplog, f"Loading file {self.file_empty_single_feed}") or pytest.fail(
            "Processing not started!"
        )
        files = ["latest_010.png"]
        wait_for_message(caplog, f"Completed file {self.file_empty_single_feed}") or pytest.fail(
            "Processing not completed!"
        )
        look_for_files_or_bust(files, STANDARD_TIMEOUT)

    def test_polling(self, caplog):
        port = get_free_tcp_port()
        self.monitor = Monitor(
            [self.datadir], config_file=self.config_file, polling=True, localhost=True, port=port
        )
        self.monitor.start()

        time.sleep(1)

        shutil.copy(self.file_empty_init, self.file_empty)
        wait_for_message(caplog, f"Loading file {self.file_empty}") or pytest.fail(
            "Processing not started!"
        )
        files = ["latest_000.png", "latest_001.png"]
        wait_for_message(caplog, f"Completed file {self.file_empty}") or pytest.fail(
            "Processing not completed!"
        )
        look_for_files_or_bust(files, STANDARD_TIMEOUT)

    def test_filemoved_event(self, caplog):
        # First we duplicate the file we want to move
        filename_init = self.file_empty_init_single_feed.replace("fits5", "fits4")

        shutil.copy(self.file_empty_init_single_feed, filename_init)

        port = get_free_tcp_port()
        self.monitor = Monitor(
            [self.datadir], config_file=self.config_file, localhost=True, port=port
        )
        self.monitor.start()

        time.sleep(1)

        # Now we move the file
        shutil.move(filename_init, self.file_empty_single_feed)
        wait_for_message(caplog, f"Loading file {self.file_empty_single_feed}") or pytest.fail(
            "Processing not started!"
        )
        files = ["latest_010.png"]
        wait_for_message(caplog, f"Completed file {self.file_empty_single_feed}") or pytest.fail(
            "Processing not completed!"
        )
        look_for_files_or_bust(files, STANDARD_TIMEOUT)

        self.removefiles.append(filename_init)

    def test_filemoved_event_source(self, caplog):
        # First we duplicate the file we want to move
        filename_init = self.file_empty_init_single_feed.replace("fits5", "fits4")
        filename = self.file_empty_single_feed.replace("fits5", "notafits")
        shutil.copy(self.file_empty_init_single_feed, filename_init)

        port = get_free_tcp_port()
        self.monitor = Monitor(
            [self.datadir], config_file=self.config_file, localhost=True, port=port
        )
        self.monitor.start()

        time.sleep(1)

        # Now we move the file changing its extension
        shutil.move(filename_init, self.file_empty_single_feed.replace("fits5", "notafits"))
        not wait_for_message(caplog, "Loading", timeout=5) or pytest.fail(
            "Processing started erroneously!"
        )

        time.sleep(1)

        self.removefiles.append(filename_init)
        self.removefiles.append(filename)

    def test_copy_same_file(self, caplog):
        port = get_free_tcp_port()
        self.monitor = Monitor([self.datadir], localhost=True, port=port)
        self.monitor.start()

        time.sleep(1)

        shutil.copy(self.file_empty_init, self.file_empty)
        wait_for_message(caplog, f"Loading file {self.file_empty}") or pytest.fail(
            "Processing not started!"
        )

        time.sleep(1.01)

        # Copy the same file
        shutil.copy(self.file_empty_init, self.file_empty)
        not wait_for_message(
            caplog, f"Loading file {self.file_empty}", clear=True, timeout=5
        ) or pytest.fail("Processing started twice!")
        files = ["latest_000.png", "latest_001.png"]
        wait_for_message(caplog, f"Completed file {self.file_empty}") or pytest.fail(
            "Processing not completed!"
        )
        look_for_files_or_bust(files, STANDARD_TIMEOUT)

    def test_workers(self, caplog):
        port = get_free_tcp_port()
        self.monitor = Monitor(
            [self.datadir], config_file=self.config_file, workers=2, localhost=True, port=port
        )
        self.monitor.start()

        time.sleep(1)

        file2 = self.file_empty_single_feed.replace("fits5", "fits4")

        shutil.copy(self.file_empty_init_single_feed, self.file_empty_single_feed)
        shutil.copy(self.file_empty_init_single_feed, file2)
        wait_for_message(caplog, f"Loading file {self.file_empty_single_feed}") or pytest.fail(
            "Processing not started!"
        )
        wait_for_message(caplog, f"Loading file {file2}") or pytest.fail("Processing not started!")
        files = ["latest_008.png", "latest_010.png"]
        wait_for_message(caplog, f"Completed file {self.file_empty_single_feed}") or pytest.fail(
            "Processing not completed!"
        )
        wait_for_message(caplog, f"Completed file {file2}") or pytest.fail(
            "Processing not completed!"
        )
        look_for_files_or_bust(files, STANDARD_TIMEOUT)

        self.removefiles += [file2]

    def test_delete_old_images(self, caplog):
        files = [f"latest_{i:03d}.png" for i in range(8)]

        for fname in files[2:]:
            with open(fname, "w") as f:
                pass

        port = get_free_tcp_port()
        self.monitor = Monitor([self.datadir], localhost=True, port=port)
        self.monitor.start()

        time.sleep(1)

        shutil.copy(self.file_empty_init, self.file_empty)
        wait_for_message(caplog, f"Loading file {self.file_empty}") or pytest.fail(
            "Processing not started!"
        )
        wait_for_message(caplog, f"Completed file {self.file_empty}") or pytest.fail(
            "Processing not completed!"
        )
        look_for_files_or_bust(files[:2], STANDARD_TIMEOUT)

        for fname in files[2:]:
            assert not os.path.exists(fname)

    def test_cancel_processing(self, caplog):
        port = get_free_tcp_port()
        self.monitor = Monitor([self.datadir], localhost=True, port=port)
        self.monitor.start()

        time.sleep(1)

        shutil.copy(self.file_empty_init, self.file_empty)
        wait_for_message(caplog, f"Loading file {self.file_empty}") or pytest.fail(
            "Processing not started!"
        )
        # Immediately exit, should cancel the processing
        # We can only verify this test from coverage

    def test_stop_while_enqueued(self, caplog):
        port = get_free_tcp_port()
        self.monitor = Monitor([self.datadir], localhost=True, port=port, workers=1)
        self.monitor.start()

        time.sleep(1)

        processed_files = []
        fname = self.file_empty.replace(".fits", f"0.fits")
        processed_files.append(fname)
        shutil.copy(self.file_empty_init, fname)
        self.removefiles.append(fname)
        # Wait for the first file to start being processed
        wait_for_message(caplog, f"Loading file {processed_files[0]}") or pytest.fail(
            "Processing not started!"
        )

        # Enqueue some other files
        for index in range(1, 7):
            fname = self.file_empty.replace(".fits", f"{index}.fits")
            processed_files.append(fname)
            shutil.copy(self.file_empty_init, fname)
            self.removefiles.append(fname)
            # Any other file processing should not be started yet
            not wait_for_message(caplog, f"Loading file {fname}", timeout=0.1) or pytest.fail(
                "Processing started anyway!"
            )

        # Now exit immediately and cancel the enqueued processes
        self.removefiles += processed_files

    def test_http_server(self, caplog):
        port = get_free_tcp_port()
        self.monitor = Monitor([self.datadir], localhost=True, port=port)
        self.monitor.start()

        time.sleep(1)

        shutil.copy(self.file_empty_init, self.file_empty)
        wait_for_message(caplog, f"Loading file {self.file_empty}") or pytest.fail(
            "Processing not started!"
        )

        # Check that index.html is provided by the HTTP server
        url = f"http://127.0.0.1:{port}/"
        r = urllib.request.urlopen(url, timeout=5)
        assert r.code == 200

        files = ["latest_000.png", "latest_001.png"]
        wait_for_message(caplog, f"Completed file {self.file_empty}") or pytest.fail(
            "Processing not completed!"
        )
        look_for_files_or_bust(files, STANDARD_TIMEOUT)

    def test_favicon(self):
        port = get_free_tcp_port()
        self.monitor = Monitor([self.datadir], localhost=True, port=port)
        self.monitor.start()

        time.sleep(1)

        # Check that the favicon is provided by the HTTP server
        url = f"http://127.0.0.1:{port}/favicon.ico"
        r = urllib.request.urlopen(url, timeout=5)
        assert r.code == 200

    def test_websocket_server(self, caplog):
        port = get_free_tcp_port()
        self.monitor = Monitor([self.datadir], localhost=True, port=port)
        self.monitor.start()

        time.sleep(1)

        dummy_image = "dummy"
        # We override the images of the monitor
        self.monitor._web_server.images[0] = dummy_image
        self.monitor._web_server.images[1] = dummy_image

        ws_url = f"ws://localhost:{port}/images"
        ws = WebSocketClient(ws_url)

        # Ask the first images
        l = ws.get_messages(2, timeout=1)
        assert len(l) == 2
        for image_string in l:
            image = json.loads(image_string)
            assert image["index"] in [0, 1]
            assert image["image"] == dummy_image

        # Trigger the process of a file
        shutil.copy(self.file_empty_init, self.file_empty)
        wait_for_message(caplog, f"Loading file {self.file_empty}") or pytest.fail(
            "Processing not started!"
        )
        files = ["latest_000.png", "latest_001.png"]
        wait_for_message(caplog, f"Completed file {self.file_empty}") or pytest.fail(
            "Processing not completed!"
        )
        look_for_files_or_bust(files, STANDARD_TIMEOUT)

        # Ask the new images
        l = ws.get_messages(2, timeout=1)
        assert len(l) == 2
        for image_string in l:
            image = json.loads(image_string)
            assert image["index"] in [0, 1]
            compare_images(image["image"], f"latest_{image['index']:03d}.png")
        ws.close()

    def test_busy_port(self):
        port = get_free_tcp_port()
        # First occupy the port with a random socket
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.bind(("localhost", port))
        with pytest.raises(OSError) as exc:
            Monitor([self.datadir], localhost=True, port=port)
        assert f"Port {port} is already being used, choose a different one!" in str(exc.value)
        s.close()

    def test_cli(self):
        port = get_free_tcp_port()
        args = ["-p", f"{port}", "--localhost", "-w", "1", "-c", self.config_file, self.datadir]
        t = threading.Thread(target=main_monitor, args=(args,))
        t.start()
        time.sleep(1)
        stop_event.set()
        t.join()

    def test_cli_busy_port(self, capfd):
        port = get_free_tcp_port()
        # First occupy the port with a random socket
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.bind(("localhost", port))
        args = ["-p", f"{port}", "--localhost", self.datadir]
        with pytest.raises(SystemExit) as exc:
            main_monitor(args)
        assert exc.value.code == 2
        _, err = capfd.readouterr()
        assert f"Port {port} is already being used, choose a different one!" in err
        s.close()

    def test_cli_random_port(self, capfd):
        # A port equal to 0 would cause to select a random port.
        # We intercept it and notify the error
        port = 0
        args = ["-p", f"{port}", "--localhost", self.datadir]
        with pytest.raises(SystemExit) as exc:
            main_monitor(args)
        assert exc.value.code == 2
        _, err = capfd.readouterr()
        assert f"Port {port} is already being used, choose a different one!" in err

    def test_cli_string_port(self, capfd):
        args = ["-p", "pippo", self.datadir]
        with pytest.raises(SystemExit) as exc:
            main_monitor(args)
        assert exc.value.code == 2
        _, err = capfd.readouterr()
        assert "Argument `port` should be an integer!" in err

    def test_cli_oor_workers(self, capfd):
        port = get_free_tcp_port()
        args = ["-p", f"{port}", "-w", "0", self.datadir]
        with pytest.raises(SystemExit) as exc:
            main_monitor(args)
        assert exc.value.code == 2
        _, err = capfd.readouterr()
        assert f"Choose a number of processes between 1 and {MAX_PROCS}" in err

    def test_cli_string_workers(self, capfd):
        port = get_free_tcp_port()
        args = ["-p", f"{port}", "-w", "pippo", self.datadir]
        with pytest.raises(SystemExit) as exc:
            main_monitor(args)
        assert exc.value.code == 2
        _, err = capfd.readouterr()
        assert f"Choose a number of processes between 1 and {MAX_PROCS}" in err

    def test_cli_nonexistent_config(self, capfd):
        port = get_free_tcp_port()
        config = "pippo"
        args = ["-p", f"{port}", "-c", f"{config}", self.datadir]
        with pytest.raises(SystemExit) as exc:
            main_monitor(args)
        assert exc.value.code == 2
        _, err = capfd.readouterr()
        assert f"Provided configuration file '{config}' does not exist!" in err


@pytest.mark.skipif("HAS_DEPENDENCIES")
class TestMonitorWithoutDeps:
    def test_dependencies_missing(self):
        with pytest.raises(ImportError) as exc:
            from srttools.monitor import Monitor
        assert np.any([f"install {dep}" in str(exc) for dep in ["tornado", "watchdog"]])

    def test_cli_dependencies_missing(self):
        with pytest.warns(UserWarning) as record:
            result = main_monitor(".")
        assert result == 1
        at_least_one_warning = False
        for string in ["watchdog", "tornado"]:
            at_least_one_warning = at_least_one_warning or np.any(
                [f"install {string}" in r.message.args[0] for r in record]
            )
        assert at_least_one_warning
