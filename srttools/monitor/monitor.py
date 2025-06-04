import glob
import os
import queue
import shutil
import signal
import sys
import threading
import time
import warnings
import re
import multiprocessing as mp
from srttools.imager import main_preprocess
from srttools.monitor.common import exit_function, log
from srttools.read_config import read_config
from srttools.scan import product_path_from_file_name

try:
    from watchdog.events import FileMovedEvent, RegexMatchingEventHandler
    from watchdog.observers import Observer
    from watchdog.observers.polling import PollingObserverVFS
except ImportError:
    raise ImportError(
        "To use SDTmonitor, you need to install watchdog: \n" "\n   > pip install watchdog"
    )

from srttools.monitor.webserver import WebServer

# Set the matplotlib backend
try:
    import matplotlib.pyplot as plt

    plt.switch_backend("Agg")
except ImportError:
    pass


def create_dummy_config(filename="monitor_config.ini", extension="png"):
    productdir = os.environ.get("TEMP", "") if "win" in sys.platform else "/tmp"
    config_str = f"""[local]\nproductdir : {productdir}\n[analysis]\n[debugging]\ndebug_file_format : {extension}"""
    with open(filename, "w") as fobj:
        print(config_str, file=fobj)
    return filename


class MyEventHandler(RegexMatchingEventHandler):
    ignore_regexes = [
        re.compile(r"^.*[/\\]tmp[/\\].*"),
        re.compile(r"^.*[/\\]tempfits[/\\].*"),
        re.compile(r"^.*[/\\][^/\\]+\.fitstemp$"),
        re.compile(r"^.*[/\\]summary[/\\].fits$"),
        re.compile(r"^.*[/\\]Sum_.*\.fits$"),
    ]
    regexes = [re.compile(r"^.*[/\\][^/\\]+\.fits(\d+)?$")]

    def __init__(self, observer):
        self._observer = observer
        self.on_modified = self._parse_filename
        self.on_created = self._parse_filename
        self.on_moved = self._parse_filename
        super().__init__()

    def _parse_filename(self, event):
        infile = ""
        if isinstance(event, FileMovedEvent):
            infile = event.dest_path
            if not any(re.match(regex, infile) for regex in self.regexes):
                return
        else:
            infile = event.src_path

        if self._observer._timers.get(infile):
            if not self._observer._timers[infile].processing:
                self._observer._timers[infile].cancel()
                del self._observer._timers[infile]
            else:
                return

        t = threading.Timer(1, self._observer._enqueue, args=(infile,))
        t.processing = False
        self._observer._timers[infile] = t
        self._observer._timers[infile].start()


class ProcessManager:
    def __init__(self, maxsize):
        self._semaphore = threading.Semaphore(maxsize)
        self._lock = threading.Lock()
        self._processes = set()
        self._terminated = threading.Event()

    def put(self, process, timeout=None):
        if not self._semaphore.acquire(timeout=timeout) or self._terminated.is_set():
            raise queue.Full
        with self._lock:
            self._processes.add(process)

    def remove(self, process):
        with self._lock:
            if not self._terminated.is_set() and process in self._processes:
                self._processes.remove(process)
                self._semaphore.release()

    def terminate(self):
        with self._lock:
            self._terminated.set()
            for process in set(self._processes):
                if process.is_alive():
                    process.terminate()
                self._processes.remove(process)
                self._semaphore.release()


class Monitor:
    def __init__(
        self,
        directories,
        config_file=None,
        workers=1,
        verbosity=0,
        polling=False,
        localhost=False,
        port=8080,
    ):
        # Save constructor parameters
        self._directories = directories
        self._config_file = config_file or create_dummy_config()
        self._workers = workers
        self._verbosity = verbosity
        self._polling = polling
        self._localhost = localhost
        self._port = port

        # Load file configuration and save needed parameters
        with warnings.catch_warnings():
            if not self._verbosity:
                warnings.simplefilter("ignore")
            configuration = read_config(self._config_file)
        self._productdir = configuration["productdir"]
        self._workdir = configuration["workdir"]
        self._basename_only = configuration["basename_only"]
        self._extension = configuration["debug_file_format"]

        # Objects needed by the file observer
        self._timers = {}
        self._processing = None
        self._files = None

        # Initialize the event handler and the file observer
        self._event_handler = MyEventHandler(self)
        if self._polling:
            self._observer = PollingObserverVFS(stat=os.stat, listdir=os.scandir)
        else:
            self._observer = Observer()

        for path in directories:
            self._observer.schedule(self._event_handler, path, recursive=True)

        self._stop = False
        self._worker_thread = None

        # Initialize the web server, this will raise a OSError if the given
        # port is already being used
        self._web_server = WebServer(self._extension, self._localhost, self._port)
        self._started = False

    def start(self):
        if self._started:
            log.info(f"SDTmonitor already running, process id: {self._pid}")
            return
        self._stop = False
        self._processing = ProcessManager(self._workers)
        self._files = queue.Queue()
        self._worker_thread = threading.Thread(target=self._worker_method)
        self._worker_thread.start()
        self._observer.start()
        self._web_server.start()
        self._pid = os.getpid()
        log.info(f"SDTmonitor started, process id: {self._pid}")
        self._started = True

    def stop(self):
        if not self._started:
            log.info("SDTmonitor is not running")
            return
        # Stop the worker thread
        self._stop = True
        if self._worker_thread:
            self._worker_thread.join()
        # Stop the web server
        self._web_server.stop()
        # Stop the observer from enqueuing newly arrived files
        self._observer.stop()
        # Terminate any running process
        self._processing.terminate()
        self._files.queue.clear()
        # Cancel and delete any pending timer
        for infile, timer in self._timers.copy().items():
            timer.cancel()
            timer.join()
            self._timers.pop(infile, None)
        self._started = False
        log.info(f"SDTmonitor stopped, process id: {self._pid}")

    def _worker_method(self):
        while not self._stop:
            try:
                to_update = []
                prefix, offset, oldfiles, prodpath = self._files.get(timeout=0.01)

                paths = {}
                originals = glob.glob(f"{prefix}_*.{self._extension}")
                for original in originals:
                    k, _ = os.path.splitext(original)
                    index = int(k.replace(f"{prefix}_", "")) + offset
                    paths[original] = f"latest_{index:03d}.{self._extension}"

                to_update = list(set(oldfiles) - set(paths.values()))
                for oldfile in to_update:
                    if os.path.exists(oldfile):
                        os.remove(oldfile)
                for oldfile, newfile in paths.items():
                    shutil.move(oldfile, newfile)
                    to_update.append(newfile)
                if prodpath:
                    empty_prodpath = not any(files for _, _, files in os.walk(prodpath))
                    shutil.rmtree(prodpath, ignore_errors=True) if empty_prodpath else None
                for image in to_update:
                    self._web_server.update(image)
            except queue.Empty:
                pass

    @staticmethod
    def _process(pp_args, verbosity):
        """Calls the main_preprocess function as a separate process, so that
        multiple processes can run concurrently speeding up the whole operation
        when receiving separate feeds files."""
        exit_code = 0
        try:
            with warnings.catch_warnings():
                if not verbosity:
                    warnings.simplefilter("ignore")
                exit_code = min(main_preprocess(pp_args), 1)
        except KeyboardInterrupt:
            exit_code = 15
        except Exception:  # pragma: no cover
            # We don't cover this scope since it's unexpected and not easily reproducible
            log.exception(sys.exc_info()[1])
            exit_code = 1
        exit_function(exit_code)

    def _enqueue(self, infile):
        self._timers[infile].processing = True
        proc_args = (
            ["--plot", "--nosave", "-c", self._config_file, "--pedantic", infile],
            self._verbosity,
        )
        p = mp.Process(target=self._process, args=proc_args)
        # The next call will stop if the queue is already full
        while not self._stop:
            try:
                self._processing.put(p, timeout=0.01)
                break
            except queue.Full:  # The queue is full, just sleep and wait for a free spot
                pass
        if self._stop:
            return
        p.start()
        log.info(f"Loading file {infile}, pid {p.pid}")

        # While the process executes, we retrieve
        # information regarding original and new files
        productdir, fname = product_path_from_file_name(
            infile,
            productdir=self._productdir,
            workdir=self._workdir,
            basename_only=self._basename_only,
        )
        root = os.path.join(productdir, fname.rsplit(".fits")[0])

        feed_idx = ""
        offset = 0
        if not infile.endswith(".fits"):
            feed_idx = infile.rsplit(".fits")[-1]
            offset = 2 * int(feed_idx)

        prefix = f"{root}{feed_idx}"

        prodpath = None
        if self._productdir and self._workdir not in self._productdir:
            prodpath = os.path.relpath(root, self._productdir)
            prodpath = prodpath.split("/")[0]
            prodpath = os.path.join(self._productdir, prodpath)

        # Retrieve the list of image files already in the page directory
        # They will be overwritten when new images come out
        oldfiles = []
        if not feed_idx:
            oldfiles = glob.glob(f"latest_*.{self._extension}")

        p.join(60)  # Wait a minute for process completion
        if p.is_alive():  # pragma: no cover
            # Process timed out, we don't have a file to simulate this scenario
            try:
                os.kill(p.pid, signal.SIGKILL)
                p.join()
            except ProcessLookupError:
                pass
        elif p.exitcode == 0:  # Completed successfully
            self._files.put((prefix, offset, oldfiles, prodpath))
            log.info(f"Completed file {infile}, pid {p.pid}")
        elif p.exitcode == 15:  # pragma: no cover
            # Forcefully terminated by KeyboardInterrupt
            log.info(f"Forcefully terminated process {p.pid}, file {infile}")
        else:  # Aborted
            log.info(f"Aborted file {infile}, pid {p.pid}")

        # Eventually notify that the queue is not full anymore
        self._processing.remove(p)
        del self._timers[infile]
