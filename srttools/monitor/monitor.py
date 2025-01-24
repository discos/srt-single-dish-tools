import glob
import multiprocessing as mp
import os
import queue
import shutil
import signal
import sys
import threading
import time
import warnings

from srttools.imager import main_preprocess
from srttools.monitor.common import MAX_FEEDS, exit_function, log
from srttools.read_config import read_config
from srttools.scan import product_path_from_file_name

try:
    from watchdog.events import FileMovedEvent, PatternMatchingEventHandler
    from watchdog.observers import Observer
    from watchdog.observers.polling import PollingObserver

    from srttools.monitor.webserver import WebServer
except ImportError:
    warnings.warn(
        "To use SDTmonitor, you need to install watchdog: \n" "\n   > pip install watchdog"
    )
    PatternMatchingEventHandler = object

# Set the matplotlib backend
try:
    import matplotlib.pyplot as plt

    plt.switch_backend("Agg")
except ImportError:
    pass


def create_dummy_config(filename="monitor_config.ini", extension="png"):
    config_str = f"""[local]\n[analysis]\n[debugging]\ndebug_file_format : {extension}"""
    with open(filename, "w") as fobj:
        print(config_str, file=fobj)
    return filename


class MyEventHandler(PatternMatchingEventHandler):
    ignore_patterns = [
        "*/tmp/*",
        "*/tempfits/*",
        "*/*.fitstemp",
        "*/summary.fits",
        "*/Sum_*.fits",
    ]
    patterns = ["*/*.fits"] + [f"*/*.fits{x}" for x in range(MAX_FEEDS)]

    def __init__(self, observer):
        self._observer = observer
        self.on_modified = self._parse_filename
        self.on_created = self._parse_filename
        self.on_moved = self._parse_filename
        super().__init__()

    def _parse_filename(self, event):
        infile = ""
        if isinstance(event, FileMovedEvent):
            for pattern in self.patterns:
                pattern = pattern.rsplit(".")[-1]
                if event.dest_path.endswith(pattern):
                    infile = event.dest_path
            if not infile:
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


class Monitor:
    def __init__(
        self,
        directories,
        config_file="",
        workers=1,
        verbosity=0,
        polling=False,
        port=8080,
    ):
        # Save constructor parameters
        self._directories = directories
        if not config_file:
            config_file = create_dummy_config()
        self._config_file = config_file
        self._workers = workers
        self._verbosity = verbosity
        self._polling = polling
        self._port = port

        # Load file configuration and save needed parameters
        with warnings.catch_warnings():
            if not self._verbosity:
                warnings.simplefilter("ignore")
            configuration = read_config(config_file)
        self._productdir = configuration["productdir"]
        self._workdir = configuration["workdir"]
        self._extension = configuration["debug_file_format"]

        # Objects needed by the file observer
        self._timers = {}
        self._processing = queue.Queue(self._workers)
        self._files = queue.Queue()

        # Initialize the event handler and the file observer
        self._event_handler = MyEventHandler(self)
        if self._polling:
            self._observer = PollingObserver()
        else:
            self._observer = Observer()

        for path in directories:
            self._observer.schedule(self._event_handler, path, recursive=True)

        self._stop = False
        self._worker_thread = None

        # Initialize the web server, this will raise a OSError if the given
        # port is already being used
        self._web_server = WebServer(self._extension, self._port)

    def start(self):
        self._stop = False
        self._worker_thread = threading.Thread(target=self._worker_method)
        self._worker_thread.start()
        self._observer.start()
        self._web_server.start()
        log.info(f"SDTmonitor started, process id: {os.getpid()}")

    def stop(self):
        # Stop the worker thread
        self._stop = True
        if self._worker_thread:
            self._worker_thread.join()
        # Stop the web server
        self._web_server.stop()
        # Stop the observer from enqueuing newly arrived files
        self._observer.stop()
        # Cancel and delete any pending timer
        for infile, timer in self._timers.copy().items():
            timer.cancel()
            timer.join()
            try:
                del self._timers[infile]
            except KeyError:
                pass
        # Terminate any running process
        while True:
            try:
                p = self._processing.get_nowait()
                p.terminate()
            except AttributeError:  # Process enqueued but not started yet
                pass
            except queue.Empty:
                break
        log.info(f"SDTmonitor stopped, process id: {os.getpid()}")

    def _worker_method(self):
        while not self._stop:
            try:
                to_update = []
                paths, oldfiles, prodpath = self._files.get_nowait()

                for key in [key for key in paths if not os.path.exists(key)]:
                    del paths[key]
                for oldfile in oldfiles:
                    if oldfile not in paths.values() and os.path.exists(oldfile):
                        os.remove(oldfile)
                        to_update.append(oldfile)
                for oldfile, newfile in paths.items():
                    shutil.move(oldfile, newfile)
                    to_update.append(newfile)
                if prodpath:
                    for dirname, _, _ in os.walk(prodpath, topdown=False):
                        try:
                            if not os.listdir(dirname):
                                os.rmdir(dirname)
                        except OSError:
                            pass
                for image in to_update:
                    self._web_server.update(image)
            except queue.Empty:
                pass
            time.sleep(0.01)

    @staticmethod
    def _process(pp_args, verbosity):
        """Calls the main_preprocess function as a separate process, so that
        multiple processes can run concurrently speeding up the whole operation
        when receiving separate feeds files.
        """
        exit_code = 0
        try:
            with warnings.catch_warnings():
                if not verbosity:
                    warnings.simplefilter("ignore")
                if main_preprocess(pp_args):
                    exit_code = 1
        except KeyboardInterrupt:
            exit_code = 15
        except Exception:
            log.exception(sys.exc_info()[1])
            exit_code = 1
        exit_function(exit_code)

    def _enqueue(self, infile):
        self._timers[infile].processing = True
        proc_args = (
            ["--plot", "--nosave", "-c", self._config_file, infile],
            self._verbosity,
        )
        p = mp.Process(target=self._process, args=proc_args)
        # The next call will stop if the queue is already full
        while not self._stop:
            try:
                self._processing.put_nowait(p)
                break
            except queue.Full:  # The queue is full, just sleep and wait for a free spot
                pass
            time.sleep(0.01)
        if self._stop:
            return
        p.start()
        log.info(f"Loading file {infile}, pid {p.pid}")

        # While the process executes, we retrieve
        # information regarding original and new files
        productdir, fname = product_path_from_file_name(
            infile, productdir=self._productdir, workdir=self._workdir
        )
        root = os.path.join(productdir, fname.rsplit(".fits")[0])

        feed_idx = ""
        offset = 0
        if not infile.endswith(".fits"):
            feed_idx = infile.rsplit(".fits")[-1]
            offset = 2 * int(feed_idx)

        paths = {
            f"{root}{feed_idx}_{i}.{self._extension}": f"latest_{i + offset}.{self._extension}"
            for i in range(2 if feed_idx else MAX_FEEDS * 2)
        }

        prodpath = None
        if self._productdir and self._workdir not in self._productdir:
            prodpath = os.path.relpath(root, self._productdir)
            prodpath = prodpath.split("/")[0]
            prodpath = os.path.join(self._productdir, prodpath)

        # Retrieve the list of image files already in the page directory
        # They will be overwritten when new images come out
        oldfiles = []
        if not feed_idx:
            oldfiles = glob.glob(f"latest*.{self._extension}")

        p.join(60)  # Wait a minute for process completion
        if p.is_alive():  # Process timed out
            try:
                os.kill(p.pid, signal.SIGKILL)
                p.join()
            except ProcessLookupError:
                pass
        elif p.exitcode == 0:  # Completed successfully
            self._files.put((paths, oldfiles, prodpath))
            log.info(f"Completed file {infile}, pid {p.pid}")
        elif p.exitcode == 1:  # Aborted
            log.info(f"Aborted file {infile}, pid {p.pid}")
        elif p.exitcode == 15:  # Forcefully terminated
            log.info(f"Forcefully terminated process {p.pid}, file {infile}")
        else:  # Unexpected code
            log.info(f"Process {p.pid} exited with unexpected code {p.exitcode}")

        # Eventually notify that the queue is not full anymore
        try:
            with self._processing.not_empty:
                self._processing.queue.remove(p)
                self._processing.not_full.notify()
        except ValueError:
            pass
        del self._timers[infile]
