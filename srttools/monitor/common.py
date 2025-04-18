import os
import sys
import time
import warnings
import signal
import argparse
import threading

# Replace the default logging configuration with a custom one
from srttools import logger as log
import logging

log.name = "SDTmonitor"
f = logging.Formatter("%(asctime)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
log.handlers[0].setFormatter(f)
log.setLevel(logging.INFO)

MAX_FEEDS = 19
# The worst case scenario is fits0 through fits18, even though I'm not sure this notation will be kept in the future.
# We set a maximum cap for the number of concurrent processes to the number of CPU cores present in the system minus one.
# If there are more CPU cores than the maximum number of feeds, we set the number of maximum concurrent processes equal to the maximum number of feeds.
MAX_PROCS = min(MAX_FEEDS, os.cpu_count() - 1)

exit_function = os._exit

if "PYTEST_CURRENT_TEST" in os.environ:
    exit_function = sys.exit


stop_event = threading.Event()


def main_monitor(args=None):
    description = "Run the SRT quicklook in a given directory."
    parser = argparse.ArgumentParser(description=description)

    min_proc = 1
    max_proc = MAX_PROCS

    def workers_count(w):
        try:
            w = int(w)
            if not (w < min_proc or w > max_proc):
                return w
            else:
                raise ValueError
        except (ValueError, TypeError):
            raise argparse.ArgumentTypeError(
                f"Choose a number of processes between {min_proc} and {max_proc}."
            )

    def config_file(filename):
        if not filename:
            return ""
        elif os.path.isfile(filename):
            return filename
        else:
            raise argparse.ArgumentTypeError(
                f"Provided configuration file '{filename}' does not exist!"
            )

    def port_available(port):
        try:
            port = int(port)
        except ValueError:
            raise argparse.ArgumentTypeError("Argument `port` should be an integer!")
        if port == 0:
            # A 0 means random open port, we want to avoid this scenario
            # If the chosen port is busy the Monitor will return the same error
            raise argparse.ArgumentTypeError(
                f"Port {port} is already being used, choose a different one!"
            )
        return port

    parser.add_argument(
        "directories",
        help="Directories to monitor",
        default=None,
        nargs="+",
        type=str,
    )
    parser.add_argument(
        "-c",
        "--config",
        help="Configuration file",
        default="",
        type=config_file,
    )
    parser.add_argument(
        "--polling",
        help="Use a platform-independent, polling watchdog",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-l",
        "--localhost",
        help="The webserver will only listen from local connections",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-p",
        "--port",
        help="The port on which the server will be listening",
        type=port_available,
        default=8080,
    )
    parser.add_argument(
        "-v",
        "--verbosity",
        action="count",
        default=0,
        help="Set the verbosity level",
    )
    parser.add_argument(
        "-w",
        "--workers",
        type=workers_count,
        default=1,
        help="The maximum number of worker processes to spawn",
    )
    args = parser.parse_args(args)

    # This block is required to translate a SIGTERM into a KeyboardInterrupt, in order to handle the process as a service
    def sigterm_received(signum, frame):  # pragma: no cover
        os.kill(os.getpid(), signal.SIGINT)

    if threading.current_thread() is threading.main_thread():
        signal.signal(signal.SIGTERM, sigterm_received)

    monitor = None
    try:
        from srttools.monitor.monitor import Monitor

        monitor = Monitor(
            args.directories,
            config_file=args.config,
            workers=args.workers,
            verbosity=args.verbosity,
            polling=args.polling,
            localhost=args.localhost,
            port=args.port,
        )
        monitor.start()
        stop_event.clear()

        while not stop_event.is_set():
            time.sleep(0.1)
    except OSError as exc:  # This happens when the given port is already busy
        parser.error(str(exc))
    except KeyboardInterrupt:
        pass
    except ImportError as exc:
        warnings.warn(str(exc))
        return 1
    if monitor is not None:
        monitor.stop()
    return 0
