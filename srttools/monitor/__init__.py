import time
import os
import warnings
import signal
import argparse
import threading

from srttools.monitor.common import MAX_PROCS
from srttools.monitor.monitor import Monitor


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
    except OSError as e:  # This happens when the given port is already busy
        parser.error(str(e))
    except KeyboardInterrupt:
        pass
    if monitor is not None:
        monitor.stop()
