import time
import os
import warnings
import signal
import argparse
import threading

try:
    from srttools.monitor.monitor import Monitor
    from srttools.monitor.common import MAX_FEEDS
except Exception:
    warnings.warn("Monitor cannot work")


def main_monitor(args=None):
    description = "Run the SRT quicklook in a given directory."
    parser = argparse.ArgumentParser(description=description)

    min_proc = 1
    max_proc = MAX_FEEDS

    def workers_count(w):
        try:
            w = int(w)
            if not (w < min_proc or w > max_proc):
                return w
            else:
                raise ValueError
        except (ValueError, TypeError):
            raise argparse.ArgumentTypeError(
                "Choose a number of processes between {} and {}.".format(min_proc, max_proc)
            )

    def config_file(filename):
        if not filename:
            return ""
        elif os.path.isfile(filename):
            return filename
        else:
            raise argparse.ArgumentTypeError(
                "Provided configuration file '{}' does not exist!".format(filename)
            )

    def port_available(port):
        try:
            port = int(port)
        except ValueError:
            raise argparse.ArgumentTypeError("Argument `port` should be an integer!")
        import socket

        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex(("localhost", port)) == 0:
                raise argparse.ArgumentTypeError(
                    "Port {} is already being used, choose a different one!".format(port)
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
        "-p",
        "--port",
        help="The port on which the server will be listening",
        type=int,
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
    def sigterm_received(signum, frame):
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
            port=args.port,
        )
        monitor.start()

        while True:
            time.sleep(0.1)
    except OSError as e:  # This happens when the given port is already busy
        parser.error(str(e))
    except KeyboardInterrupt:
        pass
    if monitor is not None:
        monitor.stop()
