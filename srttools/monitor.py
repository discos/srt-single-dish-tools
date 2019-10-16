from __future__ import (absolute_import, division,
                        print_function)
import time
import os
import shutil
import re
import sys
import signal
import argparse
try:
    from watchdog.observers import Observer
    from watchdog.observers.polling import PollingObserver
    from watchdog.events import PatternMatchingEventHandler, FileMovedEvent
    HAS_WATCHDOG = True
except ImportError:
    PatternMatchingEventHandler = object
    HAS_WATCHDOG = False

import warnings
import glob
import threading
from multiprocessing import Process, Lock

from astropy import log

try:
    from http.server import HTTPServer, SimpleHTTPRequestHandler, HTTPStatus
except ImportError:
    from BaseHTTPServer import HTTPServer
    from SimpleHTTPServer import SimpleHTTPRequestHandler
    class HTTPStatus:
        pass
    HTTPStatus.NOT_FOUND = 404

from srttools.read_config import read_config
from srttools.scan import product_path_from_file_name
from srttools.imager import main_preprocess
try:
    import matplotlib.pyplot as plt
    plt.switch_backend('Agg')
except ImportError:
    pass

MAX_FEEDS = 7

class MyEventHandler(PatternMatchingEventHandler):
    patterns = \
        ["*/*.fits"] + ["*/*.fits{}".format(x) for x in range(MAX_FEEDS)]

    def __init__(self, n_proc, nosave=False, verbosity=0, test=False):
        super().__init__()
        self.conf = getattr(sys.modules[__name__], 'conf')
        create_index_file(self.conf['debug_file_format'])

        self.max_proc = n_proc

        self.timers = {}
        self.lock = Lock()
        self.processing_queue = []
        self.waiting_queue = []

        self.nosave = nosave
        self.verbosity = verbosity

        self.on_modified = self._start_timer
        self.on_created = self._start_timer
        self.on_moved = self._start_timer

        self.test = test

    def _start_timer(self, event):
        infile = ''
        if isinstance(event, FileMovedEvent):
            for pattern in self.patterns:
                pattern = pattern.rsplit('.')[-1]
                if event.dest_path.endswith(pattern):
                    infile = event.dest_path
            if not infile:
                return
        else:
            infile = event.src_path

        if self.timers.get(infile):
            self.timers[infile].cancel()

        self.timers[infile] = threading.Timer(
            1,
            self._enqueue,
            args=[infile]
        )
        self.timers[infile].daemon = True
        self.timers[infile].start()

    def _enqueue(self, infile):
        if self.timers.get(infile):
            del self.timers[infile]
        if infile not in self.waiting_queue:
            self.waiting_queue.append(infile)
            while len(self.processing_queue) == self.max_proc:
                time.sleep(0.05)
            self.waiting_queue.remove(infile)
            self.processing_queue.append(infile)
            proc_args = (
                infile,
                self.conf,
                self.nosave,
                self.verbosity,
                self.lock
            )
            if self.test:
                p = threading.Thread(target=self.process, args=proc_args)
            else:
                p = Process(target=self.process, args=proc_args)
            p.daemon = True
            p.start()
            p.join()
            try:
                if p.exitcode not in [0, 1]:
                    log.info('Process {} exited with unexpected code {}'.format(
                        p.pid, p.exitcode
                    ))
            except AttributeError:
                pass
            self.processing_queue.remove(infile)

    @staticmethod
    def process(infile, conf, nosave, verbosity, lock):
        try:
            pid = os.getpid()
            log.info('Loading file {}, pid {}'.format(infile, pid))
            ext = conf['debug_file_format']

            productdir, fname = product_path_from_file_name(
                infile,
                productdir=conf['productdir'],
                workdir=conf['workdir']
            )
            root = os.path.join(productdir, fname.rsplit('.fits')[0])

            pp_args = ['--plot', '-c', conf['configuration_file_name']]
            if nosave:
                pp_args.append('--nosave')
            pp_args.append(infile)

            skip = False
            try:
                with warnings.catch_warnings():
                    if not verbosity:
                        warnings.simplefilter('ignore')
                    if main_preprocess(pp_args):
                        skip = True
            except KeyboardInterrupt:
                raise KeyboardInterrupt
            except:
                log.exception(sys.exc_info()[1])
                skip = True

            if skip:
                log.info('Aborted file {}, pid {}'.format(infile, pid))
                sys.exit(1)

            feed_idx = ''
            if not infile.endswith('.fits'):
                feed_idx = infile.rsplit('.fits')[-1]

            lock.acquire()

            newfiles = []
            for debugfile in glob.glob(root + '{}_*.{}'.format(feed_idx, ext)):
                if feed_idx:
                    newfile = debugfile.replace(root + infile[-1], 'latest')
                else:
                    newfile = debugfile.replace(root, 'latest')
                newfile = newfile.rstrip('.{}'.format(ext))
                if feed_idx:
                    new_index = 2 * int(infile.rsplit('.fits')[-1])
                    newfile = (
                        newfile.split('_')[0]
                        + '_'
                        + str(new_index + int(newfile.rsplit('_')[-1]))
                    )
                newfile = newfile + '.{}'.format(ext)
                newfiles.append(newfile)
                if os.path.exists(debugfile):
                    if nosave:
                        shutil.move(debugfile, newfile)
                    else:
                        shutil.copyfile(debugfile, newfile)
            if nosave and conf['productdir'] \
                    and conf['workdir'] not in conf['productdir']:
                prodpath = os.path.relpath(root, conf['productdir'])
                prodpath = prodpath.split('/')[0]
                prodpath = os.path.join(conf['productdir'], prodpath)
                for dirname, _, _ in os.walk(prodpath, topdown=False):
                    if not os.listdir(dirname):
                        os.rmdir(dirname)

            if not feed_idx:
                oldfiles = glob.glob('latest*.{}'.format(ext))
                for oldfile in oldfiles:
                    if oldfile not in newfiles and os.path.exists(oldfile):
                        os.remove(oldfile)

            lock.release()

            log.info('Completed file {}, pid {}'.format(infile, pid))
        except KeyboardInterrupt:
            pass
        sys.exit(0)


class MyRequestHandler(SimpleHTTPRequestHandler):

    def __init__(self, request, client_address, server):
        self.allowed_paths = ['index.html', 'index.htm']
        conf = getattr(sys.modules[__name__], 'conf')
        self.re_pattern = '^latest_([0-9]*).%s$' % conf['debug_file_format']
        SimpleHTTPRequestHandler.__init__(self, request, client_address, server)

    def do_GET(self):
        path = self.path.lstrip('/')
        if path == '':
            path = 'index.html'
        elif '?' in path:
            path = path.split('?')[0]
        if path not in self.allowed_paths \
                and not re.match(self.re_pattern, path):
            self.send_error(HTTPStatus.NOT_FOUND)
            return
        try:
            SimpleHTTPRequestHandler.do_GET(self)
        except (BrokenPipeError, ConnectionResetError):
            pass

    def log_message(self, format, *args):
        return


def main_monitor(args=None):
    description = ('Run the SRT quicklook in a given directory.')
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
                "Choose an number of processes between {} and {}.".format(
                    min_proc, max_proc
                )
            )

    parser.add_argument(
        "directories",
        help="Directories to monitor",
        default=None,
        nargs='+',
        type=str
    )
    parser.add_argument(
        "-c", "--config",
        help="Config file",
        default=None,
        type=str
    )
    parser.add_argument(
        "--test",
        help="Only to be used in tests!",
        action='store_true',
        default=False
    )
    parser.add_argument(
        "--nosave",
        help="Do not save the hdf5 intermediate files",
        action='store_true',
        default=False
    )
    parser.add_argument(
        "-p", "--polling",
        help="Use a platform-independent, polling watchdog",
        action='store_true',
        default=False
    )
    parser.add_argument(
        "--http-server-port",
        help="Share the results via HTTP server on given HTTP_SERVER_PORT",
        type=int
    )
    parser.add_argument(
        '-v', '--verbosity',
        action='count',
        default=0,
        help='Set the verbosity level'
    )
    parser.add_argument(
        '-w', '--workers',
        type=workers_count,
        default=1,
        help='The maximum number of worker processes to spawn'
    )
    parser.add_argument(
        '--timeout',
        type=float,
        default=None,
        help='Set a timeout for the monitor'
    )
    args = parser.parse_args(args)

    timeout = args.timeout
    if args.timeout is None:
        timeout = 1e32
        if args.test:
            timeout = 10

    if not args.test:
        del log.handlers[:]  # Needed to prevent duplicate logging entries

    from astropy.logger import logging
    logging.basicConfig(

        level=logging.INFO,
        format='%(asctime)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    log.info('Starting monitor process, pid {}'.format(os.getpid()))

    def sigterm_received(signum, frame):
        os.kill(os.getpid(), signal.SIGINT)

    if threading.current_thread() is threading.main_thread():
        signal.signal(signal.SIGTERM, sigterm_received)

    if not HAS_WATCHDOG:
        raise ImportError('To use SDTmonitor, you need to install watchdog: \n'
                          '\n   > pip install watchdog')

    if args.config is None:
        config_file = create_dummy_config()
    else:
        config_file = args.config

    with warnings.catch_warnings():
        if not args.verbosity:
            warnings.simplefilter('ignore')
        conf = read_config(config_file)

    conf['configuration_file_name'] = config_file
    setattr(sys.modules[__name__], 'conf', conf)

    event_handler = MyEventHandler(
        n_proc=args.workers,
        nosave=args.nosave,
        test=args.test
    )
    observer = None
    if args.polling:
        observer = PollingObserver()
    else:
        observer = Observer()

    for path in args.directories:
        observer.schedule(event_handler, path, recursive=True)

    observer.start()

    if args.http_server_port:
        http_server = HTTPServer(('', args.http_server_port), MyRequestHandler)
        t = threading.Thread(target=http_server.serve_forever)
        t.daemon = True
        t.start()

    try:
        t0 = time.time()
        while time.time() - t0 < timeout:
            time.sleep(1)

        raise KeyboardInterrupt
    except KeyboardInterrupt:
        pass

    if args.http_server_port:
        http_server.shutdown()

    observer.stop()
    log.info('Stopping monitor process, pid {}'.format(os.getpid()))


def create_dummy_config():
    config_str = """
    [local]
    [analysis]
    [debugging]
    debug_file_format : png
    """
    with open('monitor_config.ini', 'w') as fobj:
        print(config_str, file=fobj)

    return 'monitor_config.ini'


def create_index_file(extension, interval=500):
    """
    :param extension: the file extension of the image files to look for.
    :param interval: expressed in milliseconds, it represents the time between
        two subsequent calls to the ``updatePage`` function. Since the images
        get reloaded without any flickering (as opposed to when the whole page
        gets reloaded from the browser), its value can be set to a fraction of
        a second without any visible issue.
    """
    html_string = \
    """<html>
    <script type="text/javascript">
        window.onload = function()
        {
            document.body.innerHTML = "";

            function remove_div(index)
            {
                var div = document.getElementById("div_" + index.toString());
                if(div != null)
                {
                    div.parentElement.removeChild(div);
                    remove_div(index + 1);
                }
            }

            function update(index)
            {
                image_id = "img_" + index.toString();

                var image = document.getElementById(image_id);

                if(image == null)
                {
                    image = new Image();
                    image.id = image_id;
                    image.style.width = "100%";

                    image.addEventListener("load", function()
                    {
                        if(this.parentElement == null)
                        {
                            index = parseInt(this.id.split("_")[1]);

                            var div = document.getElementById("div_" + index.toString());

                            if(div == null)
                            {
                                div = document.createElement("DIV");
                                div.setAttribute("id", "div_" + index.toString());
                                div.setAttribute("style", "width:50%; float:left;");
                                document.body.appendChild(div);
                            }

                            while(div.firstChild)
                            {
                                div.removeChild(div.firstChild);
                            }
                            div.appendChild(this);
                        }

                        update(index + 1);
                    });

                    image.addEventListener("error", function()
                    {
                        index = parseInt(this.id.split("_")[1]);
                        remove_div(index);
                    });
                }

                image.src = "latest_" + index.toString() + ".""" + extension + """?" + new Date().getTime();
            }

            update(0);
            setInterval(update, """ + str(interval) + """, 0);
        }
    </script>\n</html>"""

    with open('index.html', 'w') as fobj:
        print(html_string, file=fobj)
