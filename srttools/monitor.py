from __future__ import (absolute_import, division,
                        print_function)
import time
import logging
import os
import shutil
import re
import sys
import signal
try:
    from watchdog.observers import Observer
    from watchdog.observers.polling import PollingObserver
    from watchdog.events import PatternMatchingEventHandler
    HAS_WATCHDOG = True
except ImportError:
    PatternMatchingEventHandler = object
    HAS_WATCHDOG = False

import warnings
import glob
import threading
from multiprocessing import Process, Queue, Manager, Lock, cpu_count
import warnings

try:
    from queue import Empty
except ImportError:
    from Queue import Empty
    warnings.warn("Monitoring interface does not work with Python < 3.5",
                  DeprecationWarning)

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


class MyEventHandler(PatternMatchingEventHandler):
    patterns = ["*/*.fits"]

    def __init__(self, n_proc, nosave=False):
        super().__init__()
        self.conf = getattr(sys.modules[__name__], 'conf')
        create_index_file(self.conf['debug_file_format'])
        self.filequeue = Queue()
        self.manager = Manager()
        self.lock = Lock()
        self.processing_queue = self.manager.list()
        proc_args = (
            self.filequeue,
            self.conf,
            nosave,
            self.lock,
            self.processing_queue
        )
        if n_proc == 1:
            p_type = threading.Thread
        else:
            p_type = Process
        self.processes = []
        for _ in range(n_proc):
            p = p_type(target=self._dequeue, args=proc_args)
            p.daemon = True
            p.start()
            self.processes.append(p)

    def __del__(self):
        for p in self.processes:
            try:
                p.terminate()
            except AttributeError:
                pass

    def _enqueue(self, infile):
        if infile not in self.processing_queue:
            self.filequeue.put(infile)
            self.processing_queue.append(infile)

    @staticmethod
    def _dequeue(filequeue, conf, nosave, lock, processing_queue):
        ext = conf['debug_file_format']
        while True:
            try:
                infile = filequeue.get()
            except (KeyboardInterrupt, EOFError):
                return

            productdir, fname = product_path_from_file_name(
                infile,
                productdir=conf['productdir'],
                workdir=conf['workdir']
            )
            root = os.path.join(productdir, fname.replace('.fits', ''))

            pp_args = ['--debug', '-c', conf['configuration_file_name']]
            if nosave:
                pp_args.append('--nosave')
            pp_args.append(infile)
            try:
                main_preprocess(pp_args)
            except:
                continue

            lock.acquire()

            newfiles = []
            for debugfile in glob.glob(root + '*.{}'.format(ext)):
                newfile = debugfile.replace(root, 'latest')
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

            oldfiles = glob.glob('latest*.{}'.format(ext))
            for oldfile in oldfiles:
                if oldfile not in newfiles and os.path.exists(oldfile):
                    os.remove(oldfile)

            lock.release()
            processing_queue.remove(infile)

    def on_modified(self, event):
        self._enqueue(event.src_path)

    def on_created(self, event):
        self._enqueue(event.src_path)

    def on_moved(self, event):
        self._enqueue(event.dest_path)


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
    import argparse

    description = ('Run the SRT quicklook in a given directory.')
    parser = argparse.ArgumentParser(description=description)

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
    args = parser.parse_args(args)

    def sigterm_received(signum, frame):
        os.kill(os.getpid(), signal.SIGINT)

    if threading.current_thread() is threading.main_thread():
        signal.signal(signal.SIGTERM, sigterm_received)

    if not HAS_WATCHDOG:
        raise ImportError('To use SDTmonitor, you need to install watchdog: \n'
                          '\n   > pip install watchdog')
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    if args.config is None:
        config_file = create_dummy_config()
    else:
        config_file = args.config
    conf = read_config(config_file)
    conf['configuration_file_name'] = config_file
    setattr(sys.modules[__name__], 'conf', conf)

    n_proc = 1
    if cpu_count() and not args.test:
        n_proc = int(min(cpu_count() / 2, 5))

    event_handler = MyEventHandler(n_proc=n_proc, nosave=args.nosave)
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
        count = 0
        while count < 10:
            time.sleep(1)
            if args.test:
                count += 1
        raise KeyboardInterrupt
    except KeyboardInterrupt:
        pass

    if args.http_server_port:
        http_server.shutdown()

    observer.stop()


def create_dummy_config():
    config_str = """
    [local]
    [analysis]
    [debugging]
    debug_file_format : jpg
    """
    with open('monitor_config.ini', 'w') as fobj:
        print(config_str, file=fobj)

    return 'monitor_config.ini'


def create_index_file(extension, max_images=50, interval=500):
    """
    :param extension: the file extension of the image files to look for.
    :param max_images: the maximum number of images to monitor. It should be
        enough to set it to twice the number of feeds of the receiver having
        the highest number of feeds (twice because of L and R channels).
        Its default value is set to 50, a number high enough to account for all
        the file images but not too high to represent a computational
        overhead. Since javascript cannot perform any client filesystem
        operation other than loading a local file from its path, the script
        below tries to access every image and just hides the not found ones,
        displaying only the images coming out of the monitor processing phase.
    :param interval: expressed in milliseconds, it represents the time between
        two subsequent calls to the `updatePage` function. Since the images
        get reloaded without any flickering (as opposed to when the whole page
        gets reloaded from the browser), its value can be set to a fraction of
        a second without any visible issue.
    """
    html_string = \
    """<html>
    <script type="text/javascript">
        window.onload = function()
        {
            var extension = '""" + extension + """';
            var maxImages = """ + str(max_images) + """;
			var interval = """ + str(interval) + """;

            document.body.innerHTML = "";

            for(i = 0; i < maxImages; i++)
            {
                document.body.innerHTML += '<div id="div_' + i + '" style="width:50%; float:left;"/></div>';
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

                        for(i = index; i < maxImages; i++)
                        {
                            var div = document.getElementById("div_" + i.toString());
                            div.innerHTML = "";
                        }
                    });
                }

                image.src = "latest_" + index.toString() + "." + extension + "?" + new Date().getTime();
            }

            update(0);
            setInterval(update, interval, 0);
        }
    </script>/n</html>"""

    with open('index.html', 'w') as fobj:
        print(html_string, file=fobj)
