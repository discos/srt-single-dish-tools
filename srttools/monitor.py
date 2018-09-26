from __future__ import (absolute_import, division,
                        print_function)
import time
import logging
import os
try:
    from watchdog.observers import Observer
    from watchdog.observers.polling import PollingObserver
    from watchdog.events import PatternMatchingEventHandler
    HAS_WATCHDOG = True
except ImportError:
    PatternMatchingEventHandler = object
    HAS_WATCHDOG = False
import warnings
import subprocess as sp
import glob

from srttools.read_config import read_config
from srttools.scan import product_path_from_file_name
global CONFIG_FILE


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


class MyHandler(PatternMatchingEventHandler):
    patterns = ["*/*.fits"]

    def __init__(self, nosave=False):
        self.nosave = nosave
        super().__init__()

    def process(self, event):
        """
        event.event_type
            'modified' | 'created' | 'moved' | 'deleted'
        event.is_directory
            True | False
        event.src_path
            path/to/observed/file
        """
        global CONFIG_FILE

        infile = event.src_path
        conf = read_config(CONFIG_FILE)
        ext = conf['debug_file_format']
        productdir, fname = \
            product_path_from_file_name(infile, productdir=conf['productdir'],
                                        workdir=conf['workdir'])

        root = os.path.join(productdir, fname.replace('.fits', ''))

        try:
            cmd_string = "SDTpreprocess --debug {} "
            if self.nosave:
                cmd_string += "--nosave "
            cmd_string += "-c {}"
            sp.check_call(cmd_string.format(infile, CONFIG_FILE).split())
        except sp.CalledProcessError:
            return

        newfiles = []
        for debugfile in glob.glob(root + '*.{}'.format(ext)):
            newfile = debugfile.replace(root, 'latest')
            newfiles.append(newfile)
            cmd_string = ''
            if self.nosave:
                cmd_string = 'mv {} {}'
            else:
                cmd_string = 'cp {} {}'
            sp.check_call(cmd_string.format(debugfile, newfile).split())

        oldfiles = glob.glob('latest*.{}'.format(ext))
        with open('index.html', "w") as fobj:
            print('<META HTTP-EQUIV="refresh" CONTENT="5">', file=fobj)
            N = len(newfiles)
            if N <= 2:
                width = "50%"
            else:
                width = "25%"
            for fname in sorted(newfiles):
                print("<div style=\"width:{}; float:left;\" />".format(width),
                      file=fobj)
                print("<img src=\"{}\" width=\"100%\"/>".format(fname),
                      file=fobj)
                print("</div>", file=fobj)
        for oldfile in oldfiles:
            if oldfile not in newfiles:
                os.remove(oldfile)

    def on_created(self, event):
        self.process(event)

    def on_modified(self, event):
        self.process(event)


def main_monitor(args=None):
    import argparse
    global CONFIG_FILE

    description = ('Run the SRT quicklook in a given directory.')
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "directory",
        help="Directory to monitor.",
        default=None,
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
        help="Do not save the hdf5 intermediate files.",
        action='store_true',
        default=False
    )
    parser.add_argument(
        "-p", "--polling",
        help="Use a platform-independent, polling watchdog.",
        action='store_true',
        default=False
    )
    args = parser.parse_args(args)

    if not HAS_WATCHDOG:
        raise ImportError('To use SDTmonitor, you need to install watchdog: \n'
                          '\n   > pip install watchdog')
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    with open('index.html', "w") as fobj:
        print('<META HTTP-EQUIV="refresh" CONTENT="5">', file=fobj)
        print("Waiting for the first observation...", file=fobj)

    path = args.directory

    if args.config is None:
        CONFIG_FILE = create_dummy_config()
    else:
        CONFIG_FILE = args.config

    event_handler = MyHandler(nosave=args.nosave)
    observer = None
    if args.polling:
        observer = PollingObserver()
    else:
        observer = Observer()
    observer.schedule(event_handler, path, recursive=True)
    observer.start()
    try:
        count = 0
        while count < 10:
            time.sleep(1)
            if args.test:
                count += 1
        raise KeyboardInterrupt
    except KeyboardInterrupt:
        pass
    observer.stop()
