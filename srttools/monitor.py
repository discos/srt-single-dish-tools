import time
import logging
try:
    from watchdog.observers import Observer
    from watchdog.events import PatternMatchingEventHandler
    HAS_WATCHDOG = True
except ImportError:
    PatternMatchingEventHandler = object
    HAS_WATCHDOG = False
import warnings
import subprocess as sp
import glob

from srttools.read_config import read_config
warnings.filterwarnings('ignore')
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
        root = infile.replace('.fits', '')
        conf = read_config(CONFIG_FILE)
        ext = conf['debug_file_format']
        try:
            sp.check_call(
                "SDTpreprocess --debug {} -c {}".format(infile,
                                                        CONFIG_FILE).split())
        except sp.CalledProcessError:
            return

        for debugfile in glob.glob(root + '*.{}'.format(ext)):
            newfile = debugfile.replace(root, 'latest')
            sp.check_call('cp {} {}'.format(debugfile, newfile).split())

        with open('index.html', "w") as fobj:
            print('<META HTTP-EQUIV="refresh" CONTENT="5">', file=fobj)
            allfiles = glob.glob('latest*.{}'.format(ext))
            N = len(allfiles)
            if N <= 2:
                width = "50%"
            else:
                width = "25%"
            for fname in sorted(allfiles):
                print("<div style=\"width:{}; float:left;\" />".format(width),
                      file=fobj)
                print("<img src=\"{}\" width=\"100%\"/>".format(fname),
                      file=fobj)
                print("</div>", file=fobj)

    def on_created(self, event):
        self.process(event)

    def on_modified(self, event):
        self.process(event)


def main_monitor(args=None):
    import argparse
    global CONFIG_FILE

    description = ('Run the SRT quicklook in a given directory.')
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("directory",
                        help="Directory to monitor",
                        default=None, type=str)
    parser.add_argument("-c", "--config",
                        help="Config file",
                        default=None, type=str)
    parser.add_argument("--test",
                        help="Only to be used in tests!",
                        action='store_true', default=False)
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

    event_handler = MyHandler()
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
