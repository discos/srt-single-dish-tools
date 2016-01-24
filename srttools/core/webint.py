from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
from .read_config import read_config, get_config_file
from .scan import list_scans, Scan
import sqlite3


def load_scans(scan_list, norefilt=True):
    '''Load the scans in the list one by ones'''
    for f in scan_list:
        yield Scan(f, norefilt=norefilt)


def create_db(dbname='example.db'):
    config = read_config(get_config_file())
    conn = sqlite3.connect(dbname)

    c = conn.cursor()
    scan_list = list_scans(config['datadir'],
                           config['list_of_directories'])

    for s in load_scans(scan_list):
#        for ch in
        x = s['time']
        y = s['time']

