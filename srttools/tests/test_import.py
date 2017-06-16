from srttools import Scan, ScanSet, CalibratorTable
from astropy.table import Table


def test_import_scan():
    s = Scan()
    assert isinstance(s, Table)


def test_import_scanset():
    s = ScanSet()
    assert isinstance(s, Table)


def test_import_calibratortable():
    s = CalibratorTable()
    assert isinstance(s, Table)
