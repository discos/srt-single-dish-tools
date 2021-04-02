from astropy.table import Table


def test_import_scan():
    from srttools.scan import Scan

    s = Scan()
    assert isinstance(s, Table)


def test_import_scanset():
    from srttools.imager import ScanSet

    s = ScanSet()
    assert isinstance(s, Table)


def test_import_calibratortable():
    from srttools.calibration import CalibratorTable

    s = CalibratorTable()
    assert isinstance(s, Table)
