from srttools.inspect_observations import split_observation_table
from srttools.inspect_observations import dump_config_files
from srttools.inspect_observations import main_inspector
from astropy.table import Table, Column
import numpy as np
import os
import glob
import pytest
import logging
import subprocess as sp

try:
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser


@pytest.fixture()
def logger():
    logger = logging.getLogger('Some.Logger')
    logger.setLevel(logging.WARNING)

    return logger


class TestInspect(object):
    @classmethod
    def setup_class(klass):

        klass.curdir = os.path.dirname(__file__)

        info = Table()
        names = ["Dir", "Sample File", "Source", "Receiver", "Backend",
                 "Time", "Frequency", "Bandwidth", 'is_skydip']

        dtype = ['S200', 'S200', 'S200', 'S200', 'S200',
                 np.double, np.float, np.float, bool]

        for n, d in zip(names, dtype):
            if n not in info.keys():
                info.add_column(Column(name=n, dtype=d))

        times = np.arange(0, 1.3, 0.1)
        files = ["f"] * 13
        receivers = ["CCB"] * 3 + ["KKG"] * 4 + ["CCB"] * 3 + ["KKG"] * 3
        backends = ["SARDARA"] * 13
        frequency = [7000] * 3 + [20000] * 4 + [7000] * 3 + [20000] * 3
        bandwidth = [1000] * 3 + [500] * 4 + [1000] * 3 + [500] * 3
        sources = ("W44,3C48,3C295,3C157,3C157,3C48,3C48,"
                   "W44,3C48,3C48,W44,3C295,3C48").split(",")
        for i, t in enumerate(times):
            info.add_row(["{}_{:.1f}_{}".format(sources[i], times[i],
                                                receivers[i]),
                          files[i], sources[i], receivers[i],
                          backends[i], times[i], frequency[i], bandwidth[i],
                          False])
        klass.info = info
        klass.groups = split_observation_table(info)
        dump_config_files(klass.info)

    def test_0_files_exist(self):
        assert os.path.exists("CCB_SARDARA_W44_Obs0.ini")
        assert os.path.exists("KKG_SARDARA_3C157_Obs0.ini")
        assert os.path.exists("KKG_SARDARA_W44_Obs0.ini")
        assert os.path.exists("CCB_SARDARA_W44_Obs1.ini")

    def test_inspect_observations01(self):
        assert self.groups["CCB,SARDARA"]["W44"]["Obs0"]["Src"] == \
               ["W44_0.0_CCB"]

    def test_inspect_observations02(self):
        assert self.groups["CCB,SARDARA"]["W44"]["Obs0"]["Cal"] == \
               ["3C48_0.1_CCB", "3C295_0.2_CCB"]

    def test_inspect_observations03(self):
        assert self.groups["KKG,SARDARA"]["3C157"]["Obs0"]["Src"] == \
               ["3C157_0.3_KKG", "3C157_0.4_KKG"]

    def test_inspect_observations04(self):
        assert self.groups["KKG,SARDARA"]["3C157"]["Obs0"]["Cal"] == \
               ["3C48_0.5_KKG", "3C48_0.6_KKG"]

    def test_inspect_observations05(self):
        assert self.groups["CCB,SARDARA"]["W44"]["Obs1"]["Src"] == \
               ["W44_0.7_CCB"]

    def test_inspect_observations06(self):
        assert self.groups["CCB,SARDARA"]["W44"]["Obs1"]["Cal"] == \
               ["3C48_0.8_CCB", "3C48_0.9_CCB"]

    def test_inspect_observations07(self):
        assert self.groups["KKG,SARDARA"]["W44"]["Obs0"]["Src"] == \
               ["W44_1.0_KKG"]

    def test_inspect_observations08(self):
        assert self.groups["KKG,SARDARA"]["W44"]["Obs0"]["Cal"] == \
               ["3C48_0.6_KKG", "3C295_1.1_KKG", "3C48_1.2_KKG"]

    def test_inspect_observations09(self):
        config = ConfigParser()
        config.read("CCB_SARDARA_W44_Obs0.ini")
        entry = config.get("analysis", "calibrator_directories")
        assert entry.strip().split("\n") == \
            self.groups["CCB,SARDARA"]["W44"]["Obs0"]["Cal"]

    def test_inspect_observations10(self):
        config = ConfigParser()
        config.read("KKG_SARDARA_3C157_Obs0.ini")
        entry = config.get("analysis", "list_of_directories")
        assert entry.strip().split("\n") == \
            self.groups["KKG,SARDARA"]["3C157"]["Obs0"]["Src"]

    @classmethod
    def teardown_class(cls):
        """Cleanup."""
        os.unlink("KKG_SARDARA_3C157_Obs0.ini")
        os.unlink("KKG_SARDARA_W44_Obs0.ini")
        os.unlink("CCB_SARDARA_W44_Obs0.ini")
        os.unlink("CCB_SARDARA_W44_Obs1.ini")


class TestRun(object):
    @classmethod
    def setup_class(klass):

        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, 'data')

    def test_script_is_installed(self):
        sp.check_call('SDTinspect -h'.split(' '))

    def test_run_dump_default(self):
        main_inspector(
            glob.glob(os.path.join(self.datadir, 'gauss_*/')) + ['-d'])
        assert os.path.exists('CCB_TP_Dummy_Obs0.ini')
        assert os.path.exists('KKG_TP_Dummy_Obs0.ini')

    def test_run_dump(self):
        main_inspector(glob.glob(os.path.join(self.datadir, 'gauss_*/')) +
                       ['-g', 'Backend', '-d'])
        assert os.path.exists('TP_Dummy_Obs0.ini')

    def test_run_nodump(self, capsys):
        main_inspector(glob.glob(os.path.join(self.datadir, 'gauss_*/')))
        out, err = capsys.readouterr()
        assert 'Dummy' in out
        assert 'gauss_dec' in out
        assert 'Skydip' in out
        assert 'gauss_skydip' in out

    def test_run_date_filter_after(self, logger, caplog):
        main_inspector(glob.glob(os.path.join(self.datadir, 'gauss_*/')) +
                       '--only-after 20000101-000000'.split(' '))
        assert 'Filter out observations before MJD 51544' in caplog.text

    def test_run_date_filter_before(self, logger, caplog):
        main_inspector(glob.glob(os.path.join(self.datadir, 'gauss_*/')) +
                       '--only-before 21000101-000000'.split(' '))
        assert 'Filter out observations after MJD 88069' in caplog.text

    @classmethod
    def teardown_class(cls):
        """Cleanup."""
        import os
        os.unlink('TP_Dummy_Obs0.ini')
        os.unlink('CCB_TP_Dummy_Obs0.ini')
        os.unlink('KKG_TP_Dummy_Obs0.ini')
        os.unlink('table.csv')
        os.unlink('sample_config_file.ini')
