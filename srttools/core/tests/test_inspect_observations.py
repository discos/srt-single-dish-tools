from ..inspect_observations import split_observation_table, dump_config_files
from astropy.table import Table, Column
import numpy as np

try:
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser


class TestInspect(object):
    @classmethod
    def setup_class(klass):
        import os
        global DEBUG_MODE
        DEBUG_MODE = True

        klass.curdir = os.path.dirname(__file__)

        info = Table()
        names = ["Dir", "Sample File", "Source", "Receiver", "Backend",
                 "Time", "Frequency", "Bandwidth"]

        dtype = ['S200', 'S200', 'S200', 'S200', 'S200',
                 np.double, np.float, np.float]

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
                          backends[i], times[i], frequency[i], bandwidth[i]])
        klass.info = info
        klass.groups = split_observation_table(info)

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
        dump_config_files(self.info)
        config = ConfigParser()
        config.read("CCB_SARDARA_W44_Obs0.ini")
        entry = config.get("analysis", "calibrator_directories")
        assert entry.strip().split("\n") == \
            self.groups["CCB,SARDARA"]["W44"]["Obs0"]["Cal"]

    def test_inspect_observations10(self):
        dump_config_files(self.info)
        config = ConfigParser()
        config.read("KKG_SARDARA_3C157_Obs0.ini")
        entry = config.get("analysis", "list_of_directories")
        assert entry.strip().split("\n") == \
            self.groups["KKG,SARDARA"]["3C157"]["Obs0"]["Src"]

    @classmethod
    def teardown_class(cls):
        """Cleanup."""
        import os
        import glob
        for ini in glob.glob('*.ini'):
            os.unlink(ini)
