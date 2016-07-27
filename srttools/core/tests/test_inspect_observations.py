from ..inspect_observations import split_observation_table
from astropy.table import Table, Column
import numpy as np


def test_inspect_observations():
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
    receivers = ["CCB"] * 3 + ["KKG"] * 4 + ["CCB"] * 3 +  ["KKG"] * 3
    backends = ["SARDARA"] * 13
    frequency = [7000] * 3 + [20000] * 4 + [7000] * 3 + [20000] * 3
    bandwidth = [1000] * 3 + [500] * 4 + [1000] * 3 + [500] * 3
    sources = "W44,3C48,3C295,3C157,3C157,3C48,3C48,W44,3C48,3C48,W44,3C295,3C48".split(",")
    for i, t in enumerate(times):
        info.add_row(["{}_{:.1f}_{}".format(sources[i], times[i], receivers[i]),
                      files[i], sources[i], receivers[i],
                      backends[i], times[i], frequency[i], bandwidth[i]])

    groups = split_observation_table(info)

    assert groups["CCB,SARDARA"]["W44"]["Obs0"]["Src"] == ["W44_0.0_CCB"]
    assert groups["CCB,SARDARA"]["W44"]["Obs0"]["Cal"] == ["3C48_0.1_CCB", "3C295_0.2_CCB"]
    assert groups["KKG,SARDARA"]["3C157"]["Obs0"]["Src"] == ["3C157_0.3_KKG", "3C157_0.4_KKG"]
    assert groups["KKG,SARDARA"]["3C157"]["Obs0"]["Cal"] == ["3C48_0.5_KKG", "3C48_0.6_KKG"]
    assert groups["CCB,SARDARA"]["W44"]["Obs1"]["Src"] == ["W44_0.7_CCB"]
    assert groups["CCB,SARDARA"]["W44"]["Obs1"]["Cal"] == ["3C48_0.8_CCB", "3C48_0.9_CCB"]
    assert groups["KKG,SARDARA"]["W44"]["Obs0"]["Src"] == ["W44_1.0_KKG"]
    assert groups["KKG,SARDARA"]["W44"]["Obs0"]["Cal"] == ["3C48_0.6_KKG", "3C295_1.1_KKG", "3C48_1.2_KKG"]
