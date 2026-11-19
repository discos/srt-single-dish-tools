import glob
import os
import shutil
import warnings

import pytest

from srttools.imager import main_preprocess
from srttools.rfistat import main_rfistat
from srttools.utils import HAS_MPL

DATADIR = os.path.join(os.path.dirname(__file__), "data")
RAW_FILES = sorted(glob.glob(os.path.join(DATADIR, "nodding_xarcos", "*.fits")))

CONFIG = """[local]
workdir = {tmp}
datadir = {tmp}/raw
productdir = {tmp}/products

[analysis]
list_of_directories =
    nodding_xarcos
"""


def _make_config(tmp, raw_files):
    """Copy only the given raw FITS files into ``tmp/raw`` and write a config file."""
    rawdir = tmp / "raw" / "nodding_xarcos"
    rawdir.mkdir(parents=True)
    for f in raw_files:
        shutil.copy(f, rawdir)
    config = tmp / "nodding.ini"
    config.write_text(CONFIG.format(tmp=tmp))
    return str(config)


@pytest.fixture(scope="module")
def preprocessed(tmp_path_factory):
    """The raw FITS files of nodding_xarcos in a temporary directory, preprocessed."""
    tmp = tmp_path_factory.mktemp("rfi")
    config = _make_config(tmp, RAW_FILES)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        main_preprocess(["-c", config])
    return tmp, config


@pytest.mark.slow
@pytest.mark.skipif("not HAS_MPL")
def test_rfi_config_finds_processed_scans(preprocessed, tmp_path, monkeypatch):
    """With -c, SDTrfistat reads the processed scans in intermediatedir and writes to productdir;
    --productdir moves the outputs, --outroot names the tables, which are then plotted."""
    tmp, config = preprocessed
    monkeypatch.chdir(tmp_path)
    main_rfistat(["-c", config])
    assert glob.glob(str(tmp / "products" / "*_rfi.hdf5"))
    assert glob.glob(str(tmp / "products" / "rfi*jpg"))
    assert os.listdir(tmp_path) == []
    assert not glob.glob(str(tmp / "raw" / "**" / "*.hdf5"), recursive=True)

    out = tmp_path / "out"
    mid = tmp / "products" / "intermediate"
    main_rfistat(
        ["-c", config, "--outroot", "bu", "--productdir", str(out), "--intermediatedir", str(mid)]
    )
    assert glob.glob(str(out / "bu_*_rfi.hdf5"))
    assert glob.glob(str(out / "rfi*jpg"))


def test_rfi_config_without_processed_scans(tmp_path):
    """With -c but no processed scans, SDTrfistat says to run SDTpreprocess first."""
    config = _make_config(tmp_path, RAW_FILES[:1])
    with pytest.warns(UserWarning, match="SDTpreprocess"):
        main_rfistat(["-c", config])


@pytest.mark.usefixtures("no_pyplot_figures")
def test_rfi_files_mpl(preprocessed, tmp_path, monkeypatch):
    """Given the processed files, SDTrfistat writes tables and plots in the current directory."""
    tmp, _ = preprocessed
    files = glob.glob(str(tmp / "products" / "intermediate" / "**" / "*.hdf5"), recursive=True)
    assert files
    monkeypatch.chdir(tmp_path)
    if HAS_MPL:
        main_rfistat(files)
        assert len(glob.glob("rfi*jpg")) > 0
    else:
        with pytest.warns(UserWarning, match="Matplotlib is required"):
            main_rfistat(files)

    assert len(glob.glob("*rfi.hdf5")) > 0
