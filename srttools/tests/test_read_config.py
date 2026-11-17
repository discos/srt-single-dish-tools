import os

import numpy as np
import pytest

import astropy.units as u
from srttools.read_config import DEFAULT_PRODUCTDIR, read_config


class TestConfig:
    @classmethod
    def setup_class(cls):
        cls.curdir = os.path.abspath(os.path.dirname(__file__))
        cls.datadir = os.path.join(cls.curdir, "data")

    def test_read_config(self):
        """Test that config file are read."""

        fname = os.path.join(self.datadir, "test_config.ini")

        config = read_config(fname)

        np.testing.assert_almost_equal(config["pixel_size"].to(u.rad).value, np.radians(0.5 / 60))
        assert config["interpolation"] == "spline"

    def test_read_incomplete_config(self):
        """Test that config file are read."""
        fname = os.path.join(self.datadir, "test_config_incomplete.ini")

        config = read_config(fname)

        assert config["pixel_size"] == "auto"
        assert config["interpolation"] == "linear"


CONFIG = "[local]\ndatadir : {datadir}\n{extra}[analysis]\nlist_of_directories :\n    bubu\n"


@pytest.mark.parametrize("extra", ["", "productdir : None\n"])
def test_missing_productdir_defaults_to_current_directory(tmp_path, monkeypatch, extra):
    """Without a productdir, products go to ./srttools_products, with a warning; never next to the data."""
    datadir = tmp_path / "raw"
    datadir.mkdir()
    workdir = tmp_path / "work"
    workdir.mkdir()
    config_file = tmp_path / "config.ini"
    config_file.write_text(CONFIG.format(datadir=datadir, extra=extra))
    monkeypatch.chdir(workdir)
    with pytest.warns(UserWarning, match="NOT next to the raw data"):
        config = read_config(str(config_file))
    assert config["productdir"] == os.path.join(os.path.realpath(workdir), DEFAULT_PRODUCTDIR)


def test_explicit_productdir_is_kept_without_warning(tmp_path, recwarn):
    """A productdir given in the config file is used as it is (relative to the config file)."""
    config_file = tmp_path / "config.ini"
    config_file.write_text(CONFIG.format(datadir=tmp_path, extra="productdir : out\n"))
    config = read_config(str(config_file))
    assert config["productdir"] == str(tmp_path / "out")
    assert not [w for w in recwarn if "productdir" in str(w.message)]


def test_intermediatedir_defaults_to_subdirectory_of_productdir(tmp_path):
    """Without an intermediatedir, intermediate products go to <productdir>/intermediate."""
    config_file = tmp_path / "config.ini"
    config_file.write_text(CONFIG.format(datadir=tmp_path, extra="productdir : out\n"))
    config = read_config(str(config_file))
    assert config["intermediatedir"] == str(tmp_path / "out" / "intermediate")


def test_explicit_intermediatedir_is_relative_to_config_file(tmp_path):
    """An intermediatedir in the config file is used as it is, relative to the config file."""
    config_file = tmp_path / "config.ini"
    extra = "productdir : out\nintermediatedir : scans\n"
    config_file.write_text(CONFIG.format(datadir=tmp_path, extra=extra))
    config = read_config(str(config_file))
    assert config["productdir"] == str(tmp_path / "out")
    assert config["intermediatedir"] == str(tmp_path / "scans")


def test_directory_overrides_on_cached_config(tmp_path, monkeypatch):
    """Overriding productdir moves the default intermediatedir too; an explicit one stays."""
    monkeypatch.chdir(tmp_path)
    config_file = tmp_path / "config.ini"
    config_file.write_text(CONFIG.format(datadir=tmp_path, extra="productdir : out\n"))
    read_config(str(config_file))
    config = read_config(str(config_file), productdir="new")
    assert config["productdir"] == str(tmp_path / "new")
    assert config["intermediatedir"] == str(tmp_path / "new" / "intermediate")
    config = read_config(str(config_file), intermediatedir="mid")
    assert config["intermediatedir"] == str(tmp_path / "mid")
    config = read_config(str(config_file), productdir="newer")
    assert config["intermediatedir"] == str(tmp_path / "mid")


def test_upgrade_config_adds_intermediatedir(tmp_path):
    """SDTupgradeconfig keeps old processed scans usable: intermediatedir = old productdir, .bak kept."""
    from srttools.read_config import main_upgrade_config

    config_file = tmp_path / "config.ini"
    text = CONFIG.format(datadir=tmp_path, extra="; my comment\nproductdir : out\n")
    config_file.write_text(text)
    main_upgrade_config([str(config_file)])
    assert (tmp_path / "config.ini.bak").read_text() == text
    assert "; my comment" in config_file.read_text()
    config = read_config(str(config_file))
    assert config["intermediatedir"] == config["productdir"] == str(tmp_path / "out")
    # A second run changes nothing
    upgraded = config_file.read_text()
    main_upgrade_config([str(config_file)])
    assert config_file.read_text() == upgraded


def test_old_product_layout_warns(tmp_path):
    """Processed scans of older versions in productdir, without intermediatedir: point to the converter."""
    (tmp_path / "out" / "scan1").mkdir(parents=True)
    (tmp_path / "out" / "scan1" / "subscan.hdf5").write_text("")
    config_file = tmp_path / "config.ini"
    config_file.write_text(CONFIG.format(datadir=tmp_path, extra="productdir : out\n"))
    with pytest.warns(UserWarning, match="SDTupgradeconfig"):
        read_config(str(config_file))
