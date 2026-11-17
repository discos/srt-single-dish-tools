import glob
import os
import shutil
from pathlib import Path

import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.logger import logging
from astropy.time import Time
from srttools.io import (
    bulk_change,
    locations,
    main_bulk_change,
    mkdir_p,
    print_obs_info_fitszilla,
    read_data_fitszilla,
)
from srttools.read_config import read_config
from srttools.scan import (
    HAS_MPL,
    Scan,
    clean_scan_using_variability,
    find_summary_file_in_dir,
    product_path_from_file_name,
)
from srttools.utils import compare_anything

try:
    import contextlib2 as contextlib

    FileNotFoundError = IOError
except ImportError:
    import contextlib


try:
    from sunpy.coordinates import frames

    HAS_SUNPY = True
except ImportError:
    HAS_SUNPY = False


@pytest.fixture
def logger():
    logger = logging.getLogger("Some.Logger")
    logger.setLevel(logging.INFO)

    return logger


class Test1_Sun:
    @classmethod
    def setup_class(klass):
        import os

        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, "data")

        klass.fname = os.path.abspath(os.path.join(klass.datadir, "sun_obs.fits"))
        klass.nosun = os.path.abspath(os.path.join(klass.datadir, "gauss_dec", "Dec0.fits"))

    @pytest.mark.skipif("not HAS_SUNPY")
    def test_read(self):
        scan = read_data_fitszilla(self.fname)
        assert "hpln" in scan.colnames

    @pytest.mark.skipif("HAS_SUNPY")
    def test_read_fails_if_no_sunpy(self):
        with pytest.raises(ValueError) as excinfo:
            _ = read_data_fitszilla(self.fname)
        assert "You need Sunpy to process Sun " in str(excinfo.value)

    @pytest.mark.skipif("not HAS_SUNPY")
    def test_read_no_sun(self):
        scan = read_data_fitszilla(self.nosun)
        assert "hpln" not in scan.colnames


class Test1_Scan:
    @classmethod
    def setup_class(klass):
        import os

        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, "data")

        klass.fname = os.path.abspath(os.path.join(klass.datadir, "gauss_dec", "Dec0.fits"))

        klass.config_file = os.path.abspath(os.path.join(klass.datadir, "test_config.ini"))

        config = read_config(klass.config_file)

        h5file = os.path.join(config["intermediatedir"], "gauss_dec", "Dec0.hdf5")
        if os.path.exists(h5file):
            os.unlink(h5file)

    @pytest.mark.usefixtures("no_pyplot_figures")
    def test_clean_large_data(self, caplog):
        """Test that large data sets needing pickling are handled correctly."""

        existing_pickle_files = glob.glob("*.p")
        clean_scan_using_variability(
            np.random.random((16000, 1000)),
            16,
            512,
            debug_file_format="jpg",
        )
        assert "data set is large" in caplog.text
        new_pickle_files = glob.glob("*.p")
        for epf in existing_pickle_files:
            if epf in new_pickle_files:
                new_pickle_files.remove(epf)
        # No new pickle files have survived the cleaning!
        assert len(new_pickle_files) == 0

    @pytest.mark.usefixtures("no_pyplot_figures")
    def test_clean_bad_thresh(self):
        """Test that an incorrect noise threshold produces no good data."""

        with pytest.warns(UserWarning, match="incorrect noise threshold"):
            clean_scan_using_variability(
                np.random.random((800, 1000)), 16, 512, debug_file_format="jpg", noise_threshold=-3
            )

    @pytest.mark.parametrize("bad_data", [0, np.nan, np.inf])
    def test_clean_bad_data(self, bad_data):
        """Test that an incorrect noise threshold produces no good data."""
        data = np.random.random((800, 1000))
        data[:, 500] = bad_data
        with pytest.warns(UserWarning, match="Bad channels in the data set"):
            clean_scan_using_variability(data, 16, 512)

    def test_bulk_change_hdr(self):
        dummyname = os.path.join(os.getcwd(), "dummyfile.fits")
        shutil.copyfile(self.fname, dummyname)
        with fits.open(dummyname) as hdul:
            header_dict1 = dict(hdul[0].header.items())
        bulk_change(dummyname, "0,header,DATE", "aaaaa")
        with fits.open(dummyname) as hdul:
            date_obs_new = hdul[0].header["DATE"]
            header_dict2 = dict(hdul[0].header.items())
        assert date_obs_new.strip() == "aaaaa"
        os.unlink(dummyname)
        # Check that all other keywords are unchanged
        header_dict1.pop("DATE")
        header_dict2.pop("DATE")
        assert compare_anything(header_dict1, header_dict2)

    def test_bulk_change_data(self):
        dummyname = os.path.join(os.getcwd(), "dummyfile.fits")
        shutil.copyfile(self.fname, dummyname)
        bulk_change(dummyname, "DATA TABLE,data,time", "0")
        with fits.open(dummyname) as hdul:
            time_new = np.array(hdul["DATA TABLE"].data["time"])
        assert np.all(time_new == 0)
        os.unlink(dummyname)

    def test_bulk_change_data_recursive(self):
        dummyname = os.path.join(os.getcwd(), "dummyfile.fits")
        shutil.copyfile(self.fname, dummyname)

        mkdir_p("blabla")
        second_dummyname = os.path.join("blabla", "dummyfile.fits")
        shutil.copyfile(self.fname, second_dummyname)
        main_bulk_change(["dummyfile.fits", "--apply-cal-mark", "--recursive"])
        with fits.open(dummyname) as hdul:
            flag_new = np.array(hdul["DATA TABLE"].data["flag_cal"])
        assert np.all(flag_new == 1)
        with fits.open(second_dummyname) as hdul:
            flag_new = np.array(hdul["DATA TABLE"].data["flag_cal"])
        assert np.all(flag_new == 1)
        os.unlink(dummyname)
        shutil.rmtree("blabla")

    def test_bulk_change_data_recursive_bad(self):
        dummyname = os.path.join(os.getcwd(), "dummyfile.fits")
        shutil.copyfile(self.fname, dummyname)
        with pytest.raises(Exception) as excinfo:
            main_bulk_change(
                [
                    dummyname,
                    "-k",
                    "DATA TABLE,data,time",
                    "-v",
                    "0",
                    "--recursive",
                ]
            )

        assert "Options recursive requires a file name" in str(excinfo.value)
        os.unlink(dummyname)

    def test_bulk_change_missing_key(self):
        dummyname = os.path.join(os.getcwd(), "dummyfile.fits")

        with pytest.raises(ValueError) as excinfo:
            main_bulk_change([dummyname])
        assert "What should I" in str(excinfo.value)

    def test_bulk_change_main(self):
        dummyname = os.path.join(os.getcwd(), "dummyfile.fits")
        shutil.copyfile(self.fname, dummyname)
        main_bulk_change([dummyname, "-k", "DATA TABLE,data,time", "-v", "0"])
        with fits.open(dummyname) as hdul:
            time_new = np.array(hdul["DATA TABLE"].data["time"])
        assert np.all(time_new == 0)
        os.unlink(dummyname)

    def test_temperatures_are_read(self):
        scan = read_data_fitszilla(os.path.join(self.datadir, "gauss_skydip", "skydip_mod.fits"))
        assert np.all(scan["Feed0_RCP-Temp"] > 0.0)

    def test_print_info(self, capsys):
        print_obs_info_fitszilla(self.fname)
        out, err = capsys.readouterr()
        assert "bandwidth" in out.lower()

    def test_repr(self):
        scan = Scan(self.fname)
        out = repr(scan)
        assert "scan from file" in out.lower()

    def test_print(self, capsys):
        scan = Scan(self.fname)
        print(scan)
        out, err = capsys.readouterr()
        assert "scan from file" not in out.lower()

    def test_scan(self):
        """Test that data are read."""

        scan = Scan(self.fname)
        scan.write("scan.hdf5", overwrite=True)
        scan2 = Scan("scan.hdf5")
        assert scan.meta == scan2.meta
        for col in scan.colnames:
            assert scan[col].meta == scan2[col].meta

    def test_scan_nofilt_executes(self):
        """Test that data are read."""

        scan = Scan(self.fname, nofilt=True)

    def test_scan_from_table(self):
        """Test that data are read."""
        from astropy.table import Table

        scan = Scan(self.fname)
        scan.write("scan.hdf5", overwrite=True)
        table = Table.read("scan.hdf5")
        scan_from_table = Scan(table)
        for c in scan.columns:
            assert np.all(scan_from_table[c] == scan[c])
        for m in scan_from_table.meta.keys():
            assert scan_from_table.meta[m] == scan.meta[m]

    @pytest.mark.parametrize("fname", ["med_data.fits", "srt_data_tp_multif.fits"])
    def test_coordinate_conversion_works(self, fname):
        scan = Scan(os.path.join(self.datadir, fname), norefilt=False)
        obstimes = Time(scan["time"] * u.day, format="mjd", scale="utc")
        idx = 1 if "_multif" in fname else 0
        ref_coords = SkyCoord(
            ra=scan["ra"][:, idx],
            dec=scan["dec"][:, idx],
            obstime=obstimes,
            location=locations[scan.meta["site"]],
        )
        altaz = ref_coords.altaz

        diff = np.abs((altaz.az.to(u.rad) - scan["az"][:, idx]).to(u.arcsec).value)
        assert np.all(diff < 1)
        diff = np.abs((altaz.alt.to(u.rad) - scan["el"][:, idx]).to(u.arcsec).value)
        assert np.all(diff < 1)

    def test_bad_nchan_detected(self):
        with pytest.raises(Exception) as excinfo:
            scan = Scan(os.path.join(self.datadir, "srt_chans_bad.fits"))

        assert "with channel subdivision:" in str(excinfo.value)

    def test_simple_in_stokes(self):
        with pytest.warns(UserWarning) as record:
            scan = Scan(os.path.join(self.datadir, "srt_pol_bad.fits"), norefilt=False)

        assert np.any(["contain polarization information" in r.message.args[0] for r in record])

    @classmethod
    def teardown_class(klass):
        """Cleanup."""
        if os.path.exists("scan.hdf5"):
            os.unlink("scan.hdf5")
        for f in glob.glob(os.path.join(klass.datadir, "*.hdf5")):
            os.unlink(f)


class Test2_Scan:
    @classmethod
    def setup_class(klass):
        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, "data")

        klass.fname = os.path.abspath(
            os.path.join(klass.datadir, "spectrum", "roach_template.fits")
        )

        klass.config_file = os.path.abspath(os.path.join(klass.datadir, "spectrum.ini"))
        config = read_config(klass.config_file)

        h5file = os.path.join(config["intermediatedir"], "spectrum", "roach_template.hdf5")
        if os.path.exists(h5file):
            os.unlink(h5file)
        klass.config = config

    @pytest.mark.usefixtures("no_pyplot_figures")
    def test_scan(self):
        """Test that data are read."""

        scan = Scan(self.fname, debug=True)

        scan.write("scan.hdf5", overwrite=True)
        scan.baseline_subtract("rough", plot=True)

    @pytest.mark.parametrize(
        "fname",
        [
            "srt_data.fits",
            "srt_data_roach_polar.fits",
            "srt_data_xarcos.fits",
            "new_sardara.fits",
            "new_sardara.fits3",
        ],
    )
    @pytest.mark.usefixtures("no_pyplot_figures")
    def test_scan_loading(self, fname):
        """Test that data are read."""

        Scan(os.path.join(self.datadir, "spectrum", fname), plot=True)
        assert os.path.exists(
            os.path.join(
                self.config["intermediatedir"],
                "spectrum",
                fname.replace(".fits", "") + ".hdf5",
            )
        )

    def test_scan_loading_badchans(self):
        """Test that data are read."""
        fname = "srt_data_roach_polar.fits0"
        with pytest.warns(UserWarning) as record:
            Scan(os.path.join(self.datadir, "spectrum", fname))
        assert np.any(
            ["No good channels found. A problem with the " in r.message.args[0] for r in record]
        )
        assert np.any(["Bad channels in the data set" in r.message.args[0] for r in record])
        assert os.path.exists(
            os.path.join(
                self.config["intermediatedir"],
                "spectrum",
                fname.replace(".fits", "") + ".hdf5",
            )
        )

    def test_scan_baseline_unknown(self):
        """Test that data are read."""

        scan = Scan(self.fname)
        scan.write("scan.hdf5", overwrite=True)
        with pytest.raises(ValueError):
            scan.baseline_subtract("asdfgh", plot=True)

    def test_scan_write_other_than_hdf5_raises(self):
        """Test that data are read."""

        scan = Scan(self.fname)
        with pytest.raises(TypeError):
            scan.write("scan.fits", overwrite=True)
        with pytest.raises(TypeError):
            scan.write("scan.json", overwrite=True)
        with pytest.raises(TypeError):
            scan.write("scan.csv", overwrite=True)

    def test_scan_clean_and_splat(self):
        """Test that data are read."""

        scan = Scan(self.fname)
        scan.meta["filtering_factor"] = 0.7
        with pytest.warns(UserWarning) as record:
            scan.clean_and_splat()
            assert np.any(
                ["Don't use filtering factors > 0.5" in r.message.args[0] for r in record]
            )

    @pytest.mark.parametrize("fname", ["srt_data.fits", "srt_data_roach_polar.fits"])
    def test_coordinate_conversion_works(self, fname):
        scan = Scan(
            os.path.join(self.datadir, "spectrum", fname),
            norefilt=False,
            debug=True,
        )
        obstimes = Time(scan["time"] * u.day, format="mjd", scale="utc")
        idx = 1 if "_multif" in fname else 0

        # Tolerance: +- 1 second of observation, or sample time, whichever is
        # larger
        sampletime = np.max([scan[ch].meta["sample_time"].value for ch in scan.chan_columns()])
        sampletime = np.max([sampletime, 1]) * u.s

        ref_coords0 = SkyCoord(
            ra=scan["ra"][:, idx],
            dec=scan["dec"][:, idx],
            obstime=obstimes - sampletime,
            location=locations[scan.meta["site"]],
        )
        ref_coords1 = SkyCoord(
            ra=scan["ra"][:, idx],
            dec=scan["dec"][:, idx],
            obstime=obstimes + sampletime,
            location=locations[scan.meta["site"]],
        )

        altaz0 = ref_coords0.altaz
        altaz1 = ref_coords1.altaz

        az0, az1 = altaz0.az.to(u.rad).value, altaz1.az.to(u.rad).value
        tol = (30 * u.arcsec).to(u.rad).value
        good0 = (scan["az"][:, idx] >= az0 - tol) & (scan["az"][:, idx] <= az1 + tol)
        good1 = (scan["az"][:, idx] <= az1 + tol) & (scan["az"][:, idx] >= az0 - tol)
        assert np.all(good0 | good1)

        el0, el1 = altaz0.alt.to(u.rad).value, altaz1.alt.to(u.rad).value
        tol = (30 * u.arcsec).to(u.rad).value
        good0 = (scan["el"][:, idx] >= el0 - tol) & (scan["el"][:, idx] <= el1 + tol)
        good1 = (scan["el"][:, idx] <= el1 + tol) & (scan["el"][:, idx] >= el0 - tol)
        assert np.all(good0 | good1)

    @classmethod
    def teardown_class(klass):
        """Cleanup."""
        with contextlib.suppress(FileNotFoundError):
            os.unlink("scan.hdf5")
            for f in glob.glob(os.path.join(klass.datadir, "spectrum", "*.jpg")):
                os.unlink(f)
            for f in glob.glob(os.path.join(klass.datadir, "spectrum", "*.hdf5")):
                os.unlink(f)
            shutil.rmtree(os.path.join(klass.datadir, "out_spectrum_test"))


class TestSummary:
    @classmethod
    def setup_class(cls):
        cls.testdir = "tmp-fake-dir"
        os.makedirs(cls.testdir)

    def test_summary_nofiles(self):
        with pytest.warns(UserWarning, match="No summary file"):
            find_summary_file_in_dir(self.testdir)

    def test_summary_onefile(self):
        fname = os.path.join(self.testdir, "summary.fits")
        Path(fname).touch()
        assert fname == find_summary_file_in_dir(self.testdir)

    def test_summary_multiple(self):
        Path(os.path.join(self.testdir, "summary.fits")).touch()
        Path(os.path.join(self.testdir, "Sum_asdfhasldfjhal.fits")).touch()
        with pytest.raises(ValueError, match="Multiple"):
            find_summary_file_in_dir(self.testdir)

    def teardown_class(cls):
        shutil.rmtree(cls.testdir)


class TestScanUtils:
    def test_product_path_mirrors_datadir(self):
        """Products of files inside datadir mirror the directory structure below datadir only."""
        d = os.path.join(os.sep + "a", "b", ".tox", "data")
        f = os.path.join(d, "obs1", "sub", "pippo.fits")
        p = os.path.abspath("prod")
        path, fname = product_path_from_file_name(f, datadir=d, productdir=p)
        assert fname == "pippo.fits"
        assert path == os.path.join(p, "obs1", "sub")
        # Relative datadir
        path, _ = product_path_from_file_name(f, datadir=os.path.relpath(d), productdir=p)
        assert path == os.path.join(p, "obs1", "sub")

    def test_product_path_outside_datadir(self):
        """Products of files outside datadir go to a folder named as their own, not to a
        copy of their whole absolute path."""
        f = os.path.join(os.sep + "x", "y", "obs2", "pippo.fits")
        d = os.path.join(os.sep + "a", "b", "data")
        p = os.path.abspath("prod")
        path, _ = product_path_from_file_name(f, datadir=d, productdir=p)
        assert path == os.path.join(p, "obs2")

    @pytest.mark.skipif(os.name != "nt", reason="Drive letters exist only on Windows")
    def test_product_path_other_drive(self):
        """Data on a different drive than datadir do not raise an error."""
        path, _ = product_path_from_file_name(
            "D:\\data\\obs3\\pippo.fits", datadir="C:\\work", productdir="C:\\prod"
        )
        assert path == "C:\\prod\\obs3"
