import glob
import os
import shutil
import subprocess as sp
import warnings

import numpy as np
import pytest

from astropy.logger import logging
from astropy.table import Table
from srttools.calibration import (
    CalibratorTable,
    _get_flux_quantity,
    _safe_fit_baseline_plus_bell,
    _treat_scan,
    main_cal,
    main_lcurve,
)
from srttools.io import mkdir_p
from srttools.read_config import output_file_name, read_config
from srttools.scan import list_scans
from srttools.simulate import _2d_gauss, sim_crossscans
from srttools.utils import HAS_MPL, compare_strings

try:
    from tqdm import tqdm
except ImportError:

    def tqdm(x):
        return x


@pytest.fixture
def logger():
    logger = logging.getLogger("Some.Logger")
    logger.setLevel(logging.DEBUG)

    return logger


np.random.seed(124137)


class TestCalibration:
    @classmethod
    def setup_class(klass):
        import os

        klass.curdir = os.path.dirname(__file__)
        klass.datadir = os.path.join(klass.curdir, "data")

        klass.config_file = os.path.abspath(os.path.join(klass.datadir, "calibrators.ini"))

        klass.config_file_empty = os.path.abspath(
            os.path.join(klass.datadir, "calibrators_nocal.ini")
        )

        klass.config = read_config(klass.config_file)
        klass.simdir = klass.caldir = os.path.join(klass.datadir, "sim")

        klass.caldir = os.path.join(klass.simdir, "calibration")
        klass.caldir2 = os.path.join(klass.simdir, "calibration2")
        klass.caldir3 = os.path.join(klass.simdir, "calibration_bad")
        klass.crossdir = os.path.join(klass.simdir, "crossscans")

        klass.scan_list = (
            list_scans(klass.caldir, ["./"])
            + list_scans(klass.caldir2, ["./"])
            + list_scans(klass.crossdir, ["./"])
        )

        klass.scan_list.sort()
        caltable = CalibratorTable()
        caltable.from_scans(klass.scan_list, debug=False)
        caltable.update()

        klass.calfile = os.path.join(klass.curdir, "test_calibration.hdf5")
        for calfile in [klass.calfile, klass.calfile.replace("hdf5", "csv")]:
            if os.path.exists(calfile):
                os.remove(calfile)
        caltable.write(klass.calfile, overwrite=True)
        caltable.write(klass.calfile.replace("hdf5", "csv"))

    def test_0_prepare(self):
        pass

    def test_script_is_installed(self):
        sp.check_call("SDTcal -h".split(" "))

    def test_check_not_empty(self):
        caltable = CalibratorTable()
        assert not caltable.check_not_empty()

    def test_check_up_to_date_empty_return_false(self):
        caltable = CalibratorTable()
        assert not caltable.check_up_to_date()

    def test_calibrate_empty_return_none(self):
        caltable = CalibratorTable()
        assert caltable.calibrate() is None

    def test_update_empty_return_none(self):
        caltable = CalibratorTable()
        assert caltable.update() is None

    def test_get_fluxes_empty_return_none(self):
        caltable = CalibratorTable()
        assert caltable.get_fluxes() is None

    @pytest.mark.usefixtures("no_pyplot_figures")
    def test_check_class(self):
        caltable = CalibratorTable()
        caltable.from_scans(list_scans(self.caldir, ["./"]), debug=True)
        with pytest.warns(UserWarning):
            caltable.check_up_to_date()

    def test_check_class_from_file(self):
        caltable = CalibratorTable.read(self.calfile)
        assert caltable.check_up_to_date()

    def test_Jy_over_counts_and_back(self):
        caltable = CalibratorTable.read(self.calfile)
        Jc, Jce = caltable.Jy_over_counts(channel="Feed0_LCP")
        Cj, Cje = caltable.counts_over_Jy(channel="Feed0_LCP")
        np.testing.assert_allclose(Jc, 1 / Cj)

    def test_zero_amplitude_fits_are_ignored(self):
        """Fits with zero amplitude (infinite Flux/Counts) neither crash nor change the calibration."""
        caltable = CalibratorTable.read(self.calfile)
        caltable.update()
        calib = np.flatnonzero((caltable["Chan"] == "Feed0_LCP") & (caltable["Flux"] > 0))
        zeroed = calib[:3]
        reference = CalibratorTable(caltable.copy())
        reference.remove_rows(zeroed)
        reference.update()

        caltable["Counts"][zeroed] = 0
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            caltable.update()
        assert np.all(np.isinf(caltable["Flux/Counts"][zeroed]))
        el = np.radians(45)
        np.testing.assert_allclose(
            caltable.Jy_over_counts(channel="Feed0_LCP", elevation=el),
            reference.Jy_over_counts(channel="Feed0_LCP", elevation=el),
        )
        np.testing.assert_allclose(
            caltable.Jy_over_counts_rough(channel="Feed0_LCP"),
            reference.Jy_over_counts_rough(channel="Feed0_LCP"),
        )

    def test_two_calibrator_fits_give_a_constant_calibration(self):
        """A channel with only two calibrator fits (e.g. one cross scan) cannot constrain a
        line in elevation: the calibration is the constant from Jy_over_counts_rough."""
        caltable = CalibratorTable.read(self.calfile)
        caltable.update()
        calib = np.flatnonzero((caltable["Chan"] == "Feed0_LCP") & (caltable["Flux"] > 0))
        caltable.remove_rows(calib[2:])
        caltable.update()
        el = np.radians(45)
        np.testing.assert_allclose(
            caltable.Jy_over_counts(channel="Feed0_LCP", elevation=el),
            caltable.Jy_over_counts_rough(channel="Feed0_LCP"),
        )

    @pytest.mark.parametrize("scatter", [0.01, 0.3])
    def test_Jy_over_counts_rough_error_includes_scatter(self, scatter):
        """With formal errors of 1%, the error of the constant calibration matches the
        scatter of the calibrator fits, also when they scatter by 30% (it used to be
        the formal error, 30x too small)."""
        caltable = CalibratorTable.read(self.calfile)
        caltable.update()
        rows = np.flatnonzero((caltable["Chan"] == "Feed0_LCP") & (caltable["Flux"] > 0))
        rng = np.random.default_rng(42)
        true = np.mean(caltable["Flux/Counts"][rows])
        caltable["Flux/Counts"][rows] = true * (1 + rng.normal(0, scatter, rows.size))
        caltable["Flux/Counts Err"][rows] = true * 0.01
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            fc, fce = caltable.Jy_over_counts_rough(channel="Feed0_LCP")
        expected = true * scatter / np.sqrt(rows.size)
        assert 0.6 * expected < fce < 2 * expected
        assert abs(fc - true) < 4 * expected

    def test_Jy_over_counts_rough_one_bad_value(self, logger, caplog):
        caltable = CalibratorTable.read(self.calfile)

        flux_quantity = _get_flux_quantity("Jy/beam")
        good = ~np.isnan(caltable[flux_quantity + "/Counts"])
        good = good & (caltable["Chan"] == "Feed0_LCP")
        assert np.count_nonzero(good) > 1
        std = np.std(np.diff(caltable[flux_quantity + "/Counts"][good]))
        assert std > 0
        firstidx = np.where(good)[0][0]
        caltable[flux_quantity + "/Counts"][firstidx] = 20000

        assert caltable[flux_quantity + "/Counts"][firstidx] == 20000

        Jc, Jce = caltable.Jy_over_counts_rough(channel="Feed0_LCP", map_unit="Jy/beam")

        assert "Outliers: " in caplog.text
        Cj, Cje = caltable.counts_over_Jy(channel="Feed0_LCP")
        np.testing.assert_allclose(Jc, 1 / Cj)

    def test_bad_file_missing_key(self, logger, caplog):
        caltable = CalibratorTable()
        caltable.from_scans([os.path.join(self.config["datadir"], "calibrators", "summary.fits")])
        assert "Missing key" in caplog.text

    def test_bad_file_generic_error(self, logger, caplog):
        caltable = CalibratorTable()

        caltable.from_scans([os.path.join(self.config["datadir"], "calibrators", "bubu.fits")])
        assert "Error while processing" in caplog.text

    def test_calibration_counts(self):
        """Simple calibration from scans."""

        caltable = CalibratorTable.read(self.calfile)
        caltable = caltable[compare_strings(caltable["Source"], "DummyCal")]
        caltable_0 = caltable[compare_strings(caltable["Chan"], "Feed0_LCP")]
        assert np.all(np.abs(caltable_0["Counts"] - 100.0) < 3 * caltable_0["Counts Err"])
        caltable_1 = caltable[compare_strings(caltable["Chan"], "Feed0_RCP")]
        assert np.all(np.abs(caltable_1["Counts"] - 80.0) < 3 * caltable_1["Counts Err"])

    def test_calibration_width(self):
        """Simple calibration from scans."""

        caltable = CalibratorTable.read(self.calfile)
        caltable0 = caltable[compare_strings(caltable["Chan"], "Feed0_LCP")]
        assert np.all(np.abs(caltable0["Width"] - 2.5 / 60.0) < 5 * caltable0["Width Err"])
        caltable1 = caltable[compare_strings(caltable["Chan"], "Feed0_RCP")]
        assert np.all(np.abs(caltable1["Width"] - 2.5 / 60.0) < 5 * caltable1["Width Err"])

        beam, beam_err = caltable.beam_width(channel="Feed0_LCP")
        assert np.all(beam - np.radians(2.5 / 60) < 3 * beam_err)

    @pytest.mark.skipif("not HAS_MPL")
    @pytest.mark.usefixtures("no_pyplot_figures")
    def test_calibration_plot_two_cols(self):
        """Simple calibration from scans."""

        caltable = CalibratorTable.read(self.calfile)
        # The values need to be positive
        caltable["RA err"] = np.abs(caltable["RA err"])
        caltable.plot_two_columns(
            "RA",
            "Flux/Counts",
            xerrcol="RA err",
            yerrcol="Flux/Counts Err",
            test=True,
        )

    @pytest.mark.skipif("not HAS_MPL")
    @pytest.mark.usefixtures("no_pyplot_figures")
    def test_calibration_show(self):
        """Simple calibration from scans."""

        caltable = CalibratorTable.read(self.calfile)

        caltable.show()

    def test_calibrated_crossscans(self):
        caltable = CalibratorTable.read(self.calfile)
        dummy_flux, dummy_flux_err = caltable.calculate_src_flux(
            source="DummySrc", channel="Feed0_LCP"
        )
        assert (dummy_flux[0] - 0.52) < dummy_flux_err[0] * 3

    def test_residual_std_is_stored(self):
        """Each fit stores the rms of its residuals, used to find RFI-contaminated subscans."""
        caltable = CalibratorTable.read(self.calfile)
        assert np.all(caltable["Resid Std"] > 0)

    def test_calculate_src_flux_ignores_outliers(self):
        """Fits flagged as outliers do not enter the source flux."""
        caltable = CalibratorTable.read(self.calfile)
        ref, _ = caltable.calculate_src_flux(source="DummySrc", channel="Feed0_LCP")
        row = np.flatnonzero(
            compare_strings(caltable["Source"], "DummySrc")
            & compare_strings(caltable["Chan"], "Feed0_LCP")
        )[0]
        caltable["Counts"][row] *= 100
        caltable["Outlier"][row] = True
        flux, _ = caltable.calculate_src_flux(source="DummySrc", channel="Feed0_LCP")
        good = np.ones(len(caltable), dtype=bool)
        good[row] = False
        ref, _ = CalibratorTable(caltable[good]).calculate_src_flux(
            source="DummySrc", channel="Feed0_LCP"
        )
        assert np.isclose(flux[0], ref[0])

    def test_check_consistency_fails_with_bad_data(self):
        scan_list = (
            list_scans(self.caldir, ["./"])
            + list_scans(self.caldir2, ["./"])
            + list_scans(self.caldir3, ["./"])
            + list_scans(self.crossdir, ["./"])
        )

        scan_list.sort()
        caltable = CalibratorTable()
        caltable.from_scans(scan_list)
        caltable.update()
        res = caltable.check_consistency(channel="Feed0_LCP")
        assert not np.all(res)

    @pytest.mark.parametrize("chan", ["Feed0_LCP", "Feed0_RCP"])
    def test_check_consistency_chan_by_chan(self, chan):
        caltable = CalibratorTable.read(self.calfile)
        res = caltable.check_consistency(channel=chan)
        assert np.all(res)

    def test_check_consistency_all(self):
        caltable = CalibratorTable.read(self.calfile)
        res = caltable.check_consistency()
        assert np.all(res)

    @pytest.mark.skipif("not HAS_MPL")
    def test_sdtcal_with_calfile(self):
        plotfile = self.calfile.replace("hdf5", "jpg")
        if os.path.exists(plotfile):
            os.unlink(plotfile)
        main_cal([self.calfile])
        assert os.path.exists(plotfile)
        os.unlink(plotfile)

    @pytest.mark.skipif("not HAS_MPL")
    def test_sdtcal_show_with_config(self):
        main_cal(("-c " + self.config_file + " --check --show").split(" "))
        outfile = output_file_name(self.config_file, "_cal.hdf5")
        plotfile = outfile.replace(".hdf5", ".jpg")
        assert os.path.exists(outfile)
        assert os.path.exists(plotfile)
        os.unlink(plotfile)

        # Reload unfiltered cal, a cache in intermediatedir
        unfilt = output_file_name(self.config_file, "_cal_unfilt.hdf5", intermediate=True)
        assert os.path.exists(unfilt)
        os.unlink(outfile)
        main_cal(("-c " + self.config_file).split(" "))
        assert os.path.exists(outfile)

    def test_sdtcal_with_sample_config(self):
        if os.path.exists("sample_config_file.ini"):
            os.unlink("sample_config_file.ini")
        main_cal(["--sample-config"])
        assert os.path.exists("sample_config_file.ini")

    def test_sdtcal_no_config(self):
        # ValueError("Please specify the config file!")
        with pytest.raises(ValueError) as excinfo:
            main_cal([])
            assert "Please specify the config file!" in str(excinfo.value)

    def test_sdtcal_no_config_dir(self):
        # ValueError("No calibrators specified in config file")
        with pytest.raises(ValueError) as excinfo:
            main_cal(["-c", self.config_file_empty])
            assert "No calibrators specified in config file" in str(excinfo.value)

    def test_lcurve_with_single_source(self):
        main_lcurve([self.calfile, "-s", "DummySrc"])
        assert os.path.exists("DummySrc.csv")
        os.unlink("DummySrc.csv")

    def test_lcurve_with_all_sources(self):
        main_lcurve(["-c", self.config_file])
        assert os.path.exists("DummySrc.csv")
        assert os.path.exists("DummyCal.csv")
        assert os.path.exists("DummyCal2.csv")

    @classmethod
    def teardown_class(klass):
        """Clean up the mess."""
        if HAS_MPL:
            for f in glob.glob("*.jpg"):
                os.unlink(f)
        # Processed scans and fit plots are all below productdir
        shutil.rmtree(klass.config["productdir"], ignore_errors=True)
        if os.path.exists(klass.calfile):
            os.remove(klass.calfile)


def test_safe_fit_baseline_plus_bell_does_not_crash():
    """When the Gaussian fit raises, the safe wrapper warns and returns None. A strong spike
    makes it raise only with some astropy/scipy versions: x and y of different lengths always do."""
    x = np.linspace(-0.24, 0.24, 475)
    y = np.zeros(474)
    with pytest.warns(UserWarning, match="Fit failed in test_scan"):
        model, fit_info = _safe_fit_baseline_plus_bell(x, y, "test_scan")
    assert model is None
    assert fit_info is None


def _fit_table(nrows=30, seed=0):
    """A CalibratorTable with only the columns used by flag_outliers, as in an old file."""
    rng = np.random.default_rng(seed)
    t = Table()
    t["Dir"] = ["obs0"] * nrows
    t["Source"] = ["Src"] * nrows
    t["Chan"] = ["Feed0_LCP"] * nrows
    t["Scan Type"] = ["RA>", "Dec>"] * (nrows // 2)
    t["Counts"] = rng.normal(100, 1, nrows)
    t["Width"] = rng.normal(0.04, 0.001, nrows)
    t["Resid Std"] = rng.normal(1, 0.03, nrows)
    return CalibratorTable(t)


def test_flag_outliers():
    """Fits with anomalous residual rms, amplitude or width are flagged, and only those."""
    caltable = _fit_table()
    caltable["Resid Std"][5] = 5
    caltable["Counts"][12] = 200
    caltable["Width"][20] = 0.2
    with pytest.warns(UserWarning, match="outlier"):
        caltable.flag_outliers(nsigma=5)
    assert sorted(np.flatnonzero(caltable["Outlier"])) == [5, 12, 20]


def test_flag_outliers_old_table_without_residuals():
    """Tables from older versions (no "Resid Std") are still checked on amplitude."""
    caltable = _fit_table()
    caltable.remove_column("Resid Std")
    caltable = CalibratorTable(caltable)
    caltable["Counts"][12] = 200
    with pytest.warns(UserWarning, match="outlier"):
        caltable.flag_outliers(nsigma=5)
    assert sorted(np.flatnonzero(caltable["Outlier"])) == [12]


def test_flag_outliers_small_groups_untouched():
    """Groups with too few fits for a robust spread are never flagged."""
    caltable = _fit_table(nrows=8)
    caltable["Counts"][2] = 200
    caltable.flag_outliers(nsigma=5)
    assert not np.any(caltable["Outlier"])


def test_treat_scan_multifeed_fits_each_channel_with_its_own_feed(tmp_path, monkeypatch):
    """On a 7-feed scan, each channel is fitted once, with its own feed's coordinates."""
    datadir = os.path.join(os.path.dirname(__file__), "data")
    fname = str(tmp_path / "srt_data_tp_multif.fits")
    shutil.copyfile(os.path.join(datadir, "srt_data_tp_multif.fits"), fname)
    # Without a config, the default ./srttools_products goes to tmp_path
    monkeypatch.chdir(tmp_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        retval, rows = _treat_scan(fname, max_offset=None)
    assert retval
    chans = [r[4] for r in rows]
    feeds = [r[5] for r in rows]
    assert len(chans) == len(set(chans))
    assert feeds == [int(ch[4]) for ch in chans]


def test_treat_scan_skips_subscans_missing_the_source(tmp_path, monkeypatch):
    """A calibrator subscan whose feeds all pass >10 beams from the source (here, a
    map subscan) gives no fits: fitting noise would give a meaningless Flux/Counts."""
    datadir = os.path.join(os.path.dirname(__file__), "data")
    fname = str(tmp_path / "srt_data_tp_multif.fits")
    shutil.copyfile(os.path.join(datadir, "srt_data_tp_multif.fits"), fname)
    # Without a config, the default ./srttools_products goes to tmp_path
    monkeypatch.chdir(tmp_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        retval, rows = _treat_scan(fname)
    assert retval
    assert rows == []


def test_calibrator_fits_write_only_below_productdir(tmp_path, monkeypatch):
    """SDTcal fits honour the config file: plots and temporary files go below
    its intermediatedir, never next to the raw data or into the current directory."""
    datadir = tmp_path / "data"
    datadir.mkdir()
    fname = str(datadir / "srt_data_tp_multif.fits")
    shutil.copyfile(
        os.path.join(os.path.dirname(__file__), "data", "srt_data_tp_multif.fits"), fname
    )
    config_file = tmp_path / "cal.ini"
    config_file.write_text(
        f"[local]\ndatadir : {datadir}\nproductdir : {tmp_path / 'out'}\n"
        "[analysis]\nlist_of_directories :\n    ./\n"
    )
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    monkeypatch.chdir(cwd)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        CalibratorTable().from_scans([fname], config_file=str(config_file))
        retval, _ = _treat_scan(fname, plot=True, max_offset=None, config_file=str(config_file))
    assert retval
    assert os.listdir(cwd) == []
    assert os.listdir(datadir) == ["srt_data_tp_multif.fits"]
    plotdirs = glob.glob(
        str(tmp_path / "out" / "intermediate" / "**" / "*_scanfit"), recursive=True
    )
    assert len(plotdirs) == 1 and os.listdir(plotdirs[0])
