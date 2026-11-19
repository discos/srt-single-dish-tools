"""End-to-end test on real data: Cas A mapped with Medicina, calibrated on NGC7027.

Medicina 32 m, 2022-12-13, 8.18 GHz, two channels (RCP, LCP). The raw
subscans are in ``data/medicina``. All products go to a temporary productdir.
"""

import glob
import os
import warnings

import numpy as np
import pytest

from srttools.calibration import CalibratorTable
from srttools.imager import ScanSet

DATADIR = os.path.join(os.path.dirname(__file__), "data", "medicina")
CALDIR = "20221213-154904-22-22-NGC7027"
MAPDIR = "20221213-162927-22-22-cas_a"
PIXEL = 1.5  # arcmin, larger than the subscan spacing (1.25 arcmin)
R_AP, R_IN, R_OUT = 12.0, 20.0, 28.0  # aperture and background annulus, arcmin
# Perley & Butler 2017, extrapolated to 8.18 GHz and to Dec 2022 (secular
# decline). Uncertain by several percent.
CASA_EXPECTED = 375.0


def _list_data_files():
    return sorted(glob.glob(os.path.join(DATADIR, "**"), recursive=True))


@pytest.fixture(scope="module")
def casa(tmp_path_factory):
    """Calibrate on NGC7027 and image Cas A in Jy/pixel, from scratch."""
    tmp = tmp_path_factory.mktemp("medicina")
    config = str(tmp / "casa.ini")
    with open(config, "w") as fobj:
        print(
            f"""[local]
workdir = {tmp}
datadir = {DATADIR}
productdir = {tmp / "products"}

[analysis]
projection = ARC
interpolation = spline
list_of_directories =
    {MAPDIR}/
calibrator_directories =
    {CALDIR}/
noise_threshold = 5
pixel_size = {PIXEL}
goodchans =
""",
            file=fobj,
        )
    files_before = _list_data_files()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        caltable = CalibratorTable()
        caltable.from_scans(
            sorted(glob.glob(os.path.join(DATADIR, CALDIR, "*.fits"))), config_file=config
        )
        caltable.update()
        calfile = str(tmp / "cal_NGC7027.hdf5")
        caltable.write(calfile, overwrite=True)
        scanset = ScanSet(config)
        images = scanset.calculate_images(calibration=calfile, map_unit="Jy/pixel")
    return dict(
        config=config,
        tmp=tmp,
        caltable=caltable,
        scanset=scanset,
        images=images,
        files_before=files_before,
    )


def _aperture_flux(img):
    """Background-subtracted flux in a circle around the brightest pixel."""
    y0, x0 = np.unravel_index(np.nanargmax(img), img.shape)
    yy, xx = np.indices(img.shape)
    r = np.hypot(xx - x0, yy - y0) * PIXEL
    good = np.isfinite(img)
    ap = good & (r < R_AP)
    bkg = np.median(img[good & (r >= R_IN) & (r < R_OUT)])
    return np.sum(img[ap] - bkg)


def test_medicina_no_files_written_next_to_data(casa):
    """With a productdir, nothing (e.g. .hdf5 caches) is written in the data directory."""
    assert _list_data_files() == casa["files_before"]


def test_medicina_casa_flux(casa):
    """Cas A flux in RCP and LCP is roughly right, and the two agree within errors."""
    fluxes, errors = [], []
    for ch in ["Feed0_RCP", "Feed0_LCP"]:
        flux = _aperture_flux(np.array(casa["images"][ch], dtype=float))
        fc, fce = casa["caltable"].Jy_over_counts_rough(channel=ch, map_unit="Jy/pixel")
        fluxes.append(flux)
        errors.append(flux * fce / fc)
        # Loose: the calibrator data are poor (see technical_details)
        assert 0.5 * CASA_EXPECTED < flux < 1.3 * CASA_EXPECTED
    assert np.abs(fluxes[0] - fluxes[1]) < 2 * np.hypot(*errors)


def test_medicina_calibration_error_includes_scatter(casa):
    """The calibration error reflects the >10% scatter of the calibrator fits, in both channels."""
    for ch in ["Feed0_RCP", "Feed0_LCP"]:
        fc, fce = casa["caltable"].Jy_over_counts_rough(channel=ch, map_unit="Jy/pixel")
        assert fce / fc > 0.05


def test_medicina_command_line_output_directories(casa):
    """SDTcal and SDTimage put final products in productdir, caches and dumps in intermediatedir."""
    from srttools.calibration import main_cal
    from srttools.imager import main_imager

    tmp = casa["tmp"]
    products = tmp / "products"
    intermediate = products / "intermediate"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        main_cal(["-c", casa["config"]])
        main_imager(["-c", casa["config"], "--noplot"])
    for fname in ["casa_cal.hdf5", "casa_cal.csv", "casa.fits"]:
        assert (products / fname).exists(), fname
    for fname in ["casa_cal_unfilt.hdf5", "casa_dump.hdf5"]:
        assert (intermediate / fname).exists(), fname
    assert sorted(os.listdir(tmp)) == ["cal_NGC7027.hdf5", "casa.ini", "products"]
