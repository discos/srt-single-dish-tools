"""Read the configuration file."""

import glob
import os
import warnings

import numpy as np

import astropy.units as u

# For Python 2 and 3 compatibility
try:
    import configparser
except ImportError:
    import ConfigParser as configparser


__all__ = ["get_config_file", "read_config", "sample_config_file"]


SRT_tools_config_file = None
SRT_tools_config = None


def sample_config_file(fname="sample_config_file.ini"):
    """Create a sample config file, to be modified by hand."""
    string = """
[local]
; the directory where the analysis will be executed.
workdir : .
; the root directory of the data repository.
datadir : .
; the root directory of the data repository.
productdir : None
; if True, output files will be saved directly inside productdir the full path.
basename_only : False

[analysis]
projection : ARC
interpolation : spline
prefix : test_
list_of_directories :
;;Two options: either a list of directories:
;        dir1
;        dir2
;; or a star symbol for all directories
;         *
calibrator_directories :
; if left empty, calibrator scans are taken from list_of_directories when
; calculating light curves, and ignored when calculating images

skydip_directories :
; if left empty, calibrator scans are taken from list_of_directories when
; calculating light curves, and ignored when calculating images

noise_threshold : 5

;; For spectral rms smoothing, in percentage of then number of spectral bins.
smooth_window : 0.05

;; Coordinates have to be specified in decimal degrees. ONLY use if different
;; from target coordinates!
;    reference_ra : 10.5
;    reference_dec : 5.3

;; Pixel size in arcminutes. If set to 'auto', the pixel size will be calculated
;; as approximately the size of the beam divided by 3.

pixel_size : auto

;; Channels to save from RFI filtering. It might indicate known strong spectral
;; lines
goodchans :

[debugging]

debug_file_format : jpg

    """
    with open(fname, "w") as fobj:
        print(string, file=fobj)
    return fname


def get_config_file():
    """Get the current config file."""
    return SRT_tools_config_file


def read_config(fname=None):
    """Read a config file and return a dictionary of all entries."""
    global SRT_tools_config_file, SRT_tools_config

    # --- If already read, use existing config ---
    if fname == SRT_tools_config_file and SRT_tools_config is not None:
        return SRT_tools_config

    if fname is None and SRT_tools_config is not None:
        return SRT_tools_config

    # ---------------------------------------------

    config_output = {}

    Config = configparser.ConfigParser()

    if fname is None:
        fname = sample_config_file()
    else:
        if not os.path.exists(fname):
            raise FileNotFoundError("Please specify an existing config file")

    SRT_tools_config_file = fname
    Config.read(fname)

    # ---------- Set default values --------------------------------------

    config_output["projection"] = "ARC"
    config_output["interpolation"] = "linear"
    config_output["workdir"] = os.path.abspath(os.curdir) + "/"
    config_output["datadir"] = os.path.abspath(os.curdir) + "/"
    config_output["productdir"] = None
    config_output["basename_only"] = False
    config_output["list_of_directories"] = "*"
    config_output["calibrator_directories"] = []
    config_output["skydip_directories"] = []
    config_output["pixel_size"] = "auto"
    config_output["goodchans"] = None
    config_output["filtering_factor"] = "0"
    config_output["noise_threshold"] = "5"
    config_output["smooth_window"] = "0.05"
    config_output["debug_file_format"] = "jpg"
    config_output["ignore_suffix"] = []
    config_output["ignore_prefix"] = []

    # --------------------------------------------------------------------

    # Read local information

    local_params = dict(Config.items("local"))

    config_output.update(local_params)
    if not config_output["workdir"].startswith("/"):
        config_output["workdir"] = os.path.abspath(
            os.path.join(os.path.split(fname)[0], config_output["workdir"])
        )

    if not config_output["datadir"].startswith("/"):
        config_output["datadir"] = os.path.abspath(
            os.path.join(os.path.split(fname)[0], config_output["datadir"])
        )

    if config_output["productdir"] is not None and config_output["productdir"].lower() == "none":
        config_output["productdir"] = None

    if config_output["productdir"] is not None and not config_output["productdir"].startswith("/"):
        config_output["productdir"] = os.path.abspath(
            os.path.join(os.path.split(fname)[0], config_output["productdir"])
        )

    try:
        # Read analysis information
        debugging_params = dict(Config.items("debugging"))

        config_output.update(debugging_params)
    except configparser.NoSectionError:
        pass

    # Read analysis information
    analysis_params = dict(Config.items("analysis"))

    config_output.update(analysis_params)

    try:
        config_output["list_of_directories"] = [
            s for s in analysis_params["list_of_directories"].splitlines() if s.strip()
        ]  # This last instruction eliminates blank lines
    except Exception:
        warnings.warn("Invalid list_of_directories in config file")

    try:
        config_output["calibrator_directories"] = [
            s for s in analysis_params["calibrator_directories"].splitlines() if s.strip()
        ]  # This last instruction eliminates blank lines
    except Exception:
        warnings.warn("Invalid calibrator_directories in config file")

    try:
        config_output["skydip_directories"] = [
            s for s in analysis_params["skydip_directories"].splitlines() if s.strip()
        ]  # This last instruction eliminates blank lines
    except Exception:
        warnings.warn("Invalid skydip_directories in config file")

    for key, val in config_output.items():
        if isinstance(val, str) and "\n" in val:
            config_output[key] = [s for s in val.splitlines() if s.strip()]

    # If the list of directories is not specified, or if a '*' symbol is used,
    # use glob in the datadir to determine the list

    if config_output["list_of_directories"] in ([], ["*"], "*"):
        config_output["list_of_directories"] = [
            os.path.split(f)[1]  # return name without path
            for f in glob.glob(os.path.join(config_output["datadir"], "*"))
            if os.path.isdir(f)
        ]  # only if it's a directory

    if config_output["pixel_size"] != "auto":
        config_output["pixel_size"] = (
            float(np.radians(float(config_output["pixel_size"]) / 60)) * u.rad
        )

    if config_output["goodchans"] is not None:
        config_output["goodchans"] = [int(n) for n in config_output["goodchans"]]

    config_output["noise_threshold"] = float(config_output["noise_threshold"])
    config_output["smooth_window"] = float(config_output["smooth_window"])
    config_output["filtering_factor"] = float(config_output["filtering_factor"])

    SRT_tools_config = config_output
    return config_output
