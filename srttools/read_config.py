"""Read the configuration file."""

import glob
import os
import shutil
import warnings

import numpy as np

import astropy.units as u

# For Python 2 and 3 compatibility
try:
    import configparser
except ImportError:
    import ConfigParser as configparser


__all__ = [
    "DEFAULT_INTERMEDIATE_SUBDIR",
    "DEFAULT_PRODUCTDIR",
    "get_config_file",
    "output_file_name",
    "read_config",
    "sample_config_file",
    "upgrade_config",
]


#: Directory (relative to the current one) where products go when the config
#: file does not specify a ``productdir``.
DEFAULT_PRODUCTDIR = "srttools_products"

#: Subdirectory of ``productdir`` where intermediate products go when the
#: config file does not specify an ``intermediatedir``.
DEFAULT_INTERMEDIATE_SUBDIR = "intermediate"

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
; the directory where the final products (images, calibration tables, light
; curves, cross scan results) are saved. If not specified, a directory called
; srttools_products is created in the directory where the scripts are launched.
; The raw data directory is never written to, unless productdir points to it.
; productdir : my_products
; the directory where the intermediate products (processed scans, diagnostic
; plots, the _dump.hdf5 files of SDTimage) are saved. Default: the subdirectory
; "intermediate" of productdir.
; intermediatedir : my_products/intermediate
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


def output_file_name(config_file, tail, outdir=None, intermediate=False):
    """Name of an output file of a tool run on ``config_file``, creating its directory.

    The name is ``<outdir>/<config file name without .ini><tail>``, e.g.
    ``products/casa_cal.hdf5`` for ``casa.ini`` and ``tail="_cal.hdf5"``.

    Parameters
    ----------
    config_file : str
        The config file. Only its name is used if ``outdir`` is given.
    tail : str
        What follows the root of the config file name.
    outdir : str, optional
        The directory. Default: ``productdir`` of the config file, for final
        products, or ``intermediatedir`` if ``intermediate`` is True.
    intermediate : bool
        The file is an intermediate product (e.g. a cache).
    """
    if outdir is None:
        outdir = read_config(config_file)["intermediatedir" if intermediate else "productdir"]
    os.makedirs(outdir, exist_ok=True)
    root = os.path.splitext(os.path.basename(config_file))[0]
    return os.path.join(outdir, root + tail)


def upgrade_config(fname):
    """Make an old config file keep using its processed scans.

    Older versions wrote the processed scans directly into ``productdir``,
    which now holds the final products, while the processed scans go to
    ``intermediatedir`` (default: ``<productdir>/intermediate``). This adds
    ``intermediatedir`` with the same value as ``productdir`` right after the
    ``productdir`` line, keeping the rest of the file (comments included)
    as it is. The original is saved as ``<fname>.bak``.

    Returns True if the file was changed. Nothing is done if the file
    already has an ``intermediatedir``, or has no ``productdir``.
    """
    with open(fname) as fobj:
        lines = fobj.readlines()

    def _key(line):
        line = line.strip()
        if not line or line[0] in ";#[":
            return None, None
        for sep in ":=":
            if sep in line:
                key, value = line.split(sep, 1)
                return key.strip().lower(), value.strip()
        return None, None

    keys = [_key(line) for line in lines]
    if any(k == "intermediatedir" for k, _ in keys):
        return False
    idx = [i for i, (k, v) in enumerate(keys) if k == "productdir" and v.lower() != "none"]
    if not idx:
        return False
    i = idx[0]
    value = keys[i][1]
    line = lines[i]
    indent = line[: len(line) - len(line.lstrip())]
    shutil.copyfile(fname, fname + ".bak")
    lines.insert(
        i + 1,
        f"{indent}; added by SDTupgradeconfig: processed scans of older versions are here\n"
        f"{indent}intermediatedir : {value}\n",
    )
    with open(fname, "w") as fobj:
        fobj.writelines(lines)
    return True


def main_upgrade_config(args=None):
    """Add ``intermediatedir`` to old config files (see :func:`upgrade_config`)."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Upgrade config files of older srttools versions: add "
        "intermediatedir, equal to productdir, so that the scans processed by "
        "older versions are used again. The originals are saved as <file>.bak"
    )
    parser.add_argument("files", nargs="+", help="Config files to upgrade")
    args = parser.parse_args(args)
    for fname in args.files:
        if upgrade_config(fname):
            print(f"{fname}: intermediatedir added (original in {fname}.bak)")
        else:
            print(f"{fname}: nothing to do")


def _warn_old_product_layout(config, fname):
    """Warn if ``productdir`` holds processed scans of an older version.

    Only when no ``intermediatedir`` was given and the default one does not
    exist yet; looks at most two directory levels deep.
    """
    productdir = config["productdir"]
    if (
        config["_explicit_intermediatedir"]
        or os.path.exists(config["intermediatedir"])
        or not os.path.isdir(productdir)
    ):
        return
    # Final products, named after the config file, are not processed scans
    final_root = os.path.splitext(os.path.basename(fname))[0]
    old = [
        f
        for f in glob.glob(os.path.join(productdir, "*.hdf5"))
        if not os.path.basename(f).startswith(final_root)
    ] + glob.glob(os.path.join(productdir, "*", "*.hdf5"))
    if old:
        warnings.warn(
            f"{productdir} contains processed scans of an older srttools version, "
            f"which are not used any more: processed scans now go to "
            f"{config['intermediatedir']}. To keep using the old ones, run "
            f"`SDTupgradeconfig {fname}`"
        )


def _set_default_intermediatedir(config):
    """Put ``intermediatedir`` under ``productdir``, unless it was chosen explicitly."""
    if not config["_explicit_intermediatedir"]:
        config["intermediatedir"] = os.path.join(config["productdir"], DEFAULT_INTERMEDIATE_SUBDIR)


def read_config(fname=None, productdir=None, intermediatedir=None):
    """Read a config file and return a dictionary of all entries.

    Parameters
    ----------
    fname : str, optional
        The config file. Default: the last one read or, if none was read
        yet, a sample one.
    productdir : str, optional
        Directory of the products, replacing the one of the config file
        (e.g. given on the command line). Relative paths start from the
        current directory. It stays valid for the later calls that read the
        same config file.
    intermediatedir : str, optional
        Directory of the intermediate products, replacing the one of the
        config file. Same rules as ``productdir``. Without it, and without
        ``intermediatedir`` in the config file, the intermediate products go
        to the subdirectory ``intermediate`` of ``productdir``.
    """
    global SRT_tools_config_file, SRT_tools_config

    # --- If already read, use existing config ---
    if SRT_tools_config is not None and (fname is None or fname == SRT_tools_config_file):
        if productdir is not None:
            SRT_tools_config["productdir"] = os.path.abspath(productdir)
        if intermediatedir is not None:
            SRT_tools_config["intermediatedir"] = os.path.abspath(intermediatedir)
            SRT_tools_config["_explicit_intermediatedir"] = True
        _set_default_intermediatedir(SRT_tools_config)
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
    config_output["intermediatedir"] = None
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

    if productdir is not None:
        config_output["productdir"] = os.path.abspath(productdir)

    if config_output["productdir"] is None:
        config_output["productdir"] = os.path.abspath(DEFAULT_PRODUCTDIR)
        warnings.warn(
            f"No productdir in the config file {fname}: the products (processed scans, "
            f"diagnostic plots) will be written to {config_output['productdir']}, "
            "NOT next to the raw data. Set productdir in the [local] section of the "
            "config file to choose another directory."
        )

    if config_output["productdir"] is not None and not config_output["productdir"].startswith("/"):
        config_output["productdir"] = os.path.abspath(
            os.path.join(os.path.split(fname)[0], config_output["productdir"])
        )

    if (
        config_output["intermediatedir"] is not None
        and config_output["intermediatedir"].lower() == "none"
    ):
        config_output["intermediatedir"] = None
    if intermediatedir is not None:
        config_output["intermediatedir"] = os.path.abspath(intermediatedir)
    elif config_output["intermediatedir"] is not None:
        config_output["intermediatedir"] = os.path.abspath(
            os.path.join(os.path.split(fname)[0], config_output["intermediatedir"])
        )
    config_output["_explicit_intermediatedir"] = config_output["intermediatedir"] is not None
    _set_default_intermediatedir(config_output)
    _warn_old_product_layout(config_output, fname)

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
