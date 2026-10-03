"""Configure Test Suite.

This file is used to configure the behavior of pytest when using the Astropy
test infrastructure. It needs to live inside the package in order for it to
get picked up when running the tests inside an interpreter using
packagename.test

"""

import os

import pytest

try:
    from pytest_astropy_header.display import (
        PYTEST_HEADER_MODULES,
        TESTED_VERSIONS,
    )

    ASTROPY_HEADER = True
except ImportError:
    ASTROPY_HEADER = False


def pytest_configure(config):
    """Configure Pytest with Astropy.

    Parameters
    ----------
    config : pytest configuration

    """
    if ASTROPY_HEADER:
        config.option.astropy_header = True

        # Customize the following lines to add/remove entries from the list of
        # packages for which version numbers are displayed when running the tests.
        PYTEST_HEADER_MODULES.pop("Pandas", None)
        PYTEST_HEADER_MODULES["scikit-image"] = "skimage"

        from . import __version__

        packagename = os.path.basename(os.path.dirname(__file__))
        TESTED_VERSIONS[packagename] = __version__


@pytest.fixture
def no_pyplot_figures(monkeypatch):
    """Make the test fail if the code creates a pyplot figure.

    Plots that are only saved to a file must not use pyplot: with an
    interactive backend, pyplot figures can open windows.
    """
    plt = pytest.importorskip("matplotlib.pyplot")

    def fail(*args, **kwargs):
        raise AssertionError("A pyplot figure was created for a plot that is only saved")

    monkeypatch.setattr(plt, "new_figure_manager", fail)
