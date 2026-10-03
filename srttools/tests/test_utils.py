import numpy as np
import pytest

from srttools.utils import HAS_MAHO, calculate_zernike_moments, look_for_files_or_bust


@pytest.mark.skipif("not HAS_MAHO")
def test_zernike_moments():
    image = np.ones((101, 101))
    res = calculate_zernike_moments(
        image, cm=[50, 50], radius=0.2, norder=8, label=None, use_log=False
    )
    assert res[1][1] < 1e-10
    assert res[3][1] < 1e-10


def test_look_for_files_or_bust():
    from pathlib import Path

    Path("blabla1.txt").touch()
    Path("blabla2.txt").touch()
    with pytest.raises(FileNotFoundError, match=".+blabla3.txt.+"):
        look_for_files_or_bust(["blabla1.txt", "blabla2.txt", "blabla3.txt"], timeout=1)


@pytest.mark.parametrize("display", [None, ":0"])
def test_matplotlib_backend_not_forced(display):
    """Importing srttools keeps the backend chosen by the user (MPLBACKEND), with or without a display."""
    import os
    import subprocess
    import sys

    import srttools

    pytest.importorskip("matplotlib")
    env = {k: v for k, v in os.environ.items() if k not in ["DISPLAY", "CI"]}
    env["MPLBACKEND"] = "svg"
    # Import this copy of srttools, not another installed one
    env["PYTHONPATH"] = os.path.dirname(os.path.dirname(srttools.__file__))
    if display is not None:
        env["DISPLAY"] = display
    code = "import srttools.utils, matplotlib; print(matplotlib.get_backend())"
    out = subprocess.run(
        [sys.executable, "-c", code], env=env, capture_output=True, text=True, check=True
    )
    assert out.stdout.strip().lower() == "svg"
