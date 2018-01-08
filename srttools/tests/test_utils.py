import pytest
import numpy as np
from ..utils import HAS_MAHO, calculate_zernike_moments


@pytest.mark.skipif('not HAS_MAHO')
def test_zernike_moments():
    image = np.ones((101, 101))
    res = calculate_zernike_moments(image, cm=[50, 50], radius=0.2,
                                    norder=8, label=None, use_log=False)
    assert res[1][1] < 1e-10
    assert res[3][1] < 1e-10
