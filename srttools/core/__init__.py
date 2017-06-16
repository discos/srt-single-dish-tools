import warnings
DEFAULT_MPL_BACKEND = 'TkAgg'
try:
    import matplotlib

    # This is necessary. Random backends might respond incorrectly.
    matplotlib.use(DEFAULT_MPL_BACKEND)
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

try:
    import statsmodels.api as sm
    version = [int(i) for i in sm.version.version.split('.')]

    # Minimum version 0.8.0
    if version < [0, 8, 0]:
        warnings.warn("Please update statsmodels")
        raise ImportError

    HAS_STATSM = True
except ImportError:
    HAS_STATSM = False

try:
    from numba import jit, vectorize
except ImportError:
    warnings.warn("Numba not installed. Faking it")

    def jit(fun):
        return fun

    def vectorize(*args, **kwargs):
        return jit

from .scan import Scan
from .imager import ScanSet
from .calibration import CalibratorTable
