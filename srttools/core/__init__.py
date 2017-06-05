import warnings
try:
    import matplotlib

    # matplotlib.use('TkAgg')
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

try:
    import statsmodels.api as sm
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

