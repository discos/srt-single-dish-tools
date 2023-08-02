# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *  # noqa

# ----------------------------------------------------------------------------
try:
    import faulthandler

    faulthandler.enable()
except ImportError:
    pass

import warnings
import numpy as np


def warning_on_one_line(message, category, filename, lineno, line=None):
    return "%s:%s: %s: %s\n" % (filename, lineno, category.__name__, message)


warnings.formatwarning = warning_on_one_line

# warnings.filterwarnings("error", category=np.VisibleDeprecationWarning)
# warnings.filterwarnings("error", ".*")
warnings.filterwarnings("once", category=UserWarning)
warnings.filterwarnings("once", category=DeprecationWarning)
warnings.filterwarnings("once", category=np.VisibleDeprecationWarning)
warnings.filterwarnings("ignore", "table path was not set via the path= ")

__all__ = []
from .logging import logging, logger

logger.setLevel(logging.INFO)
