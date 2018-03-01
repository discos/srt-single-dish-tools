# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This is an Astropy affiliated package.
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *  # NOQA
# ----------------------------------------------------------------------------
try:
    import faulthandler

    faulthandler.enable()
except ImportError:
    pass

# For egg_info test builds to pass, put package imports here.
if not _ASTROPY_SETUP_:
    from .calibration import *  # noqa: F401,F403
    from .fit import *  # noqa: F401,F403
    from .global_fit import *  # noqa: F401,F403
    from .histograms import *  # noqa: F401,F403
    from .imager import *  # noqa: F401,F403
    from .inspect_observations import *  # noqa: F401,F403
    from .interactive_filter import *  # noqa: F401,F403
    from .io import *  # noqa: F401,F403
    from .read_config import *  # noqa: F401,F403
    from .scan import *  # noqa: F401,F403
    from .simulate import *  # noqa: F401,F403
    from .utils import *  # noqa: F401,F403
