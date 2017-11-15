"""
Random utilities
"""
import sys
import numpy as np
import warnings
import logging


try:
    from mahotas.features import zernike_moments
    HAS_MAHO = True
except ImportError:
    HAS_MAHO = False

try:
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

DEFAULT_MPL_BACKEND = 'TKAgg'


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


def _generic_dummy_decorator(*args, **kwargs):
    if len(args) == 1 and len(kwargs) == 0 and callable(args[0]):
        return args[0]
    else:
        def decorator(func):
            def decorated(*args, **kwargs):
                return func(*args, **kwargs)

            return decorated

        return decorator


try:
    from numba import jit, vectorize
except ImportError:
    warnings.warn("Numba not installed. Faking it")

    jit = vectorize = _generic_dummy_decorator


__all__ = ["mad", "standard_string", "standard_byte", "compare_strings",
           "tqdm", "jit", "vectorize"]


try:
    from tqdm import tqdm
except ImportError:
    def tqdm(x):
        return x


try:
    from statsmodels.robust import mad as mad  # pylint: disable=unused-import
except ImportError:
    def mad(data, c=0.6745, axis=None):
        """Straight from statsmodels's source code, adapted"""
        data = np.asarray(data)
        if axis is not None:
            center = np.apply_over_axes(np.median, data, axis)
        else:
            center = np.median(data)
        return np.median((np.fabs(data - center)) / c, axis=axis)


def standard_string(s):
    """Standard string representation for a given Python version

    Examples
    --------
    >>> standard_string(b'a')
    'a'
    >>> standard_string(None) is None
    True
    """
    if s is None:
        return None

    if sys.version_info >= (3, 0, 0):
        # for Python 3
        # This indexing should work for both lists of strings, and strings
        if hasattr(s, 'decode'):
            s = s.decode()  # or  s = str(s)[2:-1]
        # Try to see if it's a numpy array
        elif hasattr(s, 'dtype') and s.dtype.char == 'S':
            if s.size > 1:
                s = np.array(s, dtype='U')
    else:
        # for Python 2
        if isinstance(s[0], unicode):  # NOQA
            s = str(s)
        # Try to see if it's a numpy array
        elif hasattr(s, 'dtype') and s.dtype.char == 'U':
            if s.size > 1:
                s = np.array(s, dtype='S')
    return s


def standard_byte(s):
    """Standard byte representation for a given Python version

    Examples
    --------
    >>> standard_byte(b'a') == b'a'
    True
    >>> standard_byte(np.array([u'a'])[0]) == b'a'
    True
    """
    if hasattr(s, 'encode'):
        s = s.encode()
    elif hasattr(s, 'dtype') and s.dtype.char == 'U':
        if s.size > 1:
            s = np.array(s, dtype='S')
    return s


def compare_strings(s1, s2):
    """Compare strings, that might be bytes and unicode in some cases.

    Parameters
    ----------
    s1: string, byte or array of str/bytes
    s2 : string or byte

    Examples
    --------
    >>> compare_strings(b'a', 'a')
    True
    >>> compare_strings('a', u'a')
    True
    >>> import numpy as np
    >>> compare_strings(np.array(['a', 'b'], dtype='S'), u'a')
    array([ True, False], dtype=bool)
    """

    s1 = standard_string(s1)
    s2 = standard_string(s2)
    return s1 == s2


def interpolate_invalid_points_image(array):
    '''Interpolates invalid points in an image.

    Examples
    --------
    >>> img = np.ones((3, 3))
    >>> img[1, 1] = np.nan
    >>> np.all(interpolate_invalid_points_image(img) == np.ones((3, 3)))
    True
    '''
    from scipy import interpolate
    x = np.arange(0, array.shape[1])
    y = np.arange(0, array.shape[0])
    # mask invalid values
    array = np.ma.masked_invalid(array)
    xx, yy = np.meshgrid(x, y)
    # get only the valid values
    x1 = xx[~array.mask]
    y1 = yy[~array.mask]
    newarr = array[~array.mask]

    GD1 = interpolate.griddata((x1, y1), newarr.ravel(),
                               (xx, yy),
                               method='cubic', fill_value=0)
    return GD1


def calculate_zernike_moments(im, cm=None, radius=0.3, norder=8,
                              label=None, use_log=False):
    """Calculate the Zernike moments of the image.

    These moments are useful to single out asymmetries in the image:
    for example, when characterizing the beam of the radio telescope using
    a map of a calibrator, it is useful to calculate these moments to
    understand if the beam is radially symmetric or has distorted side
    lobes.

    Parameters
    ----------
    im : 2-d array
        The image to be analyzed

    Other parameters
    ----------------
    cm : [int, int]
        'Center of mass' of the image
    radius : float
        The radius around the center of mass, in percentage of the image
        size (0 <= radius <= 0.5)
    norder : int
        Maximum order of the moments to calculate
    use_log: bool
        Rescale the image to a log scale before calculating the coefficients.
        The scale is the same documented in the ds9 docs, for consistency.
        After normalizing the image from 0 to 1, the log-rescaled image is
        log(ax + 1) / log a, with ``x`` the normalized image and ``a`` a
        constant fixed here at 1000

    Returns
    -------
    moments_dict : dict
        Dictionary containing the order, the sub-index and the moment, e.g.
        {0: {0: 0.3}, 1: {1: 1e-16}, 2: {0: 0.95, 2: 6e-19}, ...}
        Moments are symmetrical, so only the unique values are reported.

    """
    if cm is None:
        cm = np.unravel_index(im.argmax(), im.shape)

    im_to_analyze = im.copy()
    if use_log:
        vmin = im_to_analyze.min()
        vmax = im_to_analyze.max()
        rescaled_image = (im_to_analyze - vmin) / (vmax - vmin)

        im_to_analyze = np.log(1000 * rescaled_image + 1) / np.log(1000)

    im_to_analyze = interpolate_invalid_points_image(im_to_analyze)

    radius_pix = np.int(np.min(im.shape) * radius)
    moments = zernike_moments(im_to_analyze, radius_pix, norder, cm=cm)
    count = 0
    moments_dict = {}
    description_string = \
        'Zernike moments (cm: {}, radius: {}):\n'.format(cm, radius_pix)

    if HAS_MPL:
        fig = plt.figure('Zernike moments', figsize=(10, 10))
        plt.imshow(im_to_analyze, vmin=0, vmax=im_to_analyze[cm[0], cm[1]],
                   origin='lower')
        circle = plt.Circle(cm, radius_pix, color='r', fill=False)
        plt.gca().add_patch(circle)
        plt.colorbar()

    for i in range(norder + 1):
        description_string += str(i) + ': '
        moments_dict[i] = {}
        for j in range(i + 1):
            if (i - j)%2 == 0:
                description_string += "{}/{} {:.1e} ".format(i, j,
                                                             moments[count])
                moments_dict[i][j] = moments[count]
                count += 1
        description_string += '\n'

    if HAS_MPL:
        plt.text(0.05, 0.95, description_string,
                 horizontalalignment='left',
                 verticalalignment = 'top',
                 transform = plt.gca().transAxes,
                 color='white')

        if label is None:
            label = str(np.random.randint(0, 100000))
        plt.savefig('Zernike_debug_' + label +
                    '.png')
        plt.close(fig)

    logging.debug(description_string)

    moments_dict['Description'] = description_string

    return moments_dict
