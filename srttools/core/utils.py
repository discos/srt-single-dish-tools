"""
Random utilities
"""
import sys
import numpy as np


def standard_string(s):
    """Standard string representation for a given Python version
    
    Examples
    --------
    >>> standard_string(b'a')
    'a'
    """
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
        if isinstance(s[0], unicode):
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
