"""
Random utilities
"""

def standard_string(s):
    if hasattr(s, 'decode'):
        s = s.decode()
    return s


def standard_byte(s):
    if hasattr(s, 'encode'):
        s = s.encode()
    return s
