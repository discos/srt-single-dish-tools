import pytest

try:
    import tornado
    import watchdog

    MONITOR_DEPENDENCIES = True
except ImportError:
    MONITOR_DEPENDENCIES = False


def pytest_ignore_collect(collection_path):
    if MONITOR_DEPENDENCIES:
        return False
    return True
