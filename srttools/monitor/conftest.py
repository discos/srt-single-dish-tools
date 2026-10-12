from importlib.util import find_spec

MONITOR_DEPENDENCIES = find_spec("tornado") is not None and find_spec("watchdog") is not None


def pytest_ignore_collect(collection_path):
    if MONITOR_DEPENDENCIES:
        return False
    return True
