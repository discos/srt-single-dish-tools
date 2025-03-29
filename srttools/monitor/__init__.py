from srttools.monitor.common import main_monitor
import importlib

__all__ = ["Monitor", "main_monitor"]


# We try to load Monitor dynamically.
# This will allow the CLI to work and print the help regardless of the installed dependencies.
def __getattr__(name):
    if name == "Monitor":
        mod = importlib.import_module(".monitor", __package__)
        monitor = getattr(mod, "Monitor")
        globals()["Monitor"] = monitor
        return monitor
    return object.__getattribute__(globals(), name)


def __dir__():  # pragma: no cover
    # This block only serves IDEs/IPython autocompletion feature
    return __all__
