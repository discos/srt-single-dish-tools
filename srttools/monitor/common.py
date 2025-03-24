import os
import sys

# Replace the default logging configuration with a custom one
from srttools import logging

log = logging.getLogger("SDTmonitor")
log.propagate = False
sh = logging.StreamHandler()
f = logging.Formatter("%(asctime)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
sh.setFormatter(f)
log.addHandler(sh)
log.setLevel(logging.INFO)

MAX_FEEDS = 19
# The worst case scenario is fits0 through fits18, even though I'm not sure this notation will be kept in the future.
# We set a maximum cap for the number of concurrent processes to the number of CPU cores present in the system minus one.
# If there are more CPU cores than the maximum number of feeds, we set the number of maximum concurrent processes equal to the maximum number of feeds.
MAX_PROCS = min(MAX_FEEDS, os.cpu_count() - 1)

exit_function = os._exit

if "PYTEST_CURRENT_TEST" in os.environ:
    exit_function = sys.exit
