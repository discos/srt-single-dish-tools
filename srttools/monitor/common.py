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

MAX_FEEDS = 7

exit_function = os._exit

testing_environments = ["TRAVIS", "CI", "CONTINUOUS_INTEGRATION"]

for env in testing_environments:
    if env in os.environ:
        exit_function = sys.exit
        break
