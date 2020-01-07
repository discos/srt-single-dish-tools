# Replace the default logging configuration with a custom one
from astropy import log
log.handlers.clear()
from astropy.logger import logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

MAX_FEEDS = 7
