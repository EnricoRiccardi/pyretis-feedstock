# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Define common methods and variables for the tests."""
from contextlib import contextmanager
import logging


@contextmanager
def turn_on_logging():
    """Turn on logging so that tests can detect it."""
    logging.disable(logging.NOTSET)
    try:
        yield
    finally:
        logging.disable(logging.CRITICAL)
