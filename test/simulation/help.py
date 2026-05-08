# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Define common methods and variables for the tests"""
from contextlib import contextmanager
from pyretis.inout.formats.formatter import OutputFormatter
import logging


@contextmanager
def turn_on_logging():
    """Turn on logging so that tests can detect it."""
    logging.disable(logging.NOTSET)
    try:
        yield
    finally:
        logging.disable(logging.CRITICAL)


TEST_SETTINGS = {
    'system': {
        'dimensions': 3,
        'units': 'reduced',
        'temperature': 1.0
    },
    'particles': {
        'position': {
            'generate': 'fcc',
            'repeat': [2, 2, 2],
            'density': 0.9,
        }
    },
    'potential': [
        {
            'class': 'PairLennardJonesCutnp',
            'shift': True,
            'parameter': {
                0: {'sigma': 1.0, 'epsilon': 1.0, 'rcut': 2.5}
            }
        }
    ]
}


class TxtWriter(OutputFormatter):
    """A class used for testing output."""
    FMT = '{:>10d} {:>20s}'

    def format(self, step, data):
        yield self.FMT.format(step, data)
