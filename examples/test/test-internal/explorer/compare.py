# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Compare explorer output with checked-in reference results."""
import argparse
import os
import sys
import colorama
from pyretis.inout.settings import parse_settings_file
from pyretis.core.pathensemble import generate_ensemble_name
from pyretis.setup.createsimulation import create_ensembles
from pyretis.testing.simulation_comparison import (
    compare_numerical_data,
    compare_path_ensemble_data,
)
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


RESULTS = 'results'
FILES = ('pathensemble.txt', 'order.txt', 'energy.txt')


def compare_files(settings_file='retis.rst'):
    """Compare generated explorer files with the stored reference data."""
    settings = parse_settings_file(settings_file)
    status = 0

    for ensemble in create_ensembles(settings):
        ens_number = ensemble['path_ensemble'].ensemble_number
        ens_dir = generate_ensemble_name(ens_number)
        logger.info('\nComparing ensemble: %s', ens_dir)
        for filename in FILES:
            result_new = os.path.join(ens_dir, filename)
            result_old = os.path.join(RESULTS, ens_dir, filename)
            if filename == 'pathensemble.txt':
                equal, msg = compare_path_ensemble_data(result_new, result_old)
            else:
                equal, msg = compare_numerical_data(result_new, result_old)
            if equal:
                logger.info('  %s: %s', filename, msg)
            else:
                logger.error('  %s: %s', filename, msg)
                status += 1
    return status


def main():
    """Run the comparison."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--settings',
        default='retis.rst',
        help='Settings file to use when determining the ensembles.',
    )
    return compare_files(parser.parse_args().settings)


if __name__ == '__main__':
    colorama.init(autoreset=True)
    sys.exit(main())
