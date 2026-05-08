# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Result comparison for partial path simulation.

Here we compare a RETIS simulation of 50 cycles to known results.
The initial path is produced via a load sparse and the simulation
tries to generate new valid paths.
"""
import os
import sys
import colorama
from pyretis.inout.settings import parse_settings_file
from pyretis.core.pathensemble import generate_ensemble_name
from pyretis.testing.simulation_comparison import (
    compare_path_ensemble_data,
    compare_numerical_data
)
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


RESULTS = 'results'


def print_underline(msg):
    """Print a message to the screen with a text-underline.

    Parameters
    ----------
    msg : string
        The message to print.
    """
    logger.info(f'\n{msg}')
    logger.info('=' * len(msg))


def compare_files(settings):
    """Compare output files for all ensembles.

    Parameters
    ----------
    settings : dict
        The simulation settings.

    Returns
    -------
    status : bool
        True if all files are equal, False otherwise.
    """
    inter = settings['simulation']['interfaces']
    files = ('pathensemble.txt', 'order.txt')
    for i in range(len(inter)):
        ensemble_dir = generate_ensemble_name(i)
        print_underline(f'Comparing for ensemble: {ensemble_dir}')
        for file_name in files:
            logger.info(f'* Comparing {file_name} files...')
            result_old = os.path.join(RESULTS, ensemble_dir, file_name)
            result_new = os.path.join(ensemble_dir, file_name)

            if file_name == 'pathensemble.txt':
                equal, msg = compare_path_ensemble_data(result_new, result_old)
            else:
                equal, msg = compare_numerical_data(result_new, result_old)

            if not equal:
                logger.error(f'\t-> *Files differ: {msg}*')
                return False
            logger.info('\t-> Files are equal!')
    logger.info('All files are equal!')
    return True


def main():
    """Run the comparison."""
    settings = parse_settings_file('repptis.rst')
    print_underline('Comparing files for the load example.')
    compare_ok = compare_files(settings)
    if not compare_ok:
        sys.exit(1)


if __name__ == '__main__':
    colorama.init(autoreset=True)
    main()
