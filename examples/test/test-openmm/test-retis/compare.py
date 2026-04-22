# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Result comparison for OpenMM RETIS simulation.

Here we compare a RETIS simulation to known results stored in the
``openmm/results/`` directory.
"""
import os
import sys
import colorama
from pyretis.inout.settings import parse_settings_file
from pyretis.core.pathensemble import generate_ensemble_name
from pyretis.testing.simulation_comparison import (
    compare_path_ensemble_data,
    compare_numerical_data,
)
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

RUNDIR = 'openmm'
RESULTS = os.path.join(RUNDIR, 'results')


def print_underline(msg):
    logger.info(f'\n{msg}')
    logger.info('=' * len(msg))


def compare_files(settings):
    """Compare output files for all ensembles against reference results.

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
    files = ('pathensemble.txt', 'order.txt', 'energy.txt')

    for i in range(len(inter)):
        ensemble_dir = generate_ensemble_name(i)
        print_underline(f'Comparing for ensemble: {ensemble_dir}')
        for file_name in files:
            logger.info(f'* Comparing {file_name} ...')
            result_ref = os.path.join(RESULTS, ensemble_dir, file_name)
            result_new = os.path.join(RUNDIR, ensemble_dir, file_name)

            if not os.path.exists(result_ref):
                logger.error(f'\t-> *Reference missing: {result_ref}*')
                return False

            if file_name == 'pathensemble.txt':
                equal, msg = compare_path_ensemble_data(result_new, result_ref)
            else:
                equal, msg = compare_numerical_data(result_new, result_ref)

            if not equal:
                logger.error(f'\t-> *Files differ: {msg}*')
                return False
            logger.info('\t-> Files are equal!')

    logger.info('All files are equal!')
    return True


def main():
    """Run the comparison."""
    settings = parse_settings_file(os.path.join(RUNDIR, 'retis.rst'))
    print_underline('Comparing OpenMM RETIS results to reference.')
    if not compare_files(settings):
        sys.exit(1)


if __name__ == '__main__':
    colorama.init(autoreset=True)
    logging.basicConfig(level=logging.INFO)
    main()
