# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Result comparison for GROMACS REPPTIS simulation.

Here we compare a RETIS simulation of 50 cycles to known results.
The initial path is produced via a load sparse and the simulation
tries to generate new valid paths.
"""
import os
import sys
import colorama
from pyretis.inout import print_to_screen
from pyretis.inout.settings import parse_settings_file
from pyretis.core.pathensemble import generate_ensemble_name
from pyretis.testing.simulation_comparison import (
    compare_path_ensemble_data,
    compare_numerical_data
)


RESULTS = 'results'


def print_underline(msg):
    """Print a message to the screen with a text-underline.

    Parameters
    ----------
    msg : string
        The message to print.
    """
    print_to_screen(f'\n{msg}', level='message')
    print_to_screen('=' * len(msg), level='message')


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
    files = ('pathensemble.txt', 'order.txt', 'energy.txt')
    for i in range(len(inter)):
        ensemble_dir = generate_ensemble_name(i)
        print_underline(f'Comparing for ensemble: {ensemble_dir}')
        for file_name in files:
            print_to_screen(f'* Comparing {file_name} files...')
            result_old = os.path.join(RESULTS, ensemble_dir, file_name)
            result_new = os.path.join(ensemble_dir, file_name)

            if file_name == 'pathensemble.txt':
                equal, msg = compare_path_ensemble_data(result_new, result_old)
            else:
                equal, msg = compare_numerical_data(result_new, result_old)

            if not equal:
                print_to_screen(f'\t-> *Files differ: {msg}*', level='error')
                return False
            print_to_screen('\t-> Files are equal!', level='success')
    print_to_screen('All files are equal!', level='success')
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
