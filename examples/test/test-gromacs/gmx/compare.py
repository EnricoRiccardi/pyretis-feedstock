# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Simple script to compare the outcome of two simulations."""
import argparse
import os
import sys
import colorama
from pyretis.inout import print_to_screen
from pyretis.inout.settings import parse_settings_file
from pyretis.core.pathensemble import generate_ensemble_name
from pyretis.testing.simulation_comparison import (
    compare_simulation_files,
    compare_data_by_columns,
    compare_path_ensemble_data,
    compare_traj_archive,
)


TRAJ = 'traj'
RETIS_RST = 'retis.rst'


def soft_archive_comparison(dir1, dir2):
    """Compare two archives. Here we accept a difference in energies.

    For GROMACS1 and 2, the potential energy output is different.
    This is due to a difference in the way GROMACS adds dispersion
    corrections when continuing a simulation. So, here we expect
    a failure for the energy files. We check these energy files
    manually, by skipping the potential energy terms.

    Parameters
    ----------
    dir1 : string
        The path to the first archive used in the comparison.
    dir2 : string
        The path to the second archive used in the comparison.

    Returns
    -------
    errors : list of tuples
        These are the files which were found to be different.
    """
    errors = []
    archive_diffs = compare_traj_archive(dir1, dir2)
    for (file1, file2) in archive_diffs:
        # Check if it is an energy file:
        if os.path.basename(file1) == 'energy.txt':
            equal, _ = compare_data_by_columns(
                file1, file2, file_type='energy', skip=['vpot']
            )
            if not equal:
                errors.append((file1, file2))
        else:
            # For other files in archive, they must match exactly:
            equal, _ = compare_simulation_files(file1, file2, mode='cmp')
            if not equal:
                errors.append((file1, file2))
    return errors


def compare_ensemble(run1, run2, ensemble,
                     path_skip=None, energy_skip=None, traj_skip=False):
    """Run the comparison for an ensemble.

    Parameters
    ----------
    run1 : string
        The path to the first directory to use for the comparison.
    run2 : string
        The path to the second directory to use for the comparison.
    ensemble : string
        The label for the current ensemble we are investigating, e.g.
        000, 001, or similar.
    path_skip : list of integers, optional
        This selects columns to skip when comparing path ensemble
        data.
    energy_skip : list of strings, optional
        This selects columns to skip in the energy file.
    traj_skip : bool, optional
        If True, we skip comparison of trajectory archives.

    Returns
    -------
    status : bool
        True if the comparison was successful, False otherwise.
    """
    print_to_screen(f'Comparing for "{ensemble}"', level='info')
    for filei in ('energy.txt', 'order.txt', 'pathensemble.txt'):
        print_to_screen(f'\tComparing {filei} files...')
        file1 = os.path.join(run1, ensemble, filei)
        file2 = os.path.join(run2, ensemble, filei)
        if filei == 'pathensemble.txt':
            equal, msg = compare_path_ensemble_data(
                file1, file2, skip=path_skip
            )
        elif filei == 'energy.txt' and energy_skip:
            equal, msg = compare_data_by_columns(
                file1, file2, file_type='energy', skip=energy_skip
            )
        else:
            equal, msg = compare_simulation_files(
                file1, file2, mode='numerical'
            )
        lvl = 'success' if equal else 'error'
        print_to_screen(f'\t\t-> {msg}', level=lvl)
        if not equal:
            return False

    if not traj_skip:
        print_to_screen('\tComparing for trajectory archive:')
        # Check folders like 0_traj and 1_traj
        for i in {'0_', '1_'}:
            archive_errors = soft_archive_comparison(
                os.path.join(run1, ensemble, i + TRAJ),
                os.path.join(run2, ensemble, i + TRAJ),
            )
            if archive_errors:
                print_to_screen('\t\t-> Archives differ', level='error')
                return False
        print_to_screen('\t\t-> Archives are equal', level='success')
    else:
        print_to_screen('\tSkipping comparison for trajectory archive.')

    return True


def main():
    """Parse arguments and execute the comparison.

    Returns
    -------
    out : int
        0 if the comparison was successful, 1 otherwise.
    """
    colorama.init(autoreset=True)
    parser = argparse.ArgumentParser()
    parser.add_argument('directories', metavar='N', type=str, nargs=2,
                        help='The directories to use for the comparison')
    parser.add_argument('--path_skip', nargs='+', type=int,
                        help='Columns to skip for pathensemble.txt')
    parser.add_argument('--energy_skip', nargs='+', type=str,
                        help='Columns to skip for energy.txt')
    parser.add_argument('--traj_skip', action='store_true',
                        help='Skip comparison for trajectory archives')
    args = parser.parse_args()

    run1 = args.directories[0]
    run2 = args.directories[1]

    settings = parse_settings_file(os.path.join(run1, RETIS_RST))
    inter = settings['simulation']['interfaces']
    for i, _ in enumerate(inter):
        ensemble_dir = generate_ensemble_name(i)
        result = compare_ensemble(
            run1, run2, ensemble_dir,
            path_skip=args.path_skip,
            energy_skip=args.energy_skip,
            traj_skip=args.traj_skip,
        )
        if not result:
            print_to_screen(
                f'Comparison failed for {ensemble_dir}. Aborting!',
                level='error',
            )
            return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
