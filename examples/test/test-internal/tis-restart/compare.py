# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Result comparison for TIS restart simulation.

Here we compare a full simulation with one where we have stopped
and restarted after 10 steps.
"""
import os
import sys
import colorama
import numpy as np
from pyretis.inout import print_to_screen
from pyretis.inout.settings import parse_settings_file
from pyretis.core.pathensemble import generate_ensemble_name
from pyretis.inout.formats.path import PathIntFile
from pyretis.testing.simulation_comparison import (
    compare_numerical_mse,
    compare_data_by_columns,
    compare_numerical_data
)


def snapshot_difference(snap1, snap2):
    """Calculate difference between two snapshots.

    Parameters
    ----------
    snap1 : dict
        The first snapshot.
    snap2 : dict
        The second snapshot.

    Returns
    -------
    diff_pos : float
        Sum of squared position differences.
    diff_vel : float
        Sum of squared velocity differences.
    """
    diff = (snap1['pos'] - snap2['pos'])**2
    dsum = np.einsum('ij,ij -> i', diff, diff)
    diffv = (snap1['vel'] - snap2['vel'])**2
    dsumv = np.einsum('ij,ij -> i', diffv, diffv)
    return sum(dsum), sum(dsumv)


def compare_traj(traj1, traj2, tol=1e-12):
    """Compare two trajectories from PyRETIS.

    Here we calculate the mean squared error for the two
    trajectories.

    Parameters
    ----------
    traj1 : string
        A trajectory file to open.
    traj2 : string
        A trajectory file to open.
    tol : float, optional
        A tolerance for comparing numbers.

    Returns
    -------
    status : int
        0 if comparison matches, 1 otherwise.
    """
    print_to_screen('Comparing trajectories', level='info')
    print_to_screen(f'Loading files: {traj1} & {traj2}')
    print_to_screen('Checking mean squared error...')
    file1 = PathIntFile(traj1, 'r').load()
    file2 = PathIntFile(traj2, 'r').load()
    error, error_v = [], []
    for trj1, trj2 in zip(file1, file2):
        for snap1, snap2 in zip(trj1['data'], trj2['data']):
            pose, vele = snapshot_difference(snap1, snap2)
            error.append(pose)
            error_v.append(vele)
    if next(file1, False) or next(file2, False):
        print_to_screen('Number of lines are incorrect ',
                        level='error')
        return 1
    ret1 = print_error_assessment(np.mean(error), 'positions', tol)
    ret2 = print_error_assessment(np.mean(error_v), 'velocities', tol)
    return ret1 + ret2


def print_error_assessment(error, what, tol):
    """Print out some error info.

    Parameters
    ----------
    error : float
        The calculated error.
    what : string
        Description of what was checked.
    tol : float
        The tolerance.

    Returns
    -------
    status : int
        0 if within tolerance, 1 otherwise.
    """
    if abs(error) < tol:
        lev = 'success'
        ret = 0
    else:
        lev = 'error'
        ret = 1
    print_to_screen(f'Mean error - {what}: {error}',
                    level=lev)
    return ret


def compare_ensemble(ensemble):
    """Run the comparison for an ensemble.

    Parameters
    ----------
    ensemble : string
        The ensemble name.

    Returns
    -------
    status : int
        Combined status of comparisons.
    """
    print_to_screen(f'Comparing for "{ensemble}"', level='info')
    traj1 = os.path.join('run-10-20', ensemble, 'traj.txt')
    traj2 = os.path.join('run-full', ensemble, 'traj.txt')
    print_to_screen()
    ret1 = compare_traj(traj1, traj2, tol=1e-12)

    retval = ret1
    for fname, ftype, msg in [
        ('energy.txt', 'energy', 'Comparing energies'),
        ('order.txt', None, 'Comparing order parameters'),
    ]:
        print_to_screen(f'\n{msg}', level='info')
        f1 = os.path.join('run-10-20', ensemble, fname)
        f2 = os.path.join('run-full', ensemble, fname)
        if ftype == 'energy':
            equal, res_msg = compare_data_by_columns(f1, f2, ftype)
        else:
            equal, res_msg = compare_numerical_data(f1, f2)
        if not equal:
            print_to_screen(f'\t-> *Files differ: {res_msg}*', level='error')
            retval += 1
        else:
            print_to_screen('\t-> Files are equal!', level='success')

    return retval


def main():
    """Run all comparisons."""
    settings = parse_settings_file(os.path.join('run-full', 'tis-001.rst'))
    ens = generate_ensemble_name(settings['tis'].get('ensemble_number'))
    ret = compare_ensemble(ens)
    if ret == 0:
        print_to_screen('\nAll seems fine!', level='success')
    else:
        print_to_screen('\nComparison failed!', level='error')
    return ret


if __name__ == '__main__':
    colorama.init(autoreset=True)
    sys.exit(main())
