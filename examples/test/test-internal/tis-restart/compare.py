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
from pyretis.inout.settings import parse_settings_file
from pyretis.core.pathensemble import generate_ensemble_name
from pyretis.inout.formats.path import PathIntFile
from pyretis.testing.simulation_comparison import (
    compare_numerical_mse,
    compare_data_by_columns,
    compare_numerical_data
)
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


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
    logger.info('Comparing trajectories')
    logger.info(f'Loading files: {traj1} & {traj2}')
    logger.info('Checking mean squared error...')
    file1 = PathIntFile(traj1, 'r').load()
    file2 = PathIntFile(traj2, 'r').load()
    error, error_v = [], []
    for trj1, trj2 in zip(file1, file2):
        for snap1, snap2 in zip(trj1['data'], trj2['data']):
            pose, vele = snapshot_difference(snap1, snap2)
            error.append(pose)
            error_v.append(vele)
    if next(file1, False) or next(file2, False):
        logger.error('Number of lines are incorrect ')
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
    if lev == 'error':
        logger.error(f'Mean error - {what}: {error}')
    else:
        logger.info(f'Mean error - {what}: {error}')
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
    logger.info(f'Comparing for "{ensemble}"')
    traj1 = os.path.join('run-6-12', ensemble, 'traj.txt')
    traj2 = os.path.join('run-full', ensemble, 'traj.txt')
    ret1 = compare_traj(traj1, traj2, tol=1e-12)

    retval = ret1
    for fname, ftype, msg in [
        ('energy.txt', 'energy', 'Comparing energies'),
        ('order.txt', None, 'Comparing order parameters'),
    ]:
        logger.info(f'\n{msg}')
        f1 = os.path.join('run-6-12', ensemble, fname)
        f2 = os.path.join('run-full', ensemble, fname)
        if ftype == 'energy':
            equal, res_msg = compare_data_by_columns(f1, f2, ftype)
        else:
            equal, res_msg = compare_numerical_data(f1, f2)
        if not equal:
            logger.error(f'\t-> *Files differ: {res_msg}*')
            retval += 1
        else:
            logger.info('\t-> Files are equal!')

    return retval


def main():
    """Run all comparisons."""
    settings = parse_settings_file(os.path.join('run-full', 'tis-001.rst'))
    ens = generate_ensemble_name(settings['tis'].get('ensemble_number'))
    ret = compare_ensemble(ens)
    if ret == 0:
        logger.info('\nAll seems fine!')
    else:
        logger.error('\nComparison failed!')
    return ret


if __name__ == '__main__':
    colorama.init(autoreset=True)
    sys.exit(main())
