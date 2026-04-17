# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Simple script to compare the outcome of two simulations.

The outcome of the ``md_nve.py`` simulation should be independent (to numerical
precision) of the units used. This script will test that by comparing:

- the output in `thermo.txt`

- the generated trajectory, `traj.xyz`

For the energies, it will create a plot comparing the energies, the pressure
and the temperature.
"""
import sys
import colorama
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
from pyretis.inout.formats.snapshot import SnapshotFile
from pyretis.inout.settings import parse_settings_file
from pyretis.core.units import (
    create_conversion_factors,
    generate_system_conversions,
    CONVERT
)
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
plt.style.use('seaborn-colorblind')


def snapshot_difference(snap1, snap2, unit1, unit2):
    """Calculate difference between two snapshots."""
    xyz1 = np.column_stack((snap1['x'], snap1['y'], snap1['z']))
    xyz2 = np.column_stack((snap2['x'], snap2['y'], snap2['z']))
    xyz2 *= CONVERT['length'][unit1, unit2]
    diff = (xyz1 - xyz2)**2
    dsum = np.einsum('ij,ij -> i', diff, diff)
    vel1 = np.column_stack((snap1['vx'], snap1['vy'], snap1['vz']))
    vel2 = np.column_stack((snap2['vx'], snap2['vy'], snap2['vz']))
    vel2 *= CONVERT['velocity'][unit1, unit2]
    diffv = (vel1 - vel2)**2
    dsumv = np.einsum('ij,ij -> i', diffv, diffv)
    return sum(dsum), sum(dsumv)


def compare_traj(traj1, traj2, unit1, unit2, tol=1e-12):
    """Compare two trajectories from PyRETIS.

    Here we calculate the mean squared error for the two
    trajectories.

    Parameters
    ----------
    traj1 : string
        A trajectory file to open.
    traj2 : string
        A trajectory file to open.
    unit1 : string
        The system of units for the first trajectory.
    unit2 : string
        The system of units for the second trajectory.
    tol : float
        A tolerance for comparing numbers.

    Returns
    -------
    None, just prints out the result of the comparison.

    """
    retval = 1
    logger.info('Comparing trajectories')
    print('Checking mean squared error...')
    file1 = SnapshotFile(traj1, 'r').load()
    file2 = SnapshotFile(traj2, 'r').load()
    error, error_v = 0.0, 0.0
    nsnap = 0
    for snap1, snap2 in zip(file1, file2):
        pose, vele = snapshot_difference(snap1, snap2, unit1, unit2)
        error += pose
        error_v += vele
        nsnap += 1
    error /= float(nsnap)
    if abs(error) < tol:
        lev = 'success'
        retval = 0
    else:
        lev = 'error'
        return 1
    if lev == 'error':
        logger.error(f'Mean error - positions: {error}')
    else:
        logger.info(f'Mean error - positions: {error}')
    if abs(error_v) < tol:
        lev = 'success'
        retval = 0
    else:
        lev = 'error'
        return 1
    if lev == 'error':
        logger.error(f'Mean error - velocities: {error_v}')
    else:
        logger.info(f'Mean error - velocities: {error_v}')
    return retval


def run_comparison():
    """Run the comparison."""
    settings = parse_settings_file('settings.rst')
    unit = settings['system']['units']
    create_conversion_factors('lj')
    create_conversion_factors(unit)
    generate_system_conversions('lj', unit)
    ret = compare_traj('../traj.xyz', 'traj.xyz', unit, 'lj')
    return unit, 'lj', ret


def make_plot(unit1, unit2):
    """Plot for comparison."""
    ljunits = np.loadtxt('../thermo.txt')
    other_units = np.loadtxt('thermo.txt')
    # convert other_units:
    other_units[:, 1] *= CONVERT['temperature'][unit1, unit2]
    other_units[:, 2:5] *= CONVERT['energy'][unit1, unit2]
    other_units[:, 5] *= CONVERT['pressure'][unit1, unit2]
    fig1 = plt.figure(figsize=(12, 6))
    grid = gridspec.GridSpec(2, 2)
    ax1 = fig1.add_subplot(grid[:, 0])
    ljlab = f'"{unit2}"'
    unilab = f'"{unit1}"'
    ax1.plot([], [], label='Potential', lw=0, alpha=0)
    ax1.plot([], [], label='Kinetic', lw=0, alpha=0)
    ax1.plot([], [], label='Total', lw=0, alpha=0)
    ax1.plot(ljunits[:, 0], ljunits[:, 2], label=ljlab,
             ls='-', lw=7, alpha=0.8)
    ax1.plot(ljunits[:, 0], ljunits[:, 3], label=ljlab,
             ls='-', lw=7, alpha=0.8)
    ax1.plot(ljunits[:, 0], ljunits[:, 4], label=ljlab,
             ls='-', lw=7, alpha=0.8)
    ax1.plot(other_units[:, 0], other_units[:, 2], label=unilab,
             ls='--', lw=3, alpha=0.8)
    ax1.plot(other_units[:, 0], other_units[:, 3], label=unilab,
             ls='--', lw=3, alpha=0.8)
    ax1.plot(other_units[:, 0], other_units[:, 4], label=unilab,
             ls='--', lw=3, alpha=0.8)
    ax1.set_xlabel('Step no.')
    ax1.set_ylabel('Energy per particle')
    _ = ax1.legend(loc='center left', prop={'size': 'large'}, ncol=3)

    ax2 = fig1.add_subplot(grid[0, 1])
    ax2.plot(ljunits[:, 0], ljunits[:, 1], label=ljlab,
             ls='-', lw=7, alpha=0.9)
    ax2.plot(other_units[:, 0], other_units[:, 1], label=unilab,
             ls='--', lw=3, alpha=0.8, color='0.2')
    ax2.set_ylabel('Temperature')
    ax2.legend(loc='upper right', prop={'size': 'large'})
    ax3 = fig1.add_subplot(grid[1, 1])
    ax3.plot(ljunits[:, 0], ljunits[:, 5], label=ljlab,
             ls='-', lw=7, alpha=0.9)
    ax3.plot(other_units[:, 0], other_units[:, 5], label=unilab,
             ls='--', lw=3, alpha=0.8, color='0.2')
    ax3.set_xlabel('Step no.')
    ax3.set_ylabel('Pressure')
    ax3.legend(loc='lower right', prop={'size': 'large'})
    fig1.subplots_adjust(bottom=0.12, right=0.95, top=0.95,
                         left=0.08, wspace=0.2)
    fig1.savefig('compare.png')


def main(args):
    """Run the comparison."""
    uni1, uni2, ret = run_comparison()
    if 'plot' in args:
        make_plot(uni1, uni2)
        plt.show()
    return ret


if __name__ == '__main__':
    colorama.init(autoreset=True)
    sys.exit(main(sys.argv[1:]))
