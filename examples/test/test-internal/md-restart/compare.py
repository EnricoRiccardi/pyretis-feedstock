# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Result comparison for MD restart simulation.

Here we compare a full simulation with one where we have stopped
and restarted after 10 steps.
"""
import os
import sys
import itertools
import colorama
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
from pyretis.inout.formats.snapshot import SnapshotFile
from pyretis.testing.simulation_comparison import (
    compare_restarted_text_files
)
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

plt.style.use('seaborn-v0_8')


def snapshot_difference(snap1, snap2):
    """Calculate difference between two snapshots.

    Parameters
    ----------
    snap1 : dict
        A dictionary containing snapshot data (x, y, z, vx, vy, vz).
    snap2 : dict
        Identical to snap1, for comparison.

    Returns
    -------
    pos_diff : float
        Sum of squared differences in positions.
    vel_diff : float
        Sum of squared differences in velocities.
    """
    xyz1 = np.column_stack((snap1['x'], snap1['y'], snap1['z']))
    xyz2 = np.column_stack((snap2['x'], snap2['y'], snap2['z']))
    diff = (xyz1 - xyz2)**2
    dsum = np.einsum('ij,ij -> i', diff, diff)
    vel1 = np.column_stack((snap1['vx'], snap1['vy'], snap1['vz']))
    vel2 = np.column_stack((snap2['vx'], snap2['vy'], snap2['vz']))
    diffv = (vel1 - vel2)**2
    dsumv = np.einsum('ij,ij -> i', diffv, diffv)
    return sum(dsum), sum(dsumv)


def compare_traj(traj11, traj12, traj2, tol=1e-12):
    """Compare two trajectories from PyRETIS.

    Here we calculate the mean squared error for the two
    trajectories.

    Parameters
    ----------
    traj11 : string
        A trajectory to open, part 1.
    traj12 : string
        A trajectory to open, part 2.
    traj2 : string
        A trajectory file to open.
    tol : float
        A tolerance for comparing numbers.

    Returns
    -------
    status : int
        0 if comparison matches, 1 otherwise.
    """
    logger.info('Comparing trajectories')
    logger.info('Checking mean squared error...')
    file11 = SnapshotFile(traj11, 'r').load()
    file12 = SnapshotFile(traj12, 'r').load()
    # Skip the first configuration in the restart file
    # (it's the last configuration of part 1).
    next(file12, None)
    file1 = itertools.chain(file11, file12)
    file2 = SnapshotFile(traj2, 'r').load()
    error, error_v = [], []
    for snap1, snap2 in zip(file1, file2):
        pose, vele = snapshot_difference(snap1, snap2)
        error.append(pose)
        error_v.append(vele)
    if next(file1, False) or next(file2, False):
        logger.error('Number of lines are incorrect ')
        return 1
    val1 = print_error_assessment(np.mean(error), 'positions', tol)
    val2 = print_error_assessment(np.mean(error_v), 'velocities', tol)
    return val1 + val2


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
        val = 0
    else:
        lev = 'error'
        val = 1
    if lev == 'error':
        logger.error(f'Mean error - {what}: {error}')
    else:
        logger.info(f'Mean error - {what}: {error}')
    return val


def make_plots():
    """Just plot some energies for comparison."""
    fig1 = plt.figure(figsize=(12, 6))
    grid = gridspec.GridSpec(2, 2)
    ax1 = fig1.add_subplot(grid[:, 0])
    ax1.plot([], [], label='Potential', lw=0, alpha=0)
    ax1.plot([], [], label='Kinetic', lw=0, alpha=0)
    ax1.plot([], [], label='Total', lw=0, alpha=0)
    ax1.set_xlabel('Step no.')
    ax1.set_ylabel('Energy per particle')
    ax2 = fig1.add_subplot(grid[0, 1])
    ax2.set_ylabel('Temperature')
    ax3 = fig1.add_subplot(grid[1, 1])
    ax3.set_xlabel('Step no.')
    ax3.set_ylabel('Pressure')
    axes = (ax1, ax2, ax3)

    plot_in_ax(
        axes,
        os.path.join('run-full', 'md-full-thermo.txt'),
        'full',
        fat=True,
        style='-'
    )
    lines = plot_in_ax(
        axes,
        os.path.join('run-10', 'md-10-thermo.txt'),
        'restart-part1',
        style='--'
    )
    colors = [i.get_color() for i in lines]
    plot_in_ax(
        axes,
        os.path.join('run-10-100', 'md-10-100-thermo.txt'),
        'restart-part2',
        style=':', colors=colors)
    axes[0].legend(prop={'size': 'medium'}, ncol=4)
    axes[1].legend(prop={'size': 'medium'})
    axes[2].legend(prop={'size': 'medium'})
    fig1.subplots_adjust(bottom=0.12, right=0.95, top=0.95,
                         left=0.08, wspace=0.2)
    return fig1


def plot_in_ax(axes, infile, lab, fat=False, colors=None, style='-'):
    """Just do some plotting.

    Parameters
    ----------
    axes : tuple
        The axes to plot in.
    infile : string
        The file to load data from.
    lab : string
        The label for the plot.
    fat : bool, optional
        True for thicker lines.
    colors : list, optional
        Colors for the lines.
    style : string, optional
        Line style.

    Returns
    -------
    lines : list
        The plotted lines.
    """
    ax1, ax2, ax3 = axes
    data = np.loadtxt(infile)
    width = 7 if fat else 3
    lines = []
    for i, idx in enumerate((2, 3, 4)):
        color = colors[i] if colors else None
        line, = ax1.plot(data[:, 0], data[:, idx], label=lab,
                         ls=style, lw=width, alpha=0.8, color=color)
        lines.append(line)
    ax2.plot(data[:, 0], data[:, 1], label=lab, ls=style, lw=width, alpha=0.9)
    ax3.plot(data[:, 0], data[:, 5], label=lab, ls=style, lw=width, alpha=0.9)
    return lines


def main(args):
    """Run the comparison.

    Parameters
    ----------
    args : list
        Command line arguments.
    """
    val1 = compare_traj(
        os.path.join('run-10', 'md-10-traj.xyz'),
        os.path.join('run-10-100', 'md-10-100-traj.xyz'),
        os.path.join('run-full', 'md-full-traj.xyz'),
        tol=1e-12
    )

    pairs = [
        ('md-10-thermo.txt', 'md-10-100-thermo.txt', 'md-full-thermo.txt'),
        ('md-10-energy.txt', 'md-10-100-energy.txt', 'md-full-energy.txt'),
    ]

    val2 = 0
    for p1, p2, pf in pairs:
        equal, msg = compare_restarted_text_files(
            os.path.join('run-10', p1),
            os.path.join('run-10-100', p2),
            os.path.join('run-full', pf)
        )
        if not equal:
            logger.error(f'Restarted files mismatch ({p1}): {msg}')
            val2 += 1
        else:
            logger.info(f'Restarted files match ({pf})')

    if 'make_plot' in args:
        fig = make_plots()
        fig.savefig('compare.png')
        plt.show()
    return val1 + val2


if __name__ == '__main__':
    colorama.init(autoreset=True)
    sys.exit(main(sys.argv[1:]))
