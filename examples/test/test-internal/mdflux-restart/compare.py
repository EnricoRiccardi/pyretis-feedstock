# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Result comparison for MDFlux restart simulation.

Here we compare a full simulation with one where we have stopped
and restarted after 100 steps.
"""
import os
import sys
import colorama
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
from pyretis.testing.simulation_comparison import (
    compare_restarted_text_files,
    compare_restarted_cross_files
)
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

plt.style.use('seaborn-v0_8')


def make_fig():
    """Plot for comparison."""
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
    ax3.set_ylabel('Order parameter.')
    ax3.plot([], [], label=r'$\lambda$', lw=0, alpha=0)
    ax3.plot([], [], label=r'$\dot{\lambda}}$', lw=0, alpha=0)
    axes = (ax1, ax2, ax3)
    return fig1, axes


def plot_in_ax(axes, infile, lab, fat=False, colors=None, style='-'):
    """Just do some plotting."""
    ax1, ax2, _ = axes
    data = np.loadtxt(infile)
    width = 7 if fat else 3
    lines = []
    for i, idx in enumerate((1, 2, 3)):
        color = colors[i] if colors else None
        line, = ax1.plot(data[:, 0], data[:, idx], label=lab,
                         ls=style, lw=width, alpha=0.8, color=color)
        lines.append(line)
    ax2.plot(data[:, 0], data[:, 4], label=lab, ls=style, lw=width, alpha=0.9)
    return lines


def plot_in_ax_op(axes, infile, lab, fat=False, colors=None, style='-'):
    """Just do some plotting for order parameter."""
    _, _, ax3 = axes
    data = np.loadtxt(infile)
    width = 7 if fat else 3
    lines = []
    for i, idx in enumerate((1, 2)):
        color = colors[i] if colors else None
        line, = ax3.plot(data[:, 0], data[:, idx], label=lab,
                         ls=style, lw=width, alpha=0.8, color=color)
        lines.append(line)
    return lines


def make_plots():
    """Just plot some energies for comparison."""
    figure, axes = make_fig()

    plot_in_ax(
        axes,
        os.path.join('run-full', 'energy.txt'),
        'full',
        fat=True,
        style='-'
    )
    lines = plot_in_ax(
        axes,
        os.path.join('run-step1', 'energy.txt'),
        'restart-part1',
        style='--'
    )
    colors = [i.get_color() for i in lines]
    plot_in_ax(
        axes,
        os.path.join('run-step2', 'energy.txt'),
        'restart-part2',
        style=':', colors=colors
    )

    plot_in_ax_op(
        axes,
        os.path.join('run-full', 'order.txt'),
        'full',
        fat=True,
        style='-'
    )
    lines = plot_in_ax_op(
        axes,
        os.path.join('run-step1', 'order.txt'),
        'restart-part1',
        style='--'
    )
    colors = [i.get_color() for i in lines]
    plot_in_ax_op(
        axes,
        os.path.join('run-step2', 'order.txt'),
        'restart-part2',
        style=':', colors=colors
    )

    axes[0].legend(prop={'size': 'medium'}, ncol=4)
    axes[1].legend(prop={'size': 'medium'})
    axes[2].legend(prop={'size': 'medium'}, ncol=4)
    figure.subplots_adjust(bottom=0.12, right=0.95, top=0.95,
                           left=0.08, wspace=0.2)
    return figure


def main(args):
    """Run comparisons.

    Parameters
    ----------
    args : list of str
        Arguments passed to the script.

    Returns
    -------
    result : int
        The number of comparison failures.
    """
    result = 0
    logger.info('Comparing crossings...')
    equal, msg = compare_restarted_cross_files(
        os.path.join('run-step1', 'cross.txt'),
        os.path.join('run-step2', 'cross.txt'),
        os.path.join('run-full', 'cross.txt'),
    )
    if not equal:
        logger.error(f'Crossing files mismatch: {msg}')
        result += 1
    else:
        logger.info('Crossing files match!')

    pairs = [
        ('energy.txt', 'Comparing energy files'),
        ('order.txt', 'Comparing order files'),
    ]

    for fname, msg_header in pairs:
        logger.info(f'\n{msg_header}:')
        equal, msg = compare_restarted_text_files(
            os.path.join('run-step1', fname),
            os.path.join('run-step2', fname),
            os.path.join('run-full', fname)
        )
        if not equal:
            logger.error(f'Restarted {fname} mismatch: {msg}')
            result += 1
        else:
            logger.info(f'Restarted {fname} match!')

    if 'make_plot' in args:
        fig = make_plots()
        fig.savefig('compare.png')
        plt.show()

    if result != 0:
        logger.error('\nComparison failed!')
    else:
        logger.info('\nComparison was successful!')
    return result


if __name__ == '__main__':
    colorama.init(autoreset=True)
    sys.exit(main(sys.argv[1:]))
