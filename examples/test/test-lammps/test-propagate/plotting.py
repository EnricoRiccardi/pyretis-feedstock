# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Common plotting methods used by the LAMMPS engine tests."""
from math import ceil
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def plot_compare(data_sets, ylabels):
    """Plot data for comparison.

    Parameters
    ----------
    data_sets : list of list of tuples
        The list are of the form `[(x, y, label), (x, y, label)]`.
    ylabels : list of strings
        The y-labels to use for each data set.

    """
    markers = ['o', 'x', '1']
    fig = plt.figure()
    if len(data_sets) < 3:
        ncol = 1
    else:
        ncol = 2
    grid = make_grid(ncol, len(data_sets))
    for i, data in enumerate(data_sets):
        row, col = divmod(i, ncol)
        axi = fig.add_subplot(grid[row, col])
        for j, dataj in enumerate(data):
            line, = axi.plot(dataj[0], dataj[1], ls='--', lw=1, alpha=0.5)
            axi.scatter(dataj[0], dataj[1], marker=markers[j],
                        s=100, alpha=0.7, label=dataj[2],
                        color=line.get_color())
        axi.set_xlabel('Step')
        axi.set_ylabel(ylabels[i])
        axi.legend()
    fig.tight_layout()


def make_grid(ncol, nsets):
    """Create a grid for plotting some data sets.

    Parameters
    ----------
    ncol : int
        The number of columns we want to plot in.
    nsets : int
        The number of data sets to plot.

    """
    nrow = ceil(nsets / ncol)
    grid = GridSpec(nrow, ncol)
    return grid


def plot_xy(data_sets):
    """Plot x-data vs y-data.

    This is intended to plot two data sets we expect to be equal,
    in order to visualize differences between the two data sets.

    Parameters
    ----------
    data_sets : list of tuples
        The tuples are on form `(x, y, xlabel, ylabel)`.

    """
    fig = plt.figure()
    ncol = 2
    grid = make_grid(ncol, len(data_sets))
    for i, data in enumerate(data_sets):
        row, col = divmod(i, ncol)
        axi = fig.add_subplot(grid[row, col])
        axi.scatter(data[0], data[1], marker='o', s=90, alpha=0.5)
        axi.set_xlabel(data[2])
        axi.set_ylabel(data[3])
        minxy = min(min(data[0]), min(data[1]))
        maxxy = max(max(data[0]), max(data[1]))
        axi.plot([minxy, maxxy], [minxy, maxxy], ls=':',
                 color='#262626', alpha=0.5, lw=2)
        axi.set_aspect('equal')
        axi.set_xlim([minxy, maxxy])
        axi.set_ylim([minxy, maxxy])
        if not np.allclose(data[0], data[1]):
            logger.error(f'Comparison: "{data[2]}" & "{data[3]}" failed.')
        else:
            logger.info(f'Comparison: "{data[2]}" & "{data[3]}" is OK.')
    fig.tight_layout()
