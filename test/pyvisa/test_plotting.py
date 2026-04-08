# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the common methods in pyretis.pyvisa.plotting."""
import os
import logging
import pytest
import numpy as np
from matplotlib import pyplot
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D  # It is needed
from pyretis.pyvisa.plotting import (
    plot_regline,
    plot_int_plane,
    gen_surface,
    _grid_it_up)

logging.disable(logging.CRITICAL)

HERE = os.path.abspath(os.path.dirname(__file__))


class TestMethods:
    """Test some of the methods from pyretis.pyvisa.plotting."""

    def test_plot_regline(self):
        """Test of the regression plotting"""
        _, ax = pyplot.subplots()
        x = [i for i in range(10)]
        y = [i**2 - 1. for i in x]
        # Expected regression line points
        ptm = [x[0]*9. - 13., x[-1]*9. - 13.]
        line = Line2D([x[0], x[-1]], ptm)
        # Regression line
        regline = plot_regline(ax, x, y)
        # Loop over y-values and compare
        for i, j in zip(regline[0].get_ydata(), line.get_ydata()):
            assert i == j

    def test_plot_int_plane(self):
        """Test if we can make an interface plane in 3d plot"""
        fig = pyplot.figure()
        ax = fig.add_subplot(111, projection='3d')
        # Interface plane specifics
        pos = 1.0
        ymin, ymax, zmin, zmax = 0.0, 1.0, 0.0, 1.0
        visible = False
        _ = plot_int_plane(ax, pos, (ymin, ymax), (zmin, zmax),
                           visible=visible)
        # Nothing to check, only if any errors should arise from plotting
        fig.clear()

    def test_gen_surface(self):
        """Test if we can generate a surface/plot with given settings"""
        # Create data
        x, y, z = [], [], []
        for i in range(-10, 10):
            for j in range(-10, 10):
                x.append(i/10.)
                y.append(j/10.)
                z.append((i/10.)**2 + (j/10.)**2)
        # Create figure and plot+colorbar axes
        fig = pyplot.figure()
        ax = fig.add_subplot(111, projection='3d')
        cbar_ax = fig.add_axes([0.86, 0.1, 0.03, 0.8])
        # Generate surface 3d
        _ = gen_surface(x, y, z, fig, ax, cbar_ax=cbar_ax,
                        res_x=200, res_y=200, method='surface')
        ax.clear()
        _ = gen_surface(x, y, z, fig, ax, cbar_ax=cbar_ax,
                        res_x=200, res_y=200, method='contour')
        ax.clear()
        _ = gen_surface(x, y, z, fig, ax, cbar_ax=cbar_ax,
                        res_x=200, res_y=200, method='contourf')
        ax.clear()
        _ = gen_surface(x, y, z, fig, ax, cbar_ax=cbar_ax,
                        method='scatter')
        ax.clear()

        fig.clear()
        # Generate surface 2d
        ax = fig.add_subplot(111)
        cbar_ax = fig.add_axes([0.1, 0.03, 0.8, 0.1])
        _ = gen_surface(x, y, z, fig, ax, cbar_ax=cbar_ax, dim=2,
                        method='contourf')
        ax.clear()
        _ = gen_surface(x, y, z, fig, ax, cbar_ax=cbar_ax, dim=2,
                        method='contour')
        ax.clear()
        _ = gen_surface(x, y, len(x)*[0.0], fig, ax, cbar_ax=cbar_ax, dim=2,
                        method='scatter')
        ax.clear()

    def test_grid_it_up(self):
        """Test grid method"""
        x = np.array([1, 1, 1, 1, 1, 2, 9])
        y = np.array([1, 2, 2, 4, 4, 1, 9])
        z = np.array([1, 2, 7, 4, 5, 1, 1])

        mins = _grid_it_up([x, y, z], res_x=20, res_y=20, fill='max')
        maxs = _grid_it_up([x, y, z], res_x=20, res_y=20, fill='min')
        assert not np.allclose(mins, maxs)
