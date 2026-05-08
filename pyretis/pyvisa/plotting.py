# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This file contains common functions for the visualization.

It contains some functions that are used to plot regression lines
and interface planes, and generate surface plots.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gen_surface (:py:func:`.gen_surface`)
    Generates a user-defined surface/contour/etc plot with colorbar in
    given matplotlib.figure and -.axes objects.

plot_int_plane(:py:func:`.plot_regline`)
    Generates interface planes for the current span of x-values, in a
    given matplotlib.axes-object.

plot_regline (:py:func:`.plot_regline`)
    Calculates the linear regression and correlation, plots a line for the
    regression in the given matplotlib.axes-object, with info in legend.

_grid_it_up (:py:func:`._grid_it_up`)
    Maps the x,y and z data to a numpy.meshgrid using scipy interpolation
    at a user defined resolution.
"""

import numpy as np
import matplotlib as mpl  # pylint: disable=import-error
from scipy.interpolate import (  # pylint: disable=import-error
    griddata as scgriddata,
)
from scipy.stats import linregress as linreg  # pylint: disable=import-error


def _add_colorbar(fig, mappable, ax, cbar_ax=None):
    """Add a colorbar with Matplotlib's current explicit axes API."""
    if cbar_ax is not None:
        return fig.colorbar(mappable, cax=cbar_ax)
    return fig.colorbar(mappable, ax=ax)


def plot_regline(ax, x, y):
    """Plot a regression line calculated from input data in the input subplot.

    Parameters
    ----------
    ax : Matplotlib subplot, where reg.line is to be plotted.
    x, y : list
        Floats, coordinates of data regression lines
        are calculated from.

    Updates
    -------
    Regression line with values.

    """
    xplot = np.linspace(min(x, default=0), max(x, default=0), 2)
    slope, intercept, r_value, _, _ = linreg(x, y)
    rtxt = 'y={0:.2f}x + {1:.2f}, $r^2$={2:.3f}'
    rline = ax.plot(xplot, slope*xplot + intercept, '-', c='black',
                    label=rtxt.format(slope, intercept, r_value**2))
    return rline


def plot_int_plane(ax, pos, yminmax, zminmax, visible=False):
    """Generate the interface planes for 3D visualization.

    Parameters
    ----------
    ax : The matplotlib.axes object where the planes will be plotted.
    pos : float
        The x-axis position of the interface plane.
    yminmax, zminmax : list of float
        The limits of the plane in the 3D canvas.
    visible : boolean, optional
        If True, shows interface planes.

    Returns
    -------
    plane : A 3D surface at x=pos, perpendicular to the x-axis.

    """
    y, z = np.linspace(*yminmax, 2), np.linspace(*zminmax, 2)
    y, z = np.meshgrid(y, z)
    point = np.array([pos, 0.0, 0.0])
    normal = np.array([1.0, 0.0, 0.0])
    dot = -point.dot(normal)
    x = (-normal[2]*y - normal[1]*z - dot) * 1./normal[0]
    plane = ax.plot_surface(x, y, z, color='grey',
                            alpha=0.30, visible=visible)
    return plane


def gen_surface(x, y, z, fig, ax, cbar_ax=None, dim=3, method='contour',
                res_x=400, res_y=400, colormap='viridis'):
    """Generate the chosen surface/contour/scatter plot.

    Parameters
    ----------
    x, y, z : list
        Coordinates of data points. (x,y) the chosen orderP pairs, and
        z is the chosen energy value of the two combinations.
    fig : Matplotlib object
        main canvas.
    ax : Matplotlib object
        axes of the plot.
    cbar_ax : Matplotlib object, optional
        plot color-bar.
    dim: int
        dimension of the plot surface.
    method : string, optional
        Method used for plotting data, default is contour lines.
    res_x, res_y : integer, optional
        Resolution of plot, either as N*N bins in 2D histogram
        (Density plot) or as grid-points for interpolation of data
        (Surface and contour plots).
    colormap : string, optional
        Name of the colormap/color scheme to use when plotting.

    Returns
    -------
    surf : Matplotlib object
         The chosen surface/contour/plot object.
    cbar : Matplotlib object
         The chosen color-bar.

    """
    # pylint: disable=too-many-arguments, too-many-positional-arguments
    # pylint: disable=too-many-locals
    xmin, xmax = min(x, default=0), max(x, default=0)
    ymin, ymax = min(y, default=0), max(y, default=0)
    zmin, zmax = min(z, default=0), max(z, default=0)

    cmap = mpl.colormaps[colormap]

    colors = [cmap(z[i]) for i in range(len(z))]
    if not zmin == zmax:
        colors = [cmap((z[i]-zmin)/(zmax-zmin)) for i in range(len(z))]

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    if dim == 3:
        # 3d plot settings
        ax.set_zlim3d(zmin, zmax)
        ax.zaxis.set_ticklabels([])

    # Methods for plotting in 3D
    if method == 'scatter':
        if dim == 3:
            surf = ax.scatter(x, y, z, c=colors, s=res_x / 100.0)
        else:
            surf = ax.scatter(x, y, c=colors, s=res_x / 100.0)

        norm = mpl.colors.Normalize(vmin=zmin, vmax=zmax)
        cbar = _add_colorbar(
            fig,
            mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
            ax,
            cbar_ax=cbar_ax)
        return surf, cbar

    g_x, g_y, g_z = _grid_it_up([x, y, z], res_x=res_x, res_y=res_y)
    if method == 'surface':
        surf = ax.plot_surface(g_x, g_y, g_z, cmap=cmap, vmin=zmin, vmax=zmax,
                               facecolor=colors, shade=True,
                               antialiased=False)

    elif method == 'contourf':
        surf = ax.contourf(g_x, g_y, g_z, cmap=cmap)

    else:  # contour
        surf = ax.contour(g_x, g_y, g_z, cmap=cmap)

    cbar = _add_colorbar(fig, surf, ax, cbar_ax=cbar_ax)

    return surf, cbar


def _grid_it_up(xyz, res_x=200, res_y=200, fill='max'):
    """Map x, y and z data values to a numpy meshgrid by interpolation.

    Parameters
    ----------
    xyz : list of list
        Lists of data values.
    res_x, res_y : integer, optional
        Resolution (number of points in a axis range).
    fill : string, optional
        Criteria to color the un-explored regions.

    Returns
    -------
    g_x, g_y, g_z : list
        Numpy.arrays of mapped data.

    """
    # Convert 3 columns of data to grid for matplotlib"""
    x_i = np.linspace(min(xyz[0], default=0), max(xyz[0], default=0), res_x)
    y_i = np.linspace(min(xyz[1], default=0), max(xyz[1], default=0), res_y)
    g_x, g_y = np.meshgrid(x_i, y_i)

    fill_value = 0
    # Scipy griddata """ # Works
    if fill == 'max':
        fill_value = max(xyz[2], default=0)
    elif fill == 'min':
        fill_value = min(xyz[2], default=0)
    g_z = scgriddata((xyz[0], xyz[1]), np.array(xyz[2]), (g_x, g_y),
                     method='linear', fill_value=fill_value)

    return g_x, g_y, g_z
