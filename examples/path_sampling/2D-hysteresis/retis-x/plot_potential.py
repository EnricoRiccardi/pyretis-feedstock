# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Plot raw data from a simulation."""
# pylint: disable=invalid-name
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.cm import get_cmap
from pyretis.core import create_box, System, Particles
from pyretis.setup import create_force_field
from pyretis.inout.settings import parse_settings_file


def plot_potential(settings, axi, axj):
    """Plot the potential in the given axis."""
    forcefield = create_force_field(settings)
    box = create_box(periodic=[False, False])
    fakesys = System(units='reduced', box=box)
    fakesys.particles = Particles(dim=2)
    fakesys.add_particle(name='B', pos=np.zeros(2), ptype=1)
    minx, maxx = -0.5, 0.5
    miny, maxy = -1.0, 1.0
    xval = np.linspace(minx, maxx, 100)
    yval = np.linspace(miny, maxy, 100)
    xpos, ypos = np.meshgrid(xval, yval, indexing='ij')
    pot = np.zeros_like(xpos)
    for i, x in enumerate(xval):
        for j, y in enumerate(yval):
            fakesys.particles.pos[0, 0] = x
            fakesys.particles.pos[0, 1] = y
            pot[i, j] = forcefield.evaluate_potential(fakesys)
    axi.contourf(xpos, ypos, pot, 10, cmap=get_cmap('viridis'))
    cont2 = axj.contour(xpos, ypos, pot, 10, cmap=get_cmap('viridis'))
    cont2.clabel(inline=1, fmt='%3.1f', fontsize=12)
    xmin = settings['potential'][0]['parameter']['x0']
    ymin = settings['potential'][0]['parameter']['y0']
    for axx in (axi, axj):
        axx.set_xlim((minx, maxx))
        axx.set_ylim((miny, maxy))
        axx.set_xlabel((r'Position ($x$)'), fontsize='large')
        axx.set_ylabel((r'Position ($y$)'), fontsize='large')
    axi.scatter(xmin, ymin, s=50, marker='o', color='white')
    axi.scatter(-xmin, -ymin, s=50, marker='o', color='white')
    axj.scatter(xmin, ymin, s=50, marker='o', color='#262626')
    axj.scatter(-xmin, -ymin, s=50, marker='o', color='#262626')


if __name__ == '__main__':
    sim_settings = parse_settings_file('retis.rst')
    fig = plt.figure(figsize=(12, 6))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    plot_potential(sim_settings, ax1, ax2)
    fig.tight_layout()
    plt.show()
