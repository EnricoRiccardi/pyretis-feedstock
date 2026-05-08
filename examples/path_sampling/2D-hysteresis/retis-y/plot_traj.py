# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Plot raw data from a simulation."""
# pylint: disable=invalid-name
import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.cm import get_cmap
from pyretis.core import create_box, System, Particles
from pyretis.setup import create_force_field, create_orderparameter
from pyretis.inout.settings import parse_settings_file
from pyretis.inout.formats.path import PathIntFile


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
    axi.contourf(xpos, ypos, pot, 10, cmap=get_cmap('viridis'), alpha=0.8)
    # add interfaces
    for inter in settings['simulation']['interfaces']:
        axi.axhline(y=inter, lw=2, ls=':', color='#262626', alpha=0.8)
        axj.axhline(y=inter, lw=2, ls=':', color='#262626', alpha=0.8)
    extra_int = [settings['orderparameter']['inter_a'],
                 settings['orderparameter']['inter_b']]
    for inter in extra_int:
        axi.axhline(y=inter, lw=2, ls=':', color='#262626', alpha=0.5)
        axj.axhline(y=inter, lw=2, ls=':', color='#262626', alpha=0.5)
    axi.set_xlim((minx, maxx))
    axi.set_ylim((miny, maxy))
    axi.set_xlabel((r'Position ($x$)'), fontsize='large')
    axi.set_ylabel((r'Position ($y$)'), fontsize='large')
    axj.set_xlabel((r'Step number'), fontsize='large')
    axj.set_ylabel((r'Order parameter ($\lambda$)'), fontsize='large')


def plot_ensemble(settings, dirname, axi, axj, maxlines=100, minorder=None,
                  skip=1):
    """Plot trajectories from an ensemble."""
    orderp = create_orderparameter(settings)
    box = create_box(periodic=[False, False])
    fakesys = System(units='reduced', box=box)
    fakesys.particles = Particles(dim=2)
    fakesys.add_particle(name='B', pos=np.zeros(2), ptype=1)
    forcefield = create_force_field(settings)
    traj_file = os.path.join(dirname, 'traj.txt')
    iplot = 0
    all_lines = []
    all_lines2 = []
    last_point = []
    first_point = []
    order_last = []
    order_last2 = []
    with PathIntFile(traj_file, 'r') as tfile:
        for i, traj in enumerate(tfile.load()):
            if traj['comment'][0].split('status:')[-1].strip() != 'ACC':
                continue
            if i % skip != 0:
                continue
            pos = np.array([x['pos'][0] for x in traj['data']])
            order = []
            for posi in pos:
                fakesys.particles.pos[0, 0] = posi[0]
                fakesys.particles.pos[0, 1] = posi[1]
                fakesys.particles.vpot = forcefield.evaluate_potential(fakesys)
                order.append(orderp.calculate(fakesys)[0])
            if minorder is not None:
                if max(order) < minorder:
                    continue
            line, = axi.plot(pos[:, 0], pos[:, 1], lw=3, alpha=0.9)
            line2, = axj.plot(order, lw=3, alpha=0.9)
            order_last.append((len(order) - 1, order[-1]))
            order_last2.append((len(order) - 2, order[-2]))
            all_lines.append(line)
            all_lines2.append(line2)
            first_point.append((pos[0, 0], pos[0, 1]))
            last_point.append((pos[-1, 0], pos[-1, 1]))
            iplot += 1
            if iplot >= maxlines:
                break
    # Add colors now that we know how many we have created:
    cmap = get_cmap(name='coolwarm')
    colors = cmap(np.linspace(0, 1, iplot))
    for i, (line, line2) in enumerate(zip(all_lines, all_lines2)):
        line.set_color(colors[i])
        axi.scatter(first_point[i][0], first_point[i][1], s=50, marker='x',
                    color=line.get_color(), alpha=0.9)
        axi.scatter(last_point[i][0], last_point[i][1], s=50, marker='o',
                    color=line.get_color(), alpha=0.9)
        line2.set_color(colors[i])
        if order_last2[i][1] < order_last[i][1]:
            end = '^'
        else:
            end = 'v'
        axj.scatter(order_last[i][0], order_last[i][1], s=50, marker=end,
                    color=line2.get_color(), alpha=0.9)
        axj.scatter(order_last2[i][0], order_last2[i][1], s=50, marker='x',
                    color=line2.get_color(), alpha=0.9)


if __name__ == '__main__':
    ens = sys.argv[1]
    sim_settings = parse_settings_file('retis.rst')
    fig = plt.figure(figsize=(12, 6))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    plot_potential(sim_settings, ax1, ax2)
    plot_ensemble(sim_settings, ens, ax1, ax2, maxlines=25, skip=2,
                  minorder=None)
    plt.subplots_adjust(right=0.95, left=0.10, top=0.95,
                        bottom=0.10, wspace=0.3)
    plt.show()
