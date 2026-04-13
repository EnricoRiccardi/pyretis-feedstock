# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Plot raw data from a simulation."""
# pylint: disable=invalid-name
import numpy as np
import colorama
from matplotlib import pyplot as plt
from matplotlib.widgets import Button, RadioButtons
from matplotlib.cm import get_cmap
from pyretis.core import create_box, System, Particles
from pyretis.core.tis import shoot, time_reversal
from pyretis.core.retis import retis_swap
from pyretis.initiation import initiate_path_simulation
from pyretis.inout.settings import parse_settings_file
from pyretis.setup import create_simulation, create_force_field
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def plot_path(path, axi, axj, color, alphai=0.9, alphaj=0.9, lsj='-'):
    """Plot data for paths in ensembles."""
    pos = []
    order = []
    time = []
    for i, phasepoint in enumerate(path.phasepoints):
        pos.append(phasepoint.particles.get_pos()[0])
        order.append(phasepoint.order[0])
        time.append(i)
    pos = np.array(pos)
    order = np.array(order)
    axi.plot(pos[:, 0], pos[:, 1], lw=3, alpha=alphai,
             color=color)
    axj.plot(time, order, lw=3, alpha=alphaj, color=color, ls=lsj)
    axi.scatter(pos[0, 0], pos[0, 1], s=50, marker='x', color=color,
                alpha=alphai)
    axj.scatter(time[0], order[0], s=50, marker='x', color=color,
                alpha=alphaj)
    if order[-1] < order[-2]:
        end = '^'
    else:
        end = 'v'
    axi.scatter(pos[-1, 0], pos[-1, 1], s=50, marker=end, color=color,
                alpha=alphai)
    axj.scatter(time[-1], order[-1], s=50, marker=end, color=color,
                alpha=alphaj)


class PlotHelper:
    """A class for handling the plotting."""

    def __init__(self, simulation, settings, axi, axj, colors):
        """Set up the plotter."""
        self.simulation = simulation
        self.axi = axi
        self.axj = axj
        self.colors = colors
        self.potential_setup(settings)
        self.update = None
        self.last_result = None
        self.idx = 0

    def potential_setup(self, settings):
        """Set up the potential function."""
        # Set up raw data:
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
        self.pot = pot
        self.xpos = xpos
        self.ypos = ypos
        self.interfaces = settings['simulation']['interfaces']
        self.inter_a = settings['orderparameter']['inter_a']
        self.inter_b = settings['orderparameter']['inter_b']
        self.minx = minx
        self.maxx = maxx
        self.miny = miny
        self.maxy = maxy

    def plot_potential(self, axi, axj):
        """Plot the potential in the given axis."""
        axi.contourf(self.xpos, self.ypos, self.pot, 10,
                     cmap=get_cmap('viridis'), alpha=0.8)
        # Add interfaces:
        for inter in self.interfaces:
            axi.axvline(x=inter, lw=2, ls=':', color='#262626', alpha=0.8)
            axj.axhline(y=inter, lw=2, ls=':', color='#262626', alpha=0.8)
        extra_int = [self.inter_a, self.inter_b]
        for inter in extra_int:
            axi.axvline(x=inter, lw=2, ls=':', alpha=0.8, color='#BAACE6')
            axj.axhline(y=inter, lw=2, ls=':', alpha=0.8, color='#BAACE6')
        axi.set_xlim((self.minx, self.maxx))
        axi.set_ylim((self.miny, self.maxy))
        axi.set_xlabel((r'Position ($x$)'), fontsize='large')
        axi.set_ylabel((r'Position ($y$)'), fontsize='large')
        axj.set_xlabel((r'Step number'), fontsize='large')
        axj.set_ylabel((r'Order parameter ($\lambda$)'), fontsize='large')

    def plot_accepted(self, idx=None):
        """Plot the accepted paths."""
        ensembles = self.simulation.ensembles
        for i, ensemble in enumerate(ensembles):
            if idx is None:
                alpha = 0.9
            else:
                if i in idx:
                    alpha = 0.9
                else:
                    alpha = 0.25
            plot_path(ensemble['path_ensemble'].last_path, self.axi, self.axj,
                      self.colors[i], alphai=0.9, alphaj=alpha)

    def do_shoot(self, event):
        """Perform shooting moves."""
        ensembles = self.simulation.ensembles
        ensemble = ensembles[self.idx]
        accept, trial, status = shoot(
            ensemble,
            self.simulation.settings['tis'],
            ensemble['path_ensemble'].start_condition
        )
        self.axi.clear()
        self.axj.clear()
        self.plot_potential(self.axi, self.axj)
        self.plot_accepted((self.idx,))
        plot_path(trial, self.axi, self.axj, self.colors[-1], lsj='--')
        sidx = trial.generated[2]
        spoint = ensemble['path_ensemble'].last_path.phasepoints[sidx]
        pos = spoint.particles.get_pos()[0]
        order = spoint.order[0]
        self.axi.scatter(pos[0], pos[1], marker='o', s=75, alpha=0.9,
                         color=self.colors[self.idx])
        self.axj.scatter(sidx, order, marker='o', s=75, alpha=0.9,
                         color=self.colors[self.idx])
        self.last_result = ('Shoot', ensemble, accept, trial, status, self.idx)
        self.axj.relim()
        self.axj.autoscale_view()
        plt.draw()

    def do_timereversal(self, event):
        """Perform time reversal."""
        ensembles = self.simulation.ensembles
        ensemble = ensembles[self.idx]
        accept, trial, status = time_reversal(
            ensemble['path_ensemble'].last_path,
            ensemble['order_function'],
            ensemble['path_ensemble'].interfaces,
            ensemble['path_ensemble'].start_condition
        )
        self.axi.clear()
        self.axj.clear()
        self.plot_potential(self.axi, self.axj)
        self.plot_accepted((self.idx,))
        plot_path(trial, self.axi, self.axj, self.colors[-1], lsj='--')
        self.last_result = ('Time-reversal', ensemble, accept, trial,
                            status, self.idx)
        self.axj.relim()
        self.axj.autoscale_view()
        plt.draw()

    def do_accept_last_move(self, event):
        """Try to accept the last move."""
        self.accept_last_move()

    def accept_last_move(self):
        """Check if we can accept the last move."""
        if self.last_result is not None:
            move, ensemble, accept, trial, status, _ = self.last_result
            if move == 'Swap':
                # idxs = (idx, idx + 1)
                logger.info('Swap: {} <-> {}'.format(
                        ensemble[0]['path_ensemble'].ensemble_name,
                        ensemble[1]['path_ensemble'].ensemble_name
                    ))
            else:
                # idxs = (idx,)
                ensemble['path_ensemble'].add_path_data(trial, status)
                logger.info(
                    f"In ensemble: {ensemble['path_ensemble'].ensemble_name}")
            if accept:
                logger.info(f'\t{move} was accepted!')
                logger.info(f'\t{status}')
            else:
                logger.error(f'\t{move} was rejected!')
                logger.error(f'\t{status}')
            self.axi.clear()
            self.axj.clear()
            self.plot_potential(self.axi, self.axj)
            self.plot_accepted()
            self.last_result = None
            self.axj.relim()
            self.axj.autoscale_view()
            plt.draw()

    def set_idx(self, label):
        """Just to define indices for the ensembles."""
        idx = {'$[0^-]$': 0, '$[0^+]$': 1, '$[1^+]$': 2,
               '$[2^+]$': 3, '$[3^+]$': 4, '$[4^+]$': 5,
               '$[5^+]$': 6}
        self.idx = idx.get(label, 0)

    def do_swap(self, event):
        """Perform swapping when button is pressed."""
        if self.idx == len(self.simulation.ensembles) - 1:
            return
        accept, trial, status = retis_swap(
            self.simulation.ensembles,
            self.idx,
            self.simulation.settings,
            0
        )
        ensemble1 = self.simulation.ensembles[self.idx]
        ensemble2 = self.simulation.ensembles[self.idx + 1]
        self.last_result = ('Swap', (ensemble1, ensemble2), accept, trial,
                            status, self.idx)
        self.axi.clear()
        self.axj.clear()
        self.plot_potential(self.axi, self.axj)
        self.plot_accepted((self.idx, self.idx+1))
        self.axj.relim()
        self.axj.autoscale_view()
        plt.draw()


def set_up_simulation(settings):
    """Do all the set-ups to create the simulation."""
    simulation = create_simulation(settings)
    # Also, do the initialisation here:
    for i, _ in enumerate(initiate_path_simulation(simulation, settings)):
        ensemble = simulation.ensembles[i]['path_ensemble']
        name = ensemble.ensemble_name
        logger.info(f'Info about ensemble {name}:')
        logger.info(ensemble)
        logger.info('Info about the initial path:')
        logger.info(ensemble.last_path)
    return simulation


def main(inputfile='retis.rst'):
    """Run the whole interactive plotting."""
    # Setup simulation:
    sim_settings = parse_settings_file(inputfile)
    simulation = set_up_simulation(sim_settings)
    # Set up for plotting:
    fig = plt.figure(figsize=(12, 6))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    axbt = fig.add_axes([0.20, 0.12, 0.1, 0.075])
    axbt2 = fig.add_axes([0.30, 0.12, 0.1, 0.075])
    axbt3 = fig.add_axes([0.40, 0.12, 0.1, 0.075])
    axbt4 = fig.add_axes([0.30, 0.04, 0.1, 0.075])

    cmap = get_cmap(name='tab10')
    colors = cmap(np.linspace(0, 1, len(simulation.ensembles)+2))

    myplot = PlotHelper(simulation, sim_settings, ax1, ax2, colors)

    # Add buttons:
    button1 = Button(axbt, 'Shoot')
    button1.on_clicked(myplot.do_shoot)
    button2 = Button(axbt2, 'Time reversal')
    button2.on_clicked(myplot.do_timereversal)
    button3 = Button(axbt3, 'Swap')
    button3.on_clicked(myplot.do_swap)
    button4 = Button(axbt4, 'Accept?')
    button4.on_clicked(myplot.do_accept_last_move)

    # Add radio button:
    rax = fig.add_axes([0.02, 0.2, 0.1, 0.6],
                       frameon=True)
    names = []
    for ensemble in simulation.ensembles:
        names.append(f"${ensemble['path_ensemble'].ensemble_name}$")
    radio = RadioButtons(rax, names, activecolor='#262626')
    for circle in radio.circles:
        circle.set_radius(0.02)
    for i, lab in enumerate(radio.labels):
        lab.set_fontsize(12)
        lab.set_color(colors[i])
    radio.on_clicked(myplot.set_idx)

    myplot.plot_potential(ax1, ax2)
    plt.subplots_adjust(right=0.95, left=0.20, top=0.95,
                        bottom=0.3, wspace=0.3)
    for i, ensemble in enumerate(simulation.ensembles):
        plot_path(ensemble['path_ensemble'].last_path, ax1, ax2, colors[i])
    plt.show()


if __name__ == '__main__':
    colorama.init(autoreset=True)
    main(inputfile='retis-x.rst')
