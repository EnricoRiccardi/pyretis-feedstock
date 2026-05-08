# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Definition of simulation objects for Monte Carlo simulations.

This module defines some classes and functions for performing
Monte Carlo simulations.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

UmbrellaWindowSimulation (:py:class:`.UmbrellaWindowSimulation`)
    Defines a simulation for performing umbrella window simulations.
    Several umbrella window simulations can be joined to perform a
    umbrella simulation.
"""
import numpy as np
from pyretis.core.montecarlo import max_displace_step
from pyretis.simulation.simulation import Simulation
from pyretis.core.random_gen import create_random_generator

__all__ = ['UmbrellaWindowSimulation']


def mc_task(rgen, system, maxdx):
    """Perform a Monte Carlo displacement move.

    Here, a displacement step will be carried out and
    the trial move will be accepted/rejected and the
    positions and potential energy will be updated as needed.

    Parameters
    ----------
    rgen : object like :py:class`.RandomGenerator`
        This object is used for generating random numbers.
    system : object like :py:class:`.System`
        The system we act on.
    maxdx : float
        Maximum displacement step for the Monte Carlo move.

    """
    accepted_r, _, trial_r, v_trial, status = max_displace_step(rgen, system,
                                                                maxdx)
    if status:
        system.particles.pos = accepted_r
        system.particles.vpot = v_trial
    return accepted_r, v_trial, trial_r, status


class UmbrellaWindowSimulation(Simulation):
    """This class defines an Umbrella simulation.

    The Umbrella simulation is a special case of
    the simulation class with settings to simplify the
    execution of the umbrella simulation.

    Attributes
    ----------
    ensemble : dict
        It contains the simulations info
          *  system : object like :py:class:`.System`
                 The system to act on.
          *  rgen : object like :py:class:`.RandomGenerator`
                 Object to use for random number generation.
    umbrella : list, [float, float]
        The umbrella window to consider.
    overlap : float
        The position we have to cross before the simulation ends.
    maxdx : float
        Defines the maximum movement allowed in the Monte Carlo
    steps : int, optional
        The number of simulation steps to perform.
    startcycle : int, optional
        The cycle we start the simulation on, can be useful if
        restarting.

    """

    simulation_type = 'umbrella-window'

    def __init__(self, ensemble, settings, controls=None):
        """Initialise the umbrella simulation simulation.

        Parameters
        ----------
        ensemble : dict
            It contains the simulations info
              *  system : object like :py:class:`.System`
                     The system to act on.
              *  rgen : object like :py:class:`.RandomGenerator`
                     Object to use for random number generation.
        settings : dict
            Contains all the simulation settings.
        controls: dict of parameters to control the simulations. Optional
            It can contain:
              *  mincycle : int, optional
                     The *MINIMUM* number of cycles to perform. Note that in
                     base `Simulation` class this is the *MAXIMUM* number of
                     cycles to perform. The meaning is redefined in this class
                     by overriding `self.simulation_finished`.
              *  startcycle : int, optional
                     The current simulation cycle, i.e. where we start.

        """
        super().__init__(settings, controls)
        self.umbrella = settings['simulation']['umbrella']
        self.overlap = settings['simulation']['overlap']
        self.maxdx = settings['simulation']['maxdx']
        self.rgen = ensemble.get(
            'rgen',
            create_random_generator(settings['simulation']))

        self.system = ensemble['system']
        task_monte_carlo = {'func': mc_task,
                            'args': [self.rgen, self.system, self.maxdx],
                            'result': 'displace_step'}
        self.add_task(task_monte_carlo)
        self.first_step = False

        if 'restart' in settings:
            self.load_restart_info(settings['restart'])

    def is_finished(self):
        """Check if the simulation is done.

        In the umbrella simulation, the simulation is finished when we
        cycle is larger than maxcycle and all particles have
        crossed self.overlap.

        Returns
        -------
        out : boolean
            True if the simulation is finished, False otherwise.

        """
        return (self.cycle['step'] > self.cycle['endcycle'] and
                np.all(self.system.particles.pos > self.overlap))

    def __str__(self):
        """Return some info about the simulation as a string."""
        msg = ['Umbrella window simulation']
        msg += [f'Umbrella: {self.umbrella}, Overlap: {self.overlap}.']
        msg += [f"Minimum number of cycles: {self.cycle['endcycle']}"]
        return '\n'.join(msg)

    def restart_info(self):
        """Return information which can be used to restart the simulation.

        Returns
        -------
        info : dict,
            Contains all the updated simulation settings and counters.

        """
        info = super().restart_info()

        info['simulation'].update(self.rgen.get_state())
        info['umbrella'] = self.umbrella
        info['overlap'] = self.overlap
        info['type'] = self.simulation_type
        return info

    def load_restart_info(self, info):
        """Load the restart information."""
        super().load_restart_info(info)

        self.umbrella = info.get('umbrella', self.umbrella)
        self.overlap = info.get('overlap', self.overlap)
        self.rgen = create_random_generator(info['simulation'])
