# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Definitions of simulation objects for molecular dynamics simulations.

This module contains definitions of classes for performing molecular
dynamics simulations.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SimulationNVE (:py:class:`.SimulationNVE`)
    Definition of a simple NVE simulation. The engine
    used for this simulation must have dynamics equal to NVE.

SimulationMD (:py:class:`.SimulationMD`)
    Definition of a simulation for running somply MD.

SimulationMDFlux (:py:class:`.SimulationMDFlux`)
    Definition of a simulation for determining the initial flux.
    This is used for calculating rates in TIS simulations.
"""
import logging
from pyretis.simulation.simulation import Simulation
from pyretis.core.particlefunctions import calculate_thermo
from pyretis.core.path import check_crossing
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = [
    'SimulationMD',
    'SimulationNVE',
    'SimulationMDFlux'
]


class SimulationMD(Simulation):
    """A generic MD simulation.

    This class is used to define a MD simulation without whistles and bells.

    """

    simulation_type = 'md'
    simulation_output = [
        {'type': 'energy', 'name': 'md-energy-file'},
        {'type': 'thermo-file', 'name': 'md-thermo-file'},
        {'type': 'traj-xyz', 'name': 'md-traj-file'},
        {'type': 'thermo-screen', 'name': 'md-thermo-screen'},
        {'type': 'order', 'name': 'md-order-file'},
    ]

    def __init__(self, ensemble, settings=None, controls=None):
        """Only add variable.

        Parameters
        ----------
        ensemble : dict
            It contains the simulations info

            * `system` : object like :py:class:`.System`
              This is the system we are investigating.
            * `engine` : object like :py:class:`.EngineBase`
              This is the integrator that is used to propagate
              the system in time.
            * `order_function` : object like :py:class:`.OrderParameter`,
              optional. A class that can be used to calculate an order
              parameter, if needed.

        settings : dict, optional
            This dictionary contains the settings for the simulation.
        controls: dict of parameters, optional
            It contains:

            * `steps` : int, optional
              The number of simulation steps to perform.
            * `startcycle` : int, optional
              The cycle we start the simulation on, can be useful if
              restarting.

        """
        super().__init__(settings, controls)
        self.engine = ensemble['engine']
        self.system = ensemble['system']
        self.order_function = ensemble.get('order_function')
        self.system.potential_and_force()  # make sure forces are defined.
        self.ensemble = ensemble

    def run(self):
        """Run the MD simulation.

        Yields
        ------
        results : dict
            The results from a single step in the simulation.

        """
        nsteps = 1 + self.cycle['endcycle'] - self.cycle['step']
        integ = self.engine.integrate({'system': self.system,
                                       'order_function': self.order_function},
                                      nsteps,
                                      thermo='full')
        for step in integ:
            if not self.first_step:
                self.cycle['step'] += 1
                self.cycle['stepno'] += 1
            results = {'cycle': self.cycle.copy()}
            results.update(self.execute_tasks())
            if self.first_step:
                self.first_step = False
            results.update(step)
            for task in self.output_tasks:
                task.output(results)
            self.write_restart()
            if self.soft_exit():
                yield results
                break
            yield results

    def __str__(self):
        """Return a string with info about the simulation."""
        msg = ['MD simulation']
        nstep = self.cycle['endcycle'] - self.cycle['startcycle']
        msg += [f'Number of steps to do: {nstep}']
        msg += [f'Molecular Dynamics engine: {self.engine}']
        msg += [f'Time step: {self.engine.timestep}']
        return '\n'.join(msg)


class SimulationNVE(SimulationMD):
    """A MD NVE simulation class.

    This class is used to define a NVE simulation. Compared with
    the :py:class:`.SimulationMD` we here require that the engine
    supports NVE dynamics.
    """

    simulation_type = 'md-nve'
    simulation_output = [
        {'type': 'energy', 'name': 'nve-energy-file'},
        {'type': 'thermo-file', 'name': 'nve-thermo-file'},
        {'type': 'traj-xyz', 'name': 'nve-traj-file'},
        {'type': 'thermo-screen', 'name': 'nve-thermo-screen'},
        {'type': 'order', 'name': 'nve-order-file'},
    ]

    def __init__(self, ensemble, settings=None, controls=None):
        """Initialise the NVE simulation object.

        Here we will set up the tasks that are to be performed in the
        simulation, such as the integration and thermodynamics
        calculation(s).

        Parameters
        ----------
        ensemble : dict
            It contains the simulations info

            * `system` : object like :py:class:`.System`
              This is the system we are investigating.
            * `engine` : object like :py:class:`.EngineBase`
              This is the integrator that is used to propagate the
              system in time.
            * `order_function` : object like :py:class:`.OrderParameter`, opt
              A class that can be used to calculate an order parameter,
              if needed.

        settings : dict, optional
            This dictionary contains the settings for the simulation.
        controls: dict of parameters, optional
            It contains:

            * `steps` : int, optional
              The number of simulation steps to perform.
            * `startcycle` : int, optional
              The cycle we start the simulation on.

        """
        super().__init__(ensemble, settings, controls)
        dyn = getattr(self.engine, "dynamics", "")
        if dyn.lower() != 'nve':
            logger.warning(
                'MD integrator %s (%s) might not be suitable'
                ' for NVE dynamics! Please make sure this is intentional.',
                self.engine.__class__,
                self.engine.description
            )

    def step(self):
        """Run a single simulation step."""
        system = self.ensemble['system']
        if self.first_step:
            system.potential_and_force()
            self.first_step = False
        else:
            self.cycle['step'] += 1
            self.cycle['stepno'] += 1
            self.engine.integration_step(system)
        results = {'cycle': self.cycle.copy(),
                   'thermo': calculate_thermo(system),
                   'system': system}
        if self.order_function:
            results['order'] = self.engine.calculate_order(self.ensemble)
        return results

    def __str__(self):
        """Return a string with info about the simulation."""
        msg = ['NVE simulation']
        nstep = self.cycle['endcycle'] - self.cycle['startcycle']
        msg += [f'Number of steps to do: {nstep}']
        msg += [f'MD engine: {self.engine}']
        msg += [f'Time step: {self.engine.timestep}']
        return '\n'.join(msg)

    def restart_info(self):
        """Return information which can be used to restart the simulation.

        Returns
        -------
        info : dict,
            Contains all the updated simulation settings and counters.

        """
        info = super().restart_info()
        info['engine'] = self.engine.restart_info()

        return info

    def load_restart_info(self, info):
        """Load the restart information."""
        super().load_restart_info(info)
        self.engine.load_restart_info(info['engine'])


class SimulationMDFlux(SimulationMD):
    """A simulation for obtaining the initial flux for TIS.

    This class is used to define a MD simulation where the goal is
    to calculate crossings in order to obtain the initial flux for a TIS
    calculation.

    """

    simulation_type = 'md-flux'
    simulation_output = [
        {'type': 'energy', 'name': 'flux-energy-file'},
        {'type': 'traj-xyz', 'name': 'flux-traj-file'},
        {'type': 'thermo-screen', 'name': 'flux-thermo-screen'},
        {'type': 'order', 'name': 'flux-order-file'},
        {'type': 'cross', 'name': 'flux-cross-file'},
    ]

    def __init__(self, ensemble, settings=None, controls=None):
        """Initialise the MD-Flux simulation object.

        Parameters
        ----------
        ensemble : dict
            It contains the simulations info

            * `system` : object like :py:class:`.System`
              The system to act on.
            * `engine` : object like :py:class:`.EngineBase`
              This is the integrator that is used to propagate the system
              in time.
            * `order_function` : object like :py:class:`.OrderParameter`
              The class used for calculating the order parameters.

        settings : dict, optional
            This dictionary contains the settings for the simulation.
        controls: dict of parameters, optional
            It contains:

            * `steps` : int, optional
              The number of simulation steps to perform.
            * `startcycle` : int, optional
              The cycle we start the simulation on, can be useful if
              restarting.

        """
        super().__init__(ensemble, settings, controls)

        # set up for initial crossing
        self.leftside_prev = None
        self.interfaces = settings['simulation'].get('interfaces')

    def run(self):
        """Run the MD simulation.

        Yields
        ------
        results : dict
            The results from a single step in the simulation.

        """
        nsteps = 1 + self.cycle['endcycle'] - self.cycle['step']
        leftside = None
        integ = self.engine.integrate(
            {'system': self.system,
             'order_function': self.order_function},
            nsteps,
            thermo='full',
        )
        for step in integ:
            if not self.first_step:
                self.cycle['step'] += 1
                self.cycle['stepno'] += 1
            else:
                self.first_step = False
            results = {'cycle': self.cycle.copy()}
            if leftside:
                self.leftside_prev = leftside
            leftside, cross = check_crossing(self.cycle['step'],
                                             step['order'][0],
                                             self.interfaces,
                                             self.leftside_prev)
            results['cross'] = cross
            results.update(step)

            for task in self.output_tasks:
                task.output(results)
            self.write_restart()
            if self.soft_exit():
                yield results
                break
            yield results

    def __str__(self):
        """Return a string with info about the simulation."""
        msg = ['MD-flux simulation']
        nstep = self.cycle['endcycle'] - self.cycle['startcycle']
        msg += [f'Number of steps to do: {nstep}']
        msg += [f'Molecular Dynamics engine: {self.engine}']
        msg += [f'Time step: {self.engine.timestep}']
        return '\n'.join(msg)

    def restart_info(self):
        """Return information which can be used to restart the simulation.

        Returns
        -------
        info : dict
            Contains all the updated simulation settings and counters.

        """
        info = super().restart_info()
        info['engine'] = self.engine.restart_info()
        info['leftside_prev'] = self.leftside_prev
        return info

    def load_restart_info(self, info):
        """Load the restart information."""
        super().load_restart_info(info)
        self.leftside_prev = info['leftside_prev']
        self.engine.load_restart_info(info['engine'])
