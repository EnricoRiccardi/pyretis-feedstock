# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This package defines different simulations for use with PyRETIS.

The different simulations are defined as objects which inherit from
the base Simulation object defined in `simulation.py`. The simulation
object defines as simulation as a series of tasks to be executed,
typically at each step of the simulation. These tasks may produce
results which can be outputted to the user in some way.

Package structure
-----------------

Modules
~~~~~~~

md_simulation.py (:py:mod:`pyretis.simulation.md_simulation`)
    Defines simulation classes for molecular dynamics simulations.

mc_simulation.py (:py:mod:`pyretis.simulation.mc_simulation`)
    Define simulation classes for Monte Carlo simulations.

path_simulation.py (:py:mod:`pyretis.simulation.path_simulation`)
    Defines simulation classes for path simulations.

simulation.py (:py:mod:`pyretis.simulation.simulation`)
    Defines the Simulation class which is the base class for
    simulations.

simulation_task.py (:py:mod:`pyretis.simulation.simulation_task`)
    Defines classes for the handling of simulation tasks.

Important classes defined in this package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Simulation (:py:class:`.Simulation`)
    The base class for simulations.

SimulationTask (:py:class:`.SimulationTask`)
    A class for creating tasks for simulations.

SimulationSingleTIS (:py:class:`.SimulationSingleTIS`)
    A class for running a TIS simulation for a single ensemble.

SimulationRETIS (:py:class:`.SimulationRETIS`)
    A class for running a RETIS simulation for a set of ensembles.
"""
from .simulation import Simulation
from .simulation_task import SimulationTask
from .mc_simulation import UmbrellaWindowSimulation
from .md_simulation import SimulationMD, SimulationNVE, SimulationMDFlux
from .path_simulation import SimulationTIS, SimulationRETIS
