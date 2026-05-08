# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This package defines the core PyRETIS tools.

The core tools are intended to define classes which are used in simulations.

Package structure
-----------------

Modules
~~~~~~~

__init__.py
    Import core functions from the modules.

box.py (:py:mod:`pyretis.core.box`)
    Definition of the simulation box class.

common.py (:py:mod:`pyretis.core.common`)
    Some common core methods, for instance for initiating classes.

montecarlo.py (:py:mod:`pyretis.core.montecarlo`)
    This module defines methods for performing Monte Carlo moves.

particlefunctions.py (:py:mod:`pyretis.core.particlefunctions`)
    Functions that operate on (a selection of) particles, for instance
    calculation of the kinetic temperature, pressure, momentum etc.

particles.py (:py:mod:`pyretis.core.particles`)
    Definition of the particle class which is used to represent
    a collection of particles.

pathensemble.py (:py:mod:`pyretis.core.pathensemble`)
    Definition of a class for a collection of paths (i.e. a path
    ensemble).

path.py (:py:mod:`pyretis.core.path`)
    This module defines functions and classes for paths.

properties.py (:py:mod:`pyretis.core.properties`)
    This module defines a class for a generic property.

random_gen.py (:py:mod:`pyretis.core.random_gen`)
    This module defines a class for generating random numbers.

retis.py (:py:mod:`pyretis.core.retis`)
    Module defining methods for performing replica exchange transition
    interface sampling.

system.py (:py:mod:`pyretis.core.system`)
    This module defines the system class which connects different
    parts (for instance box, forcefield and particles) into a single
    structure.

tis.py (:py:mod:`pyretis.core.tis`)
    This module contains methods used in the transition
    interface sampling algorithm.

units.py (:py:mod:`pyretis.core.units`)
    This module defines conversion between units.

Important classes defined in this package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BoxBase (:py:class:`.BoxBase`)
    A base class for a simulation box. The box will also
    handle the periodic boundaries.

System (:py:class:`.System`)
    A class which defines the system we are working with. This
    class contain a lot of information and is used to group the
    information into a structure which the simulations will make use
    of. Typically the system will contain a reference to a box,
    a list of particles and also a force field.

Particles (:py:class:`.Particles`)
    A class defining a list of particles. This will contain the
    positions, velocities and forces for the particles.

Path (:py:class:`.Path`)
    A class representing a path. The path contains snapshots with
    some additional information (energies and order parameters).

PathEnsemble (:py:class:`.PathEnsemble`)
    A class representing a collection of paths. The path ensemble
    will not store the full trajectories, but only a simplified
    representation.

PathEnsembleExt (:py:class:`.PathEnsembleExt`)
    A class representing external path ensembles. This handles
    additional bookkeeping related to external paths.

RandomGenerator (:py:class:`.RandomGenerator`)
    A class for generating random numbers.
"""
from .system import System
from .box import create_box
from .particles import Particles, ParticlesExt
from .path import Path
from .pathensemble import PathEnsemble, PathEnsembleExt
from .random_gen import RandomGenerator
