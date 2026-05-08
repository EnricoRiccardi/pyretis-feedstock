# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This package handles set-up of simulations from settings.

Package structure
-----------------

Modules
~~~~~~~

common.py (:py:mod:`pyretis.setup.common`)
    Common methods for handling settings. Defines a method to dynamically
    import methods and classes from user-specified modules.

createforcefield.py (:py:mod:`pyretis.setup.createforcefield`)
    Handle creation of force fields from input simulation settings.

createsimulation.py (:py:mod:`pyretis.setup.createsimulation`)
    Handle creation of simulations from input simulation settings.

createsystem.py (:py:mod:`pyretis.setup.createsystem`)
    Handle creation of systems from input simulation settings.

__init__.py
    This file. Handles imports for PyRETIS.

Important methods defined in this package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

create_engine (:py:func:`.create_engine`)
    Create an engine from input settings.

create_force_field (:py:func:`.create_force_field`)
    Create a force field from input settings.

create_orderparameter (:py:func:`.create_orderparameter`)
    Create an order parameter from input settings.

create_simulation (:py:func:`.create_simulation`)
    Create a simulation from input settings.

create_system (:py:func:`.create_system`)
    Create a system from input settings.
"""
from .common import create_orderparameter, create_engine
from .createsystem import create_system
from .createsimulation import create_simulation
from .createforcefield import create_force_field
