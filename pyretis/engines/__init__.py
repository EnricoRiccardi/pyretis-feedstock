# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Definition of engines.

This package defines engines for PyRETIS. The engines are responsible
for carrying out dynamics for a system. This can in principle both
be molecular dynamics or Monte Carlo dynamics. Typically, with RETIS,
this will be molecular dynamics in some form in order to propagate
the equations of motion and obtain new trajectories.


Package structure
-----------------

Modules
~~~~~~~

cp2k.py (:py:mod:`pyretis.engines.cp2k`)
    Defines an engine for use with CP2K.

engine.py (:py:mod:`pyretis.engines.engine`)
    Defines the base engine class.

external.py (:py:mod:`pyretis.engines.external`)
    Defines the interface for external engines.

gromacs.py (:py:mod:`pyretis.engines.gromacs`)
    Defines an engine for use with GROMACS.

gromacs2.py (:py:mod:`pyretis.engines.gromacs2`)
    Defines an engine for use with GROMACS. This is
    an alternative implementation which does not rely on
    continuously starting and stopping the GROMACS
    executable.

internal.py (:py:mod:`pyretis.engines.internal`)
    Defines internal PyRETIS engines.

lammps.py (:py:mod:`pyretis.engines.lammps`)
    Defines and engine for use with LAMMPS.

openmm.py (:py:mod:`pyretis.engines.openmm`)
    Defines an engine for use with OpenMM.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

engine_factory (:py:func:`.engine_factory`)
    A method to create engines from settings.

get_engine_class (:py:func:`.get_engine_class`)
    Return the engine class matching a name.

resolve_engine_class (:py:func:`.resolve_engine_class`)
    Resolve the engine class from full settings.

get_default_units (:py:func:`.get_default_units`)
    Return the default unit system for an engine.
"""
import logging
import os
from pyretis.core.common import generic_factory, import_from
from .internal import MDEngine, Verlet, VelocityVerlet, Langevin
from .external import ExternalMDEngine
from .gromacs import GromacsEngine
from .gromacs2 import GromacsEngine2
from .cp2k import CP2KEngine
from .openmm import OpenMMEngine
from .lammps import LAMMPSEngine

logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())

__all__ = [
    'MDEngine',
    'Verlet',
    'VelocityVerlet',
    'Langevin',
    'ExternalMDEngine',
    'GromacsEngine',
    'GromacsEngine2',
    'CP2KEngine',
    'OpenMMEngine',
    'LAMMPSEngine',
    'engine_factory',
    'get_engine_class',
    'resolve_engine_class',
    'get_default_units',
]

ENGINE_MAP = {
    'velocityverlet': {'cls': VelocityVerlet},
    'verlet': {'cls': Verlet},
    'langevin': {'cls': Langevin},
    'gromacs': {'cls': GromacsEngine},
    'gromacs2': {'cls': GromacsEngine2},
    'cp2k': {'cls': CP2KEngine},
    'openmm': {'cls': OpenMMEngine},
    'lammps': {'cls': LAMMPSEngine},
}


def engine_factory(settings):
    """Create an engine according to the given settings.

    This function is included as a convenient way of setting up and
    selecting an engine. It will return the created engine.

    Parameters
    ----------
    settings : dict
        This defines how we set up and select the engine.

    Returns
    -------
    out : object like :py:class:`.EngineBase`
        The object representing the engine to use in a simulation.

    """
    return generic_factory(settings, ENGINE_MAP, name='engine')


def get_engine_class(name):
    """Return the class matching a given engine name."""
    try:
        return ENGINE_MAP[name.lower()]['cls']
    except (AttributeError, KeyError):
        return None


def resolve_engine_class(settings):
    """Resolve an engine class from the given settings."""
    if 'engine' in settings:
        engine_settings = settings['engine']
        exe_path = settings.get('simulation', {}).get('exe_path')
    else:
        engine_settings = settings
        exe_path = engine_settings.get('exe_path')

    engine_object = engine_settings.get('obj')
    if engine_object is not None:
        return engine_object.__class__

    klass = engine_settings.get('class')
    if klass is None:
        return None

    module = engine_settings.get('module')
    if module is None:
        cls = get_engine_class(klass)
        if cls is None:
            msg = f'Could not resolve engine class "{klass}"'
            raise ValueError(msg)
        return cls

    module_path = module
    if not os.path.isfile(module_path) and exe_path is not None:
        module_path = os.path.join(exe_path, module)
    return import_from(module_path, klass)


def get_default_units(settings):
    """Return the default unit system for the selected engine."""
    klass = resolve_engine_class(settings)
    if klass is None:
        return None
    return klass.get_default_units(settings)
