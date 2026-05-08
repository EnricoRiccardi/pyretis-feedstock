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
"""
from pyretis.core.common import generic_factory
from .internal import MDEngine, Verlet, VelocityVerlet, Langevin
from .external import ExternalMDEngine
from .gromacs import GromacsEngine
from .gromacs2 import GromacsEngine2
from .cp2k import CP2KEngine
from .openmm import OpenMMEngine
from .lammps import LAMMPSEngine


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
    engine_map = {
        'velocityverlet': {'cls': VelocityVerlet},
        'verlet': {'cls': Verlet},
        'langevin': {'cls': Langevin},
        'gromacs': {'cls': GromacsEngine},
        'gromacs2': {'cls': GromacsEngine2},
        'cp2k': {'cls': CP2KEngine},
        'openmm': {'cls': OpenMMEngine},
        'lammps': {'cls': LAMMPSEngine},
    }
    return generic_factory(settings, engine_map, name='engine')
