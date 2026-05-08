# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Example of using a integration routine implemented in FORTRAN."""
import logging
import os
import sys
from pyretis.engines import MDEngine
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
try:
    from vvintegrator import vvintegrator
except ImportError:
    MSG = ('Could not import external FORTRAN library.'
           '\nPlease compile with "make"!')
    logger.critical(MSG)
    raise ImportError(MSG)


__all__ = ['VelocityVerletF']


class VelocityVerletF(MDEngine):
    """VelocityVerletF(MDEngine).

    This class defines the Velocity Verlet integrator.

    Attributes
    ----------
    timestep : float
        The time step.
    half_timestep : float
        The half of the timestep.
    desc : string
        Description of the integrator.

    """

    def __init__(self, timestep,
                 desc='The velocity verlet integrator (FORTRAN)'):
        """Set up the Velocity Verlet integrator.

        Parameters
        ----------
        timestep : float
            The time step.
        desc : string
            Description of the integrator.

        """
        super().__init__(timestep, desc, dynamics='NVE')
        self.half_timestep = self.timestep * 0.5

    def integration_step(self, system):
        """Velocity Verlet integration, one time step.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system to integrate/act on. Assumed to have a particle
            list in `system.particles`.

        Returns
        -------
        out : None
            Does not return anything, but alters the state of the given
            `system`.

        """
        particles = system.particles
        particles.pos, particles.vel = vvintegrator.step1(
            particles.pos,
            particles.vel,
            particles.force,
            particles.imass,
            self.timestep,
            self.half_timestep
        )
        system.potential_and_force()
        particles.vel = vvintegrator.step2(
            particles.vel,
            particles.force,
            particles.imass,
            self.half_timestep
        )
