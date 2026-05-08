# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Definition of numerical integrators.

These integrators are typically used to integrate and propagate
Newtons equations of motion in time, the dynamics in molecular dynamics.
"""
import logging
from pyretis.engines import MDEngine
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


class VVIntegrator(MDEngine):
    """VVIntegrator(MDEngine).

    This class defines the Velocity Verlet integrator.

    Attributes
    ----------
    timestep : float
        The time step.
    half_timestep : float
        Half of the time step.
    description : string
        Description of the integrator.

    """

    def __init__(self, timestep,
                 description='The velocity verlet integrator'):
        """Set up the Velocity Verlet integrator.

        Parameters
        ----------
        timestep : float
            The time step in internal units.
        description : string
            Description of the integrator.

        """
        super().__init__(timestep, description, dynamics='NVE')
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
        imass = particles.imass
        particles.vel += self.half_timestep * particles.force * imass
        particles.pos += self.timestep * particles.vel
        system.potential_and_force()
        particles.vel += self.half_timestep * particles.force * imass


class Euler(MDEngine):
    """Euler(MDEngine).

    This class defines the Euler integrator.

    Attributes
    ----------
    timestep : float
        The time step.
    half_timestepsq : float
        Half of the squared time step.
    description : string
        Description of the integrator.

    """

    def __init__(self, timestep, description='The Euler integrator'):
        """Set up the Euler integrator.

        Parameters
        ----------
        timestep : float
            The time step in internal units.
        description : string
            Description of the integrator.

        """
        super().__init__(timestep, description, dynamics='NVE?')
        self.half_timestepsq = 0.5 * self.timestep**2

    def integration_step(self, system):
        """Euler integration, one time step.

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
        imass = particles.imass
        # update positions and velocities
        particles.pos += (self.timestep * particles.vel +
                          self.half_timestepsq * particles.force * imass)
        particles.vel += self.timestep * particles.force * imass
        # update forces
        system.potential_and_force()
