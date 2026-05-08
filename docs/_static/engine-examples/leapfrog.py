# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A custom engine for Leap Frog integration."""
from pyretis.engines import MDEngine


class LeapFrog(MDEngine):
    """Perform Leap Frog integration."""

    def __init__(self, timestep):
        """Initiate the engine.

        Parameters
        ----------
        timestep : float
            Set the time step for the integrator.
        """
        super().__init__(timestep, 'Leap Frog integrator', dynamics='NVE')
        self.half_timestep = self.timestep * 0.5
        self.half_timestep_sq = 0.5 * self.timestep**2

    def integration_step(self, system):
        """Perform one step for the Leap Frog integrator."""
        particles = system.particles
        imass = particles.imass

        # get current acceleration:
        acc_t = particles.force * imass
        # update positions:
        particles.pos += (self.timestep * particles.vel +
                          self.half_timestep_sq * acc_t)
        # update forces:
        system.potential_and_force()
        # update acceleration:
        acc_t2 = particles.force * imass
        # update velocities:
        particles.vel += self.half_timestep * (acc_t + acc_t2)
