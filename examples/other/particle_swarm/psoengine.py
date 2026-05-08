# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A custom engine for particle swarm optimization."""
import numpy as np
from pyretis.engines import MDEngine


class PSOEngine(MDEngine):
    """Perform particle swarm optimization."""

    def __init__(self, inertia, accp, accg):
        """Set up the engine.

        Parameters
        ----------
        intertia : float
            The intertia factor in the velocity equation
            of motion.
        accp : float
            The acceleration for the previous best term. "The congnitive term".
        accg : float
            The acceleration for the global best term. "The social term".

        """
        super().__init__(1, 'Particle Swarm Optimization')
        self.inertia = inertia
        self.accp = accp
        self.accg = accg
        self.pbest = None
        self.pbest_pot = None
        self.gbest = None
        self.gbest_pot = None

    def integration_step(self, system):
        """Perform one step for the PSO algorithm."""
        particles = system.particles
        if self.pbest is None:
            self.pbest = np.copy(particles.pos)
            self.pbest_pot = system.potential()
        if self.gbest is None:
            pot = system.potential()
            idx = np.argmin(pot)
            self.gbest = particles.pos[idx]
            self.gbest_pot = pot[idx]

        rnd1 = np.random.uniform()
        rnd2 = np.random.uniform()

        particles.vel = (self.inertia * particles.vel +
                         rnd1 * self.accp * (self.pbest - particles.pos) +
                         rnd2 * self.accg * (self.gbest - particles.pos))
        particles.pos += particles.vel
        particles.pos = system.box.pbc_wrap(particles.pos)

        pot = system.potential()

        # Update global?
        idx = np.argmin(pot)
        if pot[idx] < self.gbest_pot:
            self.gbest_pot = pot[idx]
            self.gbest = particles.pos[idx]
        # Update for individuals:
        idx = np.where(pot < self.pbest_pot)[0]
        self.pbest[idx] = np.copy(particles.pos[idx])
        self.pbest_pot[idx] = pot[idx]
        return self.gbest_pot, self.gbest
