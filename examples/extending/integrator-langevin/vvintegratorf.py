# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""The Langevin integrator implemented in FORTRAN."""
import logging
import os
import sys
import numpy as np
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


__all__ = ['LangevinF']


class LangevinF(MDEngine):
    """LangevinF(MDEngine).

    This class defines the Langevin integrator.

    Attributes
    ----------
    timestep : float
        The time step.
    gamma : float
        The gamma parameter for the Langevin integrator.
    high_friction : boolean
        Determines if we are in the high_friction limit and should
        do the over-damped version.

    """

    def __init__(self, timestep, gamma, high_friction, seed):
        """Set up the Velocity Verlet integrator.

        Parameters
        ----------
        timestep : float
            The time step.
        gamma : float
            The gamma parameter for the Langevin thermostat
        high_friction : boolean
            Determines if we should do overdamped or not.
        seed : integer
            A seed for the random number generator.

        """
        super().__init__(timestep, 'The Langevin integrator (FORTRAN)',
                         dynamics='stochastic')
        self.gamma = gamma
        self.high_friction = high_friction
        self.param_high = {'sigma': None, 'bddt': None}
        self.param_iner = {'c0': None, 'a1': None, 'a2': None, 'b1': None,
                           'b2': None, 's12os11': None, 'sqrts11': None,
                           'sqrtsos11': None}
        self.init_params = True
        size = vvintegrator.get_seed_size()
        seeds = np.array([seed + i for i in range(size)], dtype='i4')
        vvintegrator.seed_random_generator(seeds)

    def _initiate_parameters(self, system):
        """Set up Langevin parameters.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system to integrate/act on. Assumed to have a particle
            list in `system.particles`.

        """
        beta = system.temperature['beta']
        imasses = system.particles.imass
        if self.high_friction:
            self.param_high['sigma'] = np.sqrt(2.0 * self.timestep *
                                               imasses/(beta * self.gamma))
            self.param_high['bddt'] = self.timestep * imasses / self.gamma
        else:
            gammadt = self.gamma * self.timestep
            exp_gdt = np.exp(-gammadt)
            s11 = ((self.timestep * imasses / (beta * self.gamma)) *
                   (2. - (3. - 4. * exp_gdt + exp_gdt**2) / gammadt))
            sig_r = np.sqrt(s11)
            s22 = (1.0 - exp_gdt**2) * imasses / beta
            s12 = (imasses / (beta * self.gamma)) * (1. - exp_gdt)**2
            self.param_iner['s12os11'] = s12 / s11
            self.param_iner['sqrts11'] = sig_r
            self.param_iner['sqrtsos11'] = np.sqrt((s11 * s22 - s12**2) / s11)
            if self.gamma > 0.0:
                c_0 = exp_gdt
                c_1 = (1.0 - c_0) / gammadt
                c_2 = (1.0 - c_1) / gammadt
            else:
                c_0, c_1, c_2 = 1.0, 1.0, 0.5
            self.param_iner['c0'] = c_0
            self.param_iner['a1'] = c_1 * self.timestep
            self.param_iner['a2'] = c_2 * self.timestep**2 * imasses
            self.param_iner['b1'] = (c_1 - c_2) * self.timestep * imasses
            self.param_iner['b2'] = c_2 * self.timestep * imasses

    def integration_step(self, system):
        """Langevin integration, one time step.

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
        if self.init_params:
            self._initiate_parameters(system)
            self.init_params = False
        if self.high_friction:
            return self.integration_step_overdamped(system)
        return self.integration_step_inertia(system)

    def integration_step_overdamped(self, system):
        """Over damped Langevin integration, one time step.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system to integrate/act on. Assumed to have a particle
            list in ``system.particles``.

        Returns
        -------
        out : None
            Does not return anything, but alters the state of the given
            `system`.

        """
        system.force()  # update forces
        particles = system.particles
        particles.pos, particles.vel = vvintegrator.overdamped(
            particles.pos,
            particles.vel,
            particles.force,
            self.param_high['bddt'],
            self.param_high['sigma']
        )
        system.potential()

    def integration_step_inertia(self, system):
        """Langevin integration, one time step.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system to integrate/act on. Assumed to have a particle
            list in ``system.particles``.

        Returns
        -------
        out : None
            Does not return anything, but alters the state of the given
            `system`.

        """
        particles = system.particles
        particles.pos, particles.vel = vvintegrator.inertia1(
            particles.pos,
            particles.vel,
            particles.force,
            self.gamma,
            self.param_iner['c0'],
            self.param_iner['a1'],
            self.param_iner['a2'],
            self.param_iner['b1'],
            self.param_iner['s12os11'],
            self.param_iner['sqrts11'],
            self.param_iner['sqrtsos11']
        )
        system.force()  # update forces
        particles.vel = vvintegrator.inertia2(
            particles.vel,
            particles.force,
            self.param_iner['b2']
        )
        system.potential()
