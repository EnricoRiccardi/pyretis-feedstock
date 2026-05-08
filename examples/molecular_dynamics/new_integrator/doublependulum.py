# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A double pendulum potential."""
import logging
import numpy as np
from pyretis.forcefield.potential import PotentialFunction
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


class DoublePendulum(PotentialFunction):
    r"""DoublePendulum(PotentialFunction).

    Attributes
    ----------
    l1 : float
        Length of first wire.
    l2 : float
        Length of second wire.
    g : float
        The acceleration of gravity.

    """

    def __init__(self, l1=1.0, l2=1.0, g=1.0, desc='2D double pendulum'):
        """Set up the double pendulum potential.

        Parameters
        ----------
        l1 : float, optional
            Length of first wire.
        l1 : float, optional
            Length of second wire.
        g : float, optional
            The acceleration of gravity.
        desc : string, optional
            Description of the force field.

        """
        super().__init__(dim=2, desc=desc)
        self.params = {'g': g, 'l1': l1, 'l2': l2}
        self.l1 = l1
        self.l2 = l2
        self.g = g

    def potential(self, system):
        """Evaluate the potential energy.

        Parameters
        ----------
        system : object like `System`
            The system we evaluate the potential for. Here, we
            make use of the positions only.

        Returns
        -------
        out : float
            The potential energy.

        """
        costheta1 = np.cos(system.particles.pos[0][0])
        costheta2 = np.cos(system.particles.pos[1][0])
        m1 = system.particles.mass[0][0]
        m2 = system.particles.mass[1][0]
        return -self.g * ((m1 + m2) * self.l1 * costheta1 +
                          m2 * self.l2 * costheta2)

    def force(self, system):
        """Evaluate forces for the double pendulum.

        Parameters
        ----------
        system : object like `System`
            The system we evaluate the potential for. Here, we
            make use of the positions only.

        Returns
        -------
        out[0] : numpy.array
            The calculated force.
        out[1] : numpy.array
            The virial, currently not implemented for this potential.

        """
        g = self.g
        l1 = self.l1
        l2 = self.l2
        pos = system.particles.pos
        vel = system.particles.vel
        m1 = system.particles.mass[0][0]
        m2 = system.particles.mass[1][0]
        theta1 = pos[0][0]
        theta2 = pos[1][0]
        dtheta1 = vel[0][0]
        dtheta2 = vel[1][0]
        dth = theta1 - theta2
        M = m2 / (m1 + m2)
        ll = l2 / l1
        wsq = g / l1
        denom = ll*(1.0 - M*np.cos(dth)**2)
        atheta1 = wsq*ll*(-np.sin(theta1) + M*np.cos(dth)*np.sin(theta2))
        atheta1 -= M*ll*np.sin(dth)*(dtheta1**2*np.cos(dth) +
                                     ll*dtheta2**2)
        atheta1 /= denom

        atheta2 = wsq*np.cos(dth)*np.sin(theta1) - wsq*np.sin(theta2)
        atheta2 += (dtheta1**2 + M*ll*dtheta2**2*np.cos(dth))*np.sin(dth)
        atheta2 /= denom

        forces = np.zeros((2, self.dim))
        forces[0, 0] = atheta1
        forces[1, 0] = atheta2
        virial = np.zeros((self.dim, self.dim))  # just return zeros here
        return forces, virial

    def potential_and_force(self, system):
        """Evaluate the potential and the force.

        Parameters
        ----------
        system : object like `System`
            The system we evaluate the potential for. Here, we
            make use of the positions only.

        Returns
        -------
        out[0] : float
            The potential energy as a float.
        out[1] : numpy.array
            The force as a numpy.array of the same shape as the
            positions in `particles.pos`.
        out[2] : numpy.array
            The virial, currently not implemented for this potential.

        """
        forces, virial = self.force(system)
        pot = self.potential(system)
        return pot, forces, virial
