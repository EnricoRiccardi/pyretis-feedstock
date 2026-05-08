# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
# adapted by An Ghysels, April 10, 2019
"""This is a 1D flat potential"""
import logging
import numpy as np
from pyretis.forcefield.potential import PotentialFunction
logger = logging.getLogger(__name__)  # pylint: disable=C0103
logger.addHandler(logging.NullHandler())


class FlatWall1D(PotentialFunction):
    """Flat, 1D potential with 2 walls."""
    def __init__(self, desc='1D flat plus 2 walls'):
        """Initiate the potential."""

        super().__init__(dim=1, desc=desc)
        self.params = {}

    def potential(self, _):
        """Evaluate the potential.

        Returns
        -------
        out : float
            The potential energy.
        """
        return 0.

    def force(self, _):
        """Evaluate forces.

        Returns
        -------
        out[0] : numpy.array
            The calculated force.
        out[1] : numpy.array
            The virial, currently not implemented for this potential
        """
        # x = system.particles.pos[:, 0]

        forces = np.zeros(1)
        virial = np.zeros((1, 1))

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
        v_pot = self.potential(system)
        forces, virial = self.force(system)

        return v_pot, forces, virial
