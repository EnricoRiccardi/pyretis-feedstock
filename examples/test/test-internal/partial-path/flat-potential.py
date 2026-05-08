# Copyright (c) 2015, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
# adapted by An Ghysels, March 12, 2019
# adapted by An Ghysels, April 10, 2019
"""This is a 1D example potential"""
import logging
import numpy as np
from pyretis.forcefield.potential import PotentialFunction
logger = logging.getLogger(__name__)  # pylint: disable=C0103
logger.addHandler(logging.NullHandler())


class FlatPotential1D(PotentialFunction):
    r"""Flat potential (PotentialFunction).

    This class defines a 1D potential which is simply flat.
    The potential energy (:math:`V_\text{pot}`) is given by

    .. math::

       V_\text{pot}(x) = 0.

    where :math:`x` gives the position.
    """

    def __init__(self, desc='1D flat potential no walls'):
        """Initiate the potential.

        Parameters
        ----------
        desc : string, optional
            Description of the force field.
        """
        super().__init__(dim=1, desc=desc)
        # these are the default parameters
        self.params = {}

    def potential(self, system):
        """Evaluate the potential.

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
        return 0.
        # x = system.particles.pos  # this is array, assume 1D (x)
        # v_pot = np.zeros_like(x)
        # return v_pot.sum()

    def force(self, system):
        """Evaluate forces.

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
            The virial, currently not implemented for this potential
        """
        x = system.particles.pos[:, 0]

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
