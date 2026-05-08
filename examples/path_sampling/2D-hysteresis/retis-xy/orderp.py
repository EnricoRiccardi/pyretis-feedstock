# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This file defines the order parameter used for the WCA example."""
import logging
import numpy as np
from pyretis.orderparameter import OrderParameter
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


class OrderXY(OrderParameter):
    """OrderXY(OrderParameter).

    This class defines a 2D order parameter such that x + y = constant.

    Attributes
    ----------
    index : integer
        This selects the particle to use for the order parameter.
    inter_a : float
        An interface such that we are in state A for postions < inter_a.
    inter_b : float
        An interface such that we are in state B for postions > inter_b.
    energy_a : float
        An energy such that we are in state A for potential energy < energy_a.
    energy_b : float
        An energy such that we are in state A for potential energy < energy_b.

    """

    def __init__(self, index, inter_a, inter_b, energy_a, energy_b):
        """Set up the order parameter.

        Parameters
        ----------
        index : tuple of ints
            This is the indices of the atom we will use the position of.
        inter_a : float
            An interface such that we are in state A for postions < inter_a.
        inter_b : float
            An interface such that we are in state B for postions > inter_b.
        energy_a : float
            An energy such that we are in state A for
            potential energy < energy_a.
        energy_b : float
            An energy such that we are in state A for
            potential energy < energy_b.

        """
        super().__init__(description='2D->1D projection')
        self.index = index
        x_1 = 0.2
        y_1 = 0.4
        x_0 = -0.2
        y_0 = -0.4
        self.x_0 = x_0
        self.y_0 = y_0
        self.origin = np.array([x_0, y_0])
        self.vec = np.array([x_1 - x_0, y_1 - y_0])
        self.vec /= np.sqrt(np.dot(self.vec, self.vec))
        self.inter_a = inter_a
        self.inter_b = inter_b
        self.energy_a = energy_a
        self.energy_b = energy_b

    def calculate(self, system):
        """Calculate the order parameter.

        Here, the order parameter is just the distance between two
        particles.

        Parameters
        ----------
        system : object like :py:class:`.System`
            This object is used for the actual calculation, typically
            only `system.particles.pos` and/or `system.particles.vel`
            will be used. In some cases `system.forcefield` can also be
            used to include specific energies for the order parameter.

        Returns
        -------
        out : list of floats
            The order parameter.

        """
        pos = system.particles.pos[self.index]
        vec = np.array([pos[0] - self.x_0, pos[1] - self.y_0])
        proj = np.dot(vec, self.vec)
        if proj < self.inter_a:
            if system.particles.vpot > self.energy_a:
                proj = self.inter_a
        elif proj > self.inter_b:
            if system.particles.vpot > self.energy_b:
                proj = self.inter_b
        return [proj]
