# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This file contains examples for order parameters.

This file is distributed as part of the documentation of PyRETIS.
"""
import numpy as np
from pyretis.orderparameter import OrderParameter


class PlaneDistanceX(OrderParameter):
    """A positional order parameter.

    This class defines a very simple order parameter which is
    the distance from a plane for a given particle.
    """

    def __init__(self, index, plane_pos):
        """Initialise the order parameter.

        Parameters
        ----------
        index : integer
            Selects the particle to use.

        plane_pos : float
            The location of the plane, along the x-axis.

        """
        txt = f'Distance from particle {index} to the plane at {plane_pos}'
        super().__init__(description=txt)
        self.index = index
        self.plane_position = plane_pos

    def calculate(self, system):
        """Calculate the order parameter."""
        pos = system.particles.pos[self.index]
        dist = -np.abs(pos[0] - self.plane_position)
        return [dist]
