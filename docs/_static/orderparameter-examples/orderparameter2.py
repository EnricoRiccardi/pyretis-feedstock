# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This file contains examples for order parameters.

This file is distributed as part of the documentation of PyRETIS.

Note: mdtraj needs to be installed on top of the regular PyRETIS installation.

"""
import mdtraj
from pyretis.orderparameter import OrderParameter


def some_function_to_obtain_order(_):
    """A always right function test."""
    return 42


class PlaneDistanceX(OrderParameter):
    """A positional order parameter.

    This class defines a very simple order parameter which is
    the distance from a plane for a given particle.
    """

    def __init__(self, index, plane_position):
        """Initialise the order parameter.

        Parameters
        ----------
        index : integer
            Selects the particle to use.

        plane_position : float
            The location of the plane, along the x-axis.
        """
        txt = f'Distance from particle {index} to a plane at {plane_position}'
        super().__init__(description=txt)
        self.index = index
        self.plane_position = plane_position

    def calculate(self, system):
        """Calculate the order parameter."""
        filename = system.particles.config[0]
        traj = mdtraj.load(filename)
        return some_function_to_obtain_order(traj)
