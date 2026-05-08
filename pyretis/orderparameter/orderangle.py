# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This file contains classes to represent angle order parameters.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Angle (:py:class:`.Angle`)
    An angle defined by three atoms.
"""
import logging
from numpy import dot
from numpy import arccos  # pylint: disable=no-name-in-module
from numpy import clip
from numpy.linalg import norm
from pyretis.orderparameter.orderparameter import OrderParameter
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = ['Angle']


class Angle(OrderParameter):
    """An angle order parameter.

    This class defines an order parameter which is an angle
    ABC for 3 particles A, B and C. The angle is defined as the
    angle given by the two vectors BA and BC.

    Attributes
    ----------
    index : list of integers
        These are the indices of atoms to be used for the angle,
        i.e. system.particles.pos[index] will be used.
    periodic : boolean
        This determines if periodic boundaries should be applied to
        the positions/distances.

    """

    def __init__(self, index, periodic=False):
        """Initialise the order parameter.

        Parameters
        ----------
        index : list/tuple of integers
            The indices for the atoms defining the angle.
        periodic : boolean, optional
            This determines if periodic boundary conditions should be
            applied to the distance vectors.

        """
        try:
            if len(index) != 3:
                msg = ('Wrong number of atoms for angle definition. '
                       f'Expected 3 got {len(index)}')
                logger.error(msg)
                raise ValueError(msg)
        except TypeError as err:
            msg = 'Angle should be defined as a tuple/list of integers!'
            logger.error(msg)
            raise TypeError(msg) from err
        txt = f'Angle between particles {index[0]}, {index[1]} and {index[2]}'
        super().__init__(description=txt)
        self.periodic = periodic
        self.index = [int(i) for i in index]

    def calculate(self, system):
        """Calculate the angle.

        Parameters
        ----------
        system : object like :py:class:`.System`
            Object containing the positions and box info we use for
            the calculation.

        Returns
        -------
        out : list of floats
            The order parameters.

        """
        pos = system.particles.pos
        vector_ba = pos[self.index[1]] - pos[self.index[0]]
        vector_bc = pos[self.index[1]] - pos[self.index[2]]
        if self.periodic:
            vector_ba = system.box.pbc_dist_coordinate(vector_ba)
            vector_bc = system.box.pbc_dist_coordinate(vector_bc)
        vector_ba /= norm(vector_ba)
        vector_bc /= norm(vector_bc)
        angle = arccos(clip(dot(vector_ba, vector_bc), -1, 1))
        return [angle]
