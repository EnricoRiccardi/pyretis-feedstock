# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Definition of a dihedral order parameter.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dihedral (:py:class:`.Dihedral`)
    A dihedral angle order parameter.
"""
import logging
from numpy import dot, cross, arctan2
from numpy.linalg import norm
from pyretis.orderparameter.orderparameter import OrderParameter
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = ['Dihedral']


class Dihedral(OrderParameter):
    """Calculates the dihedral angle defined by 4 atoms.

    The angle definition is given by Blondel and Karplus,
    J. Comput. Chem., vol. 17, 1996, pp. 1132--1141. If we
    label the 4 atoms A, B, C and D, then the angle is given by
    the vectors u = A - B, v = B - C, w = D - C

    Attributes
    ----------
    index : list/tuple of integers
        These are the indices for the atoms to use in the
        definition of the dihedral angle.
    periodic : boolean
        This determines if periodic boundaries should be applied to
        the position or not.

    """

    def __init__(self, index, periodic=False):
        """Initialise the order parameter.

        Parameters
        ----------
        index : list/tuple of integers
            This list gives the indices for the atoms to use in the
            definition of the dihedral angle.
        periodic : boolean, optional
            This determines if periodic boundary conditions should be
            applied to the distance vectors.

        """
        try:
            if len(index) != 4:
                msg = ('Wrong number of atoms for dihedral definition. '
                       f'Expected 4 got {len(index)}')
                logger.error(msg)
                raise ValueError(msg)
        except TypeError as err:
            msg = 'Dihedral should be defined as a tuple/list of integers!'
            logger.error(msg)
            raise TypeError(msg) from err
        self.index = [int(i) for i in index]
        txt = ('Dihedral angle between particles '
               f'{index[0]}, {index[1]}, {index[2]} and {index[3]}')
        super().__init__(description=txt)
        self.periodic = periodic

    def calculate(self, system):
        """Calculate the dihedral angle.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The object containing the information we need to calculate
            the order parameter.

        Returns
        -------
        out : list of float
            The order parameter.

        """
        pos = system.particles.pos
        vector1 = pos[self.index[0]] - pos[self.index[1]]
        vector2 = pos[self.index[1]] - pos[self.index[2]]
        vector3 = pos[self.index[3]] - pos[self.index[2]]
        if self.periodic:
            vector1 = system.box.pbc_dist_coordinate(vector1)
            vector2 = system.box.pbc_dist_coordinate(vector2)
            vector3 = system.box.pbc_dist_coordinate(vector3)
        # Norm to simplify formulas:
        vector2 /= norm(vector2)
        denom = (dot(vector1, vector3) -
                 dot(vector1, vector2) * dot(vector2, vector3))
        numer = dot(cross(vector1, vector2), vector3)
        angle = arctan2(numer, denom)
        return [angle]
