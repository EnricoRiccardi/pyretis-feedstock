# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This file defines the order parameter used for the WCA example."""
import logging
import numpy as np
from pyretis.orderparameter import OrderParameter
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


class OrderParameterWCAJCP1(OrderParameter):
    """OrderParameterWCAJCP1(OrderParameter).

    This class defines a very simple order parameter which is just
    the scalar distance between two particles.

    Attributes
    ----------
    index : tuple of integers
        These are the indices used for the two particles.
        `system.particles.pos[index[0]]` and
        `system.particles.pos[index[1]]` will be used.
    periodic : boolean
        This determines if periodic boundaries should be applied to
        the position or not.

    """

    def __init__(self, index, periodic=True):
        """Set uå the order parameter.

        Parameters
        ----------
        index : tuple of ints
            This is the indices of the atom we will use the position of.
        periodic : boolean, optional
            This determines if periodic boundary conditions should be
            applied to the position.

        """
        pbc = 'Periodic' if periodic else 'Non-periodic'
        txt = '{} distance particles {} and {}'.format(pbc, index[0],
                                                       index[1])
        super().__init__(description=txt)
        self.periodic = periodic
        self.index = index

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
        # pylint: disable=invalid-name
        particles = system.particles
        delta = particles.pos[self.index[1]] - particles.pos[self.index[0]]
        if self.periodic:
            delta = system.box.pbc_dist_coordinate(delta)
        r = np.sqrt(np.dot(delta, delta))
        dx = delta
        dv = particles.vel[self.index[1]] - particles.vel[self.index[0]]
        dxdv = np.dot(dx, dv) / r
        m1 = particles.mass[self.index[0]]
        m2 = particles.mass[self.index[1]]
        m = m1 * m2 / (m1 + m2)
        potential_func = system.forcefield.potential[0]
        if potential_func is None:
            return [float(r)]
        if r < 1.2:
            pot = potential_func.potential_well(system) + 0.5 * m * (dxdv)**2
            orderp = 1.19
            if pot < 1.5:
                orderp = 1.18 - (1.5 - pot) / 0.5 * 0.02
        elif r > 1.42:
            pot = potential_func.potential_well(system) + 0.5 * m * (dxdv)**2
            orderp = 1.43
            if pot < 5.0:
                orderp = 1.44 + (5.0 - pot) / 0.5 * 0.02
        else:
            orderp = r
        return [float(orderp), dxdv]
