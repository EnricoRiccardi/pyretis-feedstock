# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""The order parameter for the hysteresis example."""
import logging
from pyretis.orderparameter.orderparameter import OrderParameter
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


class OrderX(OrderParameter):
    """A positional order parameter.

    Order parameter for the hysteresis example. In addition to using
    the position, we also use the energy to tell if we are in states A/B.


    Attributes
    ----------
    index : integer
        This is the index of the atom which will be used, i.e.
        system.particles.pos[index] will be used.
    inter_a : float
        An interface such that we are in state A for postions < inter_a.
    inter_b : float
        An interface such that we are in state B for postions > inter_b.
    energy_a : float
        An energy such that we are in state A for potential energy < energy_a.
    energy_b : float
        An energy such that we are in state A for potential energy < energy_b.
    dim : integer
        This is the dimension of the coordinate to use.
        0, 1 or 2 for 'x', 'y' or 'z'.
    periodic : boolean
        This determines if periodic boundaries should be applied to
        the position or not.

    """

    def __init__(self, index, inter_a, inter_b, energy_a, energy_b,
                 dim='x', periodic=False):
        """Initialise the order parameter.

        Parameters
        ----------
        index : int
            This is the index of the atom we will use the position of.
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
        dim : string
            This select what dimension we should consider,
            it should equal 'x', 'y' or 'z'.
        periodic : boolean, optional
            This determines if periodic boundary conditions should be
            applied to the position.

        """
        txt = f'Position of particle {index} (dim: {dim})'
        super().__init__(description=txt)
        self.inter_a = inter_a
        self.inter_b = inter_b
        self.energy_a = energy_a
        self.energy_b = energy_b
        self.periodic = periodic
        self.index = index
        dims = {'x': 0, 'y': 1, 'z': 2}
        try:
            self.dim = dims[dim]
        except KeyError:
            msg = f'Unknown dimension {dim} requested'
            logger.critical(msg)
            raise

    def calculate(self, system):
        """Calculate the order parameter.

        Parameters
        ----------
        system : object like :py:class:`.System`
            This object is used for the actual calculation, typically
            only `system.particles.pos` and/or `system.particles.vel`
            will be used. In some cases `system.forcefield` can also be
            used to include specific energies for the order parameter.

        Returns
        -------
        out : float
            The order parameter.

        """
        particles = system.particles
        pos = particles.pos[self.index]
        lmb = pos[self.dim]
        if lmb < self.inter_a:
            if particles.vpot > self.energy_a:
                lmb = self.inter_a
        elif lmb > self.inter_b:
            if particles.vpot > self.energy_b:
                lmb = self.inter_b
        if self.periodic:
            return [system.box.pbc_coordinate_dim(lmb, self.dim)]
        return [lmb]
