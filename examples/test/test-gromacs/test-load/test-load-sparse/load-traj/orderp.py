# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This file defines the order parameter used for the GROMACS example."""
import logging
import numpy as np
from pyretis.orderparameter import OrderParameter
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


class Distance(OrderParameter):
    """Distance(OrderParameter).

    This class computes the distance between two atoms.

    Attributes
    ----------
    name : string
        A human-readable name for the order parameter.

    """

    def __init__(self, index):
        """Set up the order parameter.

        Parameters
        ----------
        index : tuple of ints
            Indices of the two atoms whose distance is computed.

        """
        super().__init__(description='Water molecule distance')
        self.idx1 = index[0]
        self.idx2 = index[1]

    def calculate(self, system):
        """Calculate the order parameter.

        Here, the order parameter is just the distance between two
        particles.

        Parameters
        ----------
        system : object like `System` from `pyretis.core.system`
            Provides ``system.particles.pos`` with current coordinates.

        Returns
        -------
        out : list of floats
            The order parameter for the current frame only.

        """
        pos = system.particles.pos
        box = system.box.length
        dr = pos[self.idx2] - pos[self.idx1]
        if box is not None:
            dr -= np.rint(dr / box) * box
        return [float(np.linalg.norm(dr))]
