# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This file defines the order parameter used for the GROMACS example."""
import os
import logging
import mdtraj as md
import numpy as np
from pyretis.inout.formats.gromacs import read_gromos96_file
from pyretis.orderparameter import OrderParameter
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


class Distance(OrderParameter):
    """Distance(OrderParameter).

    This class computes the distance between the O of water using mdtraj.

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
            This is the indices of the atom we will use the position of.

        """
        super().__init__(description='Water molecule distance')
        self.idx1 = index[0]
        self.idx2 = index[1]
        self.top = 'gromacs_input/conf.gro'
        self.top_path = os.path.join(
            os.path.dirname(__file__), self.top
        )
        self.topology = md.load(self.top_path).topology

    def calculate(self, system):
        """Calculate the order parameter.

        Here, the order parameter is just the distance between two
        particles.

        Parameters
        ----------
        system : object like `System` from `pyretis.core.system`
            This object is used to import the file names required
            for mdtraj.

        Returns
        -------
        out : list of floats
            The order parameter for the current frame only.

        """
        filename, frame_index = system.particles.config
        atom_pair = [(self.idx1, self.idx2)]
        if frame_index is None:
            frame_index = 0

        extension = os.path.splitext(filename)[1]
        if extension in ('.trr', '.xtc'):
            trj = md.load_frame(filename, frame_index, top=self.top_path)
        elif extension == '.g96':
            _, xyz, _, box = read_gromos96_file(filename)
            kwargs = {}
            if box is not None:
                kwargs['unitcell_lengths'] = np.asarray(
                    box, dtype=float
                )[None, :]
                kwargs['unitcell_angles'] = np.array([[90.0, 90.0, 90.0]])
            trj = md.Trajectory(
                np.asarray(xyz, dtype=float)[None, :, :],
                self.topology,
                **kwargs,
            )
        else:
            trj = md.load(filename, top=self.top_path)

        orderp = md.compute_distances(trj, atom_pair, periodic=True)
        return orderp[0].tolist()
