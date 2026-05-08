# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This file tests the Python implementation of the order parameter."""
import logging
import unittest
from pyretis.engines.gromacs import read_gromos96_file
from pyretis.core import System, Particles
from orderp import RingDiffusion
from ordermod import ordermod
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


class TestOrder(unittest.TestCase):
    """Just test that we can reproduce the FORTRAN implementation."""

    def test_orderp(self):
        """Test the FORTRAN implemented order parameter."""
        testo = RingDiffusion()
        _, pos, _, _ = read_gromos96_file('gromacs_input/conf.g96')
        system = System()
        system.particles = Particles(dim=3)
        system.particles.pos = pos
        lambpy = testo.calculate(system)[0]
        lambf = ordermod.calculate(pos)
        self.assertAlmostEqual(lambpy, lambf)


def print_idx():
    """Just print out the indices."""
    ordr = RingDiffusion()
    print('Indexes for the two groups:')
    print('Idx1', ordr.idx1)
    print('Idx2', ordr.idx2)


if __name__ == '__main__':
    unittest.main()
