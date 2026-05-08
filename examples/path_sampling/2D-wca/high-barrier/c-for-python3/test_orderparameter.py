"""Test the calculation of the order parameter.

This test is checking that the order parameter is calculated
correctly.

1) The order parameter from PyRETIS.

2) A Python implementation.

3) A C implementation.
"""
# pylint: disable=invalid-name
import unittest
import numpy as np
from pyretis.core import Particles, create_box, System
from pyretis.core.units import create_conversion_factors
from pyretis.orderparameter import DistanceVelocity
from wcafunctions import WCAOrderParameter


class WCAOrderTest(unittest.TestCase):
    """Run the tests for the C potential class."""

    def test_wca_orderp(self):
        """Test evaluation of the order parameter."""
        box = create_box(low=[0., 0.], high=[3., 3.])
        particles = Particles(dim=2)
        particles.add_particle(np.array([1.0, 1.0]), np.zeros(2), np.zeros(2),
                               mass=1.0, name='A', ptype=0)
        particles.add_particle(np.array([1.0, 2.0]), np.zeros(2), np.zeros(2),
                               mass=1.0, name='A', ptype=0)
        create_conversion_factors('lj')
        system = System(box=box, units='lj')
        system.particles = particles

        order1 = WCAOrderParameter((0, 1))
        order2 = DistanceVelocity((0, 1), periodic=True)

        for i in np.arange(0.001, 5.0, 0.1):
            particles.pos[1] = np.array([1.0, i])
            particles.vel = np.random.random(particles.vel.shape)
            out1 = order1.calculate(system)
            out2 = order2.calculate(system)
            self.assertEqual(len(out1), len(out2))
            for orderi, orderj in zip(out1, out2):
                self.assertAlmostEqual(orderi, orderj)


if __name__ == '__main__':
    unittest.main()
