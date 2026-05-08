# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the C implementation of the WCA potential.

This test is comparing the C implementation with the PyRETIS
implementation.
"""
# pylint: disable=invalid-name
import unittest
import itertools
import numpy as np
from pyretis.core import Particles, create_box, System
from pyretis.core.units import create_conversion_factors
from pyretis.forcefield import ForceField
from pyretis.forcefield.potentials import DoubleWellWCA
from pyretis.forcefield.potentials import PairLennardJonesCutnp
from pyretis.tools import generate_lattice
from wcafunctions import WCAPotential


def set_up_initial_state():
    """Create particles for the test."""
    create_conversion_factors('lj')
    lattice, size = generate_lattice('sq', [3, 3], density=0.6)
    lattice += np.random.random((lattice.shape)) * 0.05
    size = np.array(size)
    box = create_box(low=size[:, 0], high=size[:, 1], periodic=[True, True])
    system = System(temperature=1.0, units='lj', box=box)
    system.particles = Particles(dim=2)
    for i, pos in enumerate(lattice):
        if i < 2:
            system.add_particle(name='Ar', pos=pos, mass=1.0, ptype=1)
        else:
            system.add_particle(name='Ar', pos=pos, mass=1.0, ptype=0)
    msg = 'Created lattice with {} atoms.'
    print(msg.format(system.particles.npart))
    return system


def set_up_python_forcefield():
    """Create internal Python WCA force field."""
    wca_pot = PairLennardJonesCutnp(dim=2, shift=True, mixing='geometric')
    wca_params = {0: {'sigma': 1.0, 'epsilon': 1.0, 'factor': 2.**(1./6.)},
                  1: {'sigma': 1.0, 'epsilon': 1.0, 'factor': 2.**(1./6.)}}
    dwca_pot = DoubleWellWCA(dim=2)
    dwca_params = {'types': [(1, 1)],
                   'rzero': 1.0 * (2.0**(1.0/6.0)),
                   'height': 15.0, 'width': 0.5}
    forcefield = ForceField('Double well force field, PyRETIS',
                            potential=[wca_pot, dwca_pot],
                            params=[wca_params, dwca_params])
    return forcefield


def set_up_c_forcefield():
    """Create external C WCA force field."""
    wca_pot_c = WCAPotential()
    wca_paramsc = {'sigma': 1.0, 'epsilon': 1.0, 'rcut': 2.**(1./6.),
                   'idxi': 0, 'idxj': 1, 'rzero': 2.0**(1./6),
                   'height': 15.0, 'width': 0.5}
    forcefield = ForceField('Double well force field, PyRETIS + c',
                            potential=[wca_pot_c],
                            params=[wca_paramsc])
    return forcefield


class WCATest(unittest.TestCase):
    """Run the tests for the C potential class."""

    def test_wca_pot(self):
        """Test evaluation of potential."""
        system = set_up_initial_state()
        forcefield_python = set_up_python_forcefield()
        forcefield_c = set_up_c_forcefield()

        pot, force_virial, force_and_pot = [], [], []
        for ffield in [forcefield_python, forcefield_c]:
            system.forcefield = ffield
            vpot1 = system.potential()
            force1, virial1 = system.force()
            vpot2, force2, virial2 = system.potential_and_force()

            self.assertAlmostEqual(vpot1, vpot2)

            forceok = np.allclose(force1, force2)
            self.assertTrue(forceok)

            virialok = np.allclose(virial1, virial2)
            self.assertTrue(virialok)

            pot.append(vpot1)
            force_virial.append((force1, virial1))
            force_and_pot.append((vpot2, force2, virial2))
        for pair in itertools.combinations(range(len(pot)), 2):
            i, j = pair
            self.assertAlmostEqual(pot[i], pot[j])
            virialok = np.allclose(force_virial[i][1], force_virial[j][1])
            self.assertTrue(virialok)
            forceok = np.allclose(force_virial[i][0], force_virial[j][0])
            self.assertTrue(forceok)


if __name__ == '__main__':
    unittest.main()
