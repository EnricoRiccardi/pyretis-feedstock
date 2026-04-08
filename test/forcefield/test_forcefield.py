# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the ForceField and PotentialFunction classes."""
import logging
import pytest
import numpy as np
from pyretis.core.system import System
from pyretis.core.particles import Particles
from pyretis.forcefield import ForceField, PotentialFunction
logging.disable(logging.CRITICAL)


class PotentialTest(PotentialFunction):
    """A potential function to use in tests."""

    def __init__(self, desc='Test potential'):
        super().__init__(dim=1, desc=desc)
        self.params = {'a': 10}

    def potential(self, system):
        """Evaluate the potential."""
        pos = system.particles.pos
        return np.sum(pos * pos * self.params['a'])

    def force(self, system):
        """Evaluate force and virial."""
        pos = system.particles.pos
        return pos * self.params['a'], 2.0 * pos

    def potential_and_force(self, system):
        """Evaluate potential, force and virial."""
        pot = self.potential(system)
        force, virial = self.force(system)
        return pot, force, virial


class TestForceField:
    """Test set-up of force fields."""

    def test_forcefield_class(self, caplog):
        """Test functionality of the ForceField class."""
        system = System()
        system.particles = Particles(dim=system.get_dim())
        system.add_particle(1.0)
        forcefield = ForceField('Generic testing force field')
        param1 = {'a': 1.0}
        pot1 = PotentialTest()
        forcefield.add_potential(pot1, parameters=param1)

        force, virial = forcefield.evaluate_force(system)
        assert 1.0 == force
        assert 2.0 == virial

        vpot = forcefield.evaluate_potential(system)
        assert 1.0 == vpot

        vpot, force, virial = forcefield.evaluate_potential_and_force(system)
        assert 1.0 == force
        assert 2.0 == virial
        assert 1.0 == vpot

        param2 = {'a': 2.0}
        forcefield.update_potential_parameters(pot1, param2)

        vpot, force, virial = forcefield.evaluate_potential_and_force(system)
        assert 2.0 == force
        assert 2.0 == virial
        assert 2.0 == vpot

        potr, paramr = forcefield.remove_potential(pot1)
        assert pot1 is potr
        assert param2 is paramr

        pot2 = PotentialTest()
        potr, paramr = forcefield.remove_potential(pot2)
        assert potr is None
        assert paramr is None
        logging.disable(logging.INFO)
        with caplog.at_level(logging.WARNING,
                             logger='pyretis.forcefield.forcefield'):
            forcefield.update_potential_parameters(pot2, param2)
        logging.disable(logging.CRITICAL)

    def test_add_potential(self):
        """Test the forcefield add potential method in more detail."""
        forcefield = ForceField('Generic testing force field')
        param1 = {'a': 1.0}
        pot1 = PotentialTest()
        ret = forcefield.add_potential(pot1, parameters=param1)
        assert ret
        ret = forcefield.add_potential(None, parameters=param1)
        assert not ret
        pot2 = PotentialTest()
        forcefield = ForceField('Generic testing force field',
                                potential=[pot1, pot2],
                                params=[param1])
        assert pot1 in forcefield.potential
        assert pot2 in forcefield.potential
        assert param1 == forcefield.params[0]
        assert forcefield.params[1] is None

    def test_evaluation(self):
        """Test evaluation of the force field."""
        system = System()
        system.particles = Particles(dim=system.get_dim())
        system.add_particle(1.0)
        param1 = {'a': 1.0}
        pot1 = PotentialTest()
        param2 = {'a': 2.0}
        pot2 = PotentialTest()
        forcefield = ForceField('Generic testing force field',
                                potential=[pot1, pot2],
                                params=[param1, param2])
        vpot = forcefield.evaluate_potential(system)
        assert vpot == 3.0
        _, force, virial = forcefield.evaluate_potential_and_force(system)
        assert np.allclose(force, np.array([[3.0]]))
        assert np.allclose(virial, np.array([[4.0]]))
        force, virial = forcefield.evaluate_force(system)
        assert np.allclose(force, np.array([[3.0]]))
        assert np.allclose(virial, np.array([[4.0]]))

    def test_print_potentials(self):
        """Test the printing of force field information."""
        param1 = {'a': 1.0}
        pot1 = PotentialTest()
        param2 = {'a': 2.0}
        pot2 = PotentialTest()
        forcefield = ForceField('Generic testing force field',
                                potential=[pot1, pot2],
                                params=[param1, param2])
        txt = forcefield.print_potentials()
        correct_txt = ['Force field: Generic testing force field',
                       '1: Test potential', '2: Test potential']
        for i, j in zip(txt.split('\n'), correct_txt):
            assert i.strip() == j.strip()
