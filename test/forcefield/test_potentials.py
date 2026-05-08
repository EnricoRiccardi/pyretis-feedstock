# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A test for the simple potentials."""
import logging
import pytest
import numpy as np
from pyretis.core.system import System
from pyretis.core.particles import Particles
from pyretis.forcefield.potentials import DoubleWell, RectangularWell
logging.disable(logging.CRITICAL)


def evaluate_potential(positions, system, function):
    """Just evaluate the given potential"""
    all_pot, all_pot2 = [], []
    all_force, all_force2 = [], []
    all_vir, all_vir2 = [], []
    for i in positions:
        system.particles.pos = i
        pot = function.potential(system)
        all_pot.append(pot)
        force, vir = function.force(system)
        all_force.append(force)
        all_vir.append(vir)
        pot2, force2, vir2 = function.potential_and_force(system)
        all_pot2.append(pot2)
        all_force2.append(force2)
        all_vir2.append(vir2)
    out = []
    for i in (all_pot, all_pot2, all_force, all_force2, all_vir, all_vir2):
        out.append(np.array(i))
    return out


class TestDoubleWell:
    """Test the double well potential."""

    def test_doublewell_potential(self):
        """Test the double well potential."""
        params = {'a': 1.0, 'b': 2.0, 'c': 3.0}
        function = DoubleWell(**params)
        system = System()
        system.particles = Particles(dim=system.get_dim())
        system.add_particle(0.0)
        pos2 = np.linspace(-10, 10, 100)
        pos = np.reshape(np.linspace(-10, 10, 100), (100, 1, 1))
        out = evaluate_potential(pos, system, function)
        all_pot2 = (params['a'] * pos2**4 -
                    params['b'] * (pos2 - params['c'])**2)
        assert np.allclose(out[0], out[1])
        assert np.allclose(out[2], out[3])
        assert np.allclose(out[4], out[5])
        assert np.allclose(out[0], all_pot2)
        all_force2 = (-4.0 * params['a'] * pos**3 +
                      2.0 * params['b'] * (pos - params['c']))
        assert np.allclose(out[2], all_force2)

    def test_more_particles(self):
        """Test what happens if we have more particles."""
        params = {'a': 1.0, 'b': 2.0, 'c': 3.0}
        npart = 10
        function = DoubleWell(**params)
        system = System()
        system.particles = Particles(dim=system.get_dim())
        for _ in range(npart):
            system.add_particle(0.0)
        pos = np.linspace(-10, 10, 100)
        all_pot = []
        for i in pos:
            for j in range(npart):
                system.particles.pos[j][0] = i
            pot = function.potential(system)
            all_pot.append(pot)
        all_pot = np.array(all_pot)
        all_pot2 = params['a'] * pos**4 - params['b'] * (pos - params['c'])**2
        assert np.allclose(all_pot, all_pot2 * npart)


class TestRectangularWell:
    """Test the rectangular well potential."""

    def test_initiation(self, caplog):
        """Test the way we give input."""
        params = {'left': -1.0, 'right': 1.0}
        RectangularWell(**params)
        params = {'left': 1.0, 'right': -1.0}
        logging.disable(logging.INFO)
        with caplog.at_level(
                logging.WARNING,
                logger='pyretis.forcefield.potentials.potentials'
        ):
            RectangularWell(**params)
        logging.disable(logging.CRITICAL)

    def test_rectwell_potential(self):
        """Test the double well potential."""
        params = {'left': -1.0, 'right': 1.0}
        function = RectangularWell(**params)
        system = System()
        system.particles = Particles(dim=system.get_dim())
        system.add_particle(0.0)
        pos = np.linspace(-10, 10, 100)
        for i in pos:
            system.particles.pos[0][0] = i
            pot = function.potential(system)
            if params['left'] < i < params['right']:
                assert pot == 0.0
            else:
                assert pot == float('inf')
        for i in (params['left'], params['right']):
            system.particles.pos[0][0] = i
            pot = function.potential(system)
            assert pot == float('inf')
