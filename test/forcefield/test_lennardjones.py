# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A test of the Lennard-Jones pair potential."""
import logging
import pytest
import numpy as np
from pyretis.core.box import create_box
from pyretis.core.particles import Particles
from pyretis.core.system import System
from pyretis.forcefield.potentials.pairpotentials.lennardjones import (
    PairLennardJonesCut,
    PairLennardJonesCutnp,
)
logging.disable(logging.CRITICAL)


CORRECT_FORCE = np.array([[-393024., -2880., 0.],
                          [390144., -390144., 0.],
                          [2880., 393024., 0.]])
CORRECT_VIRIAL = np.array([[196512., 1440., 0.0],
                           [1440., 196512., 0.0],
                           [0.0, 0.0, 0.0]])
CORRECT_VPOT = 32480.0489507


def create_system_box(use_numpy=False):
    """Create a simple test system."""
    box = create_box(cell=[10, 10, 10])
    system = System(box=box)
    system.particles = Particles(system.get_dim())
    system.add_particle(name='Ar', pos=np.array([1.0, 1.0, 1.0]),
                        mass=1.0, ptype=0)
    system.add_particle(name='Ar', pos=np.array([1.5, 1.0, 1.0]),
                        mass=1.0, ptype=0)
    system.add_particle(name='Ar', pos=np.array([1.5, 1.5, 1.0]),
                        mass=1.0, ptype=0)
    if use_numpy:
        pot = PairLennardJonesCutnp()
    else:
        pot = PairLennardJonesCut()
    parameters = {
        0: {'sigma': 1, 'epsilon': 1, 'factor': 2.5},
        1: {'sigma': 2, 'epsilon': 2, 'factor': 2.5},
    }
    pot.set_parameters(parameters)
    return system, box, pot


class TestLennardJonesCutPairPotential:
    """Run tests for the Lennard-Jones pair potential."""

    def test_parameters(self):
        """Test setup with parameters."""
        pot = PairLennardJonesCut()
        parameters = {
            0: {'sigma': 1, 'epsilon': 1, 'factor': 2.5},
            1: {'sigma': 2, 'epsilon': 2, 'factor': 2.5},
        }
        pot.set_parameters(parameters)
        correct_pair = {
            (0, 0): {'epsilon': 1, 'rcut': 2.5, 'sigma': 1.},
            (0, 1): {'epsilon': np.sqrt(2.), 'rcut': 2.5 * np.sqrt(2.),
                     'sigma': np.sqrt(2.)},
            (1, 1): {'sigma': 2, 'epsilon': 2, 'rcut': 2.5 * 2.},
        }
        correct_pair[(1, 0)] = correct_pair[(0, 1)]
        for key, val in correct_pair.items():
            for key2, val2 in val.items():
                assert val2 == pytest.approx(pot.params[key][key2])
        parameters = {
            0: {'sigma': 1, 'epsilon': 1, 'factor': 0.0},
        }
        pot.set_parameters(parameters)
        assert pot.params[(0, 0)]['rcut'] == pytest.approx(0.0)

    def test_potential(self):
        """Test evaluation of the Lennard-Jones potential."""
        system, _, pot = create_system_box()
        vpot = pot.potential(system)
        assert vpot == pytest.approx(CORRECT_VPOT)

    def test_force(self):
        """Test evaluation of the Lennard-Jones force."""
        system, _, pot = create_system_box()
        force, virial = pot.force(system)
        assert np.allclose(force, CORRECT_FORCE)
        assert np.allclose(virial, CORRECT_VIRIAL)

    def test_potential_and_force(self):
        """Test evaluate of Lennard-Jones potential and force."""
        system, _, pot = create_system_box()
        vpot, force, virial = pot.potential_and_force(system)
        assert np.allclose(force, CORRECT_FORCE)
        assert np.allclose(virial, CORRECT_VIRIAL)
        assert vpot == pytest.approx(CORRECT_VPOT)


class TestLennardJonesCutPairPotentialnp:
    """Run tests for the Lennard-Jones pair potential (numpy)."""

    def test_potential(self):
        """Test evaluation of the Lennard-Jones potential (numpy)."""
        system, _, pot = create_system_box(use_numpy=True)
        vpot = pot.potential(system)
        assert vpot == pytest.approx(CORRECT_VPOT)

    def test_force(self):
        """Test evaluation of the Lennard-Jones force (numpy)."""
        system, _, pot = create_system_box(use_numpy=True)
        force, virial = pot.force(system)
        assert np.allclose(force, CORRECT_FORCE)
        assert np.allclose(virial, CORRECT_VIRIAL)

    def test_potential_and_force(self):
        """Test evaluate of Lennard-Jones potential and force (numpy)."""
        system, _, pot = create_system_box(use_numpy=True)
        vpot, force, virial = pot.potential_and_force(system)
        assert np.allclose(force, CORRECT_FORCE)
        assert np.allclose(virial, CORRECT_VIRIAL)
        assert vpot == pytest.approx(CORRECT_VPOT)
