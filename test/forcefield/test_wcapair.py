# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A test of the WCA pair potential."""
import logging
import pytest
import numpy as np
from pyretis.core.box import create_box
from pyretis.core.particles import Particles
from pyretis.core.system import System
from pyretis.forcefield.potentials.pairpotentials.wca import (
    DoubleWellWCA,
)
logging.disable(logging.CRITICAL)


CORRECT_FORCE = np.array([-2416.67712532,
                          -2416.67712532,
                          -2416.67712532])
CORRECT_VIRIAL = np.ones((3, 3)) * 483.33542506
CORRECT_VPOT = 873.18173451


def create_system_box(in_max=True):
    """Create a simple test system."""
    box = create_box(cell=[10, 10, 10])
    system = System(box=box)
    system.particles = Particles(system.get_dim())
    system.add_particle(name='Ar1', pos=np.array([1.0, 1.0, 1.0]),
                        mass=1.0, ptype=0)
    if in_max:
        system.add_particle(name='Ar2', pos=np.array([2.25, 1.0, 1.0]),
                            mass=1.0, ptype=0)
    else:
        system.add_particle(name='Ar2', pos=np.array([1.2, 1.2, 1.2]),
                            mass=1.0, ptype=0)
    system.add_particle(name='X', pos=np.array([1.0, 2.0, 1.0]),
                        mass=1.0, ptype=1)
    system.add_particle(name='X', pos=np.array([1.0, 3.25, 1.0]),
                        mass=1.0, ptype=1)
    pot = DoubleWellWCA()
    parameters = {
        'types': {(0, 0)},
        'rzero': 1.0,
        'width': 0.25,
        'height': 6.0,
    }
    pot.set_parameters(parameters)
    return system, box, pot


class TestWCAPairPotential:
    """Run tests for the WCA pair potential."""

    def test_parameters(self, caplog):
        """Test setup with parameters."""
        pot = DoubleWellWCA()
        parameters = {
            'types': {(0, 0)},
            'rzero': 1.0,
            'width': 0.25,
            'height': 6.0,
        }
        pot.set_parameters(parameters)
        for key, val in parameters.items():
            assert val == pot.params[key]
        # Check derived parameters:
        assert parameters['width']**2 == pot.params['width2']
        assert (parameters['rzero'] + parameters['width'] ==
                pot.params['rwidth'])
        assert parameters['height'] * 4.0 == pot.params['height4']
        # Check that extra parameters are ignored:
        parameters['junk'] = 'should-not-be-added!'
        logging.disable(logging.INFO)
        mod = 'pyretis.forcefield.potentials.pairpotentials.wca'
        with caplog.at_level(logging.WARNING, logger=mod):
            pot.set_parameters(parameters)
        logging.disable(logging.CRITICAL)
        assert 'junk' not in pot.params

    def test_activate(self):
        """Test that the WCA potential is activated for correct pairs."""
        pot = DoubleWellWCA()
        parameters = {
            'types': {(0, 0)},
            'rzero': 1.0,
            'width': 0.25,
            'height': 6.0,
        }
        pot.set_parameters(parameters)
        # pylint: disable=protected-access
        assert pot._activate(0, 0)
        assert not pot._activate(0, 1)
        parameters['types'] = {(0, 0), (0, 1)}
        pot.set_parameters(parameters)
        assert pot._activate(0, 0)
        assert pot._activate(1, 0)
        assert pot._activate(0, 1)
        assert not pot._activate(1, 1)
        parameters['types'] = None
        pot.set_parameters(parameters)
        assert pot._activate(0, 0)
        assert pot._activate(1, 0)
        assert pot._activate(0, 1)
        assert pot._activate(1, 1)

    def test_min_max(self):
        """Test that the min/max is correctly set for the WCA potential."""
        pot = DoubleWellWCA()
        parameters = {
            'types': {(0, 0)},
            'rzero': 1.0,
            'width': 0.25,
            'height': 6.0,
        }
        pot.set_parameters(parameters)
        min1, min2, max1 = pot.min_max()
        assert min1 == 1.0
        assert min2 == 1.5
        assert max1 == 1.25

    def test_potential(self):
        """Test that the potential is evaluated correctly."""
        system, _, pot = create_system_box()
        vpot = pot.potential(system)
        assert vpot == 6.0

    def test_force(self):
        """Test that the potential and virial is correctly evaluated."""
        system, _, pot = create_system_box(in_max=False)
        force, virial = pot.force(system)
        assert np.allclose(force[0], CORRECT_FORCE)
        assert np.allclose(force[1], -1.0 * CORRECT_FORCE)
        assert np.allclose(force[2], np.zeros(3))
        assert np.allclose(force[3], np.zeros(3))
        assert np.allclose(virial, CORRECT_VIRIAL)

    def test_potential_and_force(self):
        """Test that the potential and virial is correctly evaluated."""
        system, _, pot = create_system_box(in_max=False)
        vpot, force, virial = pot.potential_and_force(system)
        assert np.allclose(force[0], CORRECT_FORCE)
        assert np.allclose(force[1], -1.0 * CORRECT_FORCE)
        assert np.allclose(force[2], np.zeros(3))
        assert np.allclose(force[3], np.zeros(3))
        assert np.allclose(virial, CORRECT_VIRIAL)
        assert vpot == pytest.approx(CORRECT_VPOT)
