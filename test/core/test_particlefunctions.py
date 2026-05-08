# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the particle functions from pyretis.core.particlefunctions"""
import logging
import pytest
import numpy as np
from pyretis.core.box import create_box
from pyretis.core.system import System
from pyretis.core.particles import Particles
from pyretis.core.particlefunctions import (
    _get_vel_mass,
    atomic_kinetic_energy_tensor,
    calculate_kinetic_energy,
    calculate_kinetic_energy_tensor,
    calculate_kinetic_temperature,
    kinetic_temperature,
    calculate_linear_momentum,
    calculate_pressure_from_temp,
    calculate_pressure_tensor,
    calculate_scalar_pressure,
    calculate_thermo,
    calculate_thermo_path,
    reset_momentum,
)


def empty_particles(masses):
    """Set up a particle object for testing."""
    particles = Particles(dim=3)
    for i in masses:
        particles.add_particle(np.zeros(3), np.zeros(3), np.zeros(3), mass=i)
    particles.virial = np.zeros((3, 3))
    return particles


def test_get_selection():
    """Test the creation a particle list."""
    particles = empty_particles([1, 2, 1, 2])
    _, mass = _get_vel_mass(particles, selection=[0, 1])
    for i, j in zip(mass, (1, 2)):
        assert i[0] == pytest.approx(j)
    _, mass = _get_vel_mass(particles, selection=[0, 2])
    for i in mass:
        assert i[0] == pytest.approx(1.0)
    _, mass = _get_vel_mass(particles, selection=[1, 3])
    for i in mass:
        assert i[0] == pytest.approx(2.0)


def test_get_atomic_ek():
    """Test the atomic ek tensor calculation."""
    particles = empty_particles([1, 2, 1, 2])
    particles.vel = np.array([[1., 1., 1.], [2., 2., 2.],
                              [3., 3., 3.], [4., 4., 4.]])
    kin = atomic_kinetic_energy_tensor(particles)
    kin1 = atomic_kinetic_energy_tensor(particles, selection=[0, 2])
    kin2 = atomic_kinetic_energy_tensor(particles, selection=[1, 3])
    assert np.allclose(kin[0], kin1[0])
    assert np.allclose(kin[2], kin1[1])
    assert np.allclose(kin[1], kin2[0])
    assert np.allclose(kin[3], kin2[1])
    kin3 = atomic_kinetic_energy_tensor(particles, selection=[3])
    assert np.allclose(kin[3], kin3[0])


def test_calculate_ek():
    """Test calculation of kinetic energy."""
    particles = empty_particles([1, 2, 1, 2])
    particles.vel = np.array([[1., 1., 1.], [2.1, 2.1, 2.1],
                              [3.1, 3.1, 3.1], [4.1, 4.1, 4.1]])
    ekin, tensor = calculate_kinetic_energy(particles)
    tensor2 = calculate_kinetic_energy_tensor(particles)
    assert np.allclose(tensor, tensor2)
    assert ekin == pytest.approx(79.575)


def test_calculate_temp():
    """Test calculation of temperature."""
    particles = empty_particles([1, 2, 1, 0.01])
    particles.vel = np.ones_like(particles.pos)
    temp1, temp2, _ = calculate_kinetic_temperature(particles, 1.0)
    assert temp2 == pytest.approx(1.0025)
    for i in temp1:
        assert i == pytest.approx(1.0025)


def test_kinetic_temp():
    """Test calculation of kinetic temperature."""
    particles = empty_particles([1, 2, 1, 0.01])
    particles.vel = np.ones_like(particles.pos)
    vel = particles.vel
    mass = particles.mass
    temp1, temp2, _ = kinetic_temperature(vel, mass, 1.0)
    assert temp2 == pytest.approx(1.0025)
    for i in temp1:
        assert i == pytest.approx(1.0025)
    temp1, temp2, _ = kinetic_temperature(vel, mass, 1.0,
                                          dof=[1.0, 0.0, 0.0])
    assert temp2 == pytest.approx(1.11388888889)
    for i, j in zip(temp1, (1.33666667, 1.0025, 1.0025)):
        assert i == pytest.approx(j)


def test_linear_momentum():
    """Test calculation of linear momentum."""
    particles = empty_particles([1, 0.5, 2, 1])
    particles.vel = np.array([[1., 1., 1.], [2., 2., 2.],
                              [0.5, 0.5, 0.5], [1., 1., 1.]])
    mom = calculate_linear_momentum(particles)
    for i in mom:
        assert i == pytest.approx(4.)


def test_pressure_from_temp():
    """Test calculation of pressure."""
    particles = empty_particles([1, 0.5, 2, 1])
    particles.vel = np.array([[1., 1., 1.], [2., 2., 2.],
                              [0.5, 0.5, 0.5], [1., 1., 1.]])
    pvol, press = calculate_pressure_from_temp(particles, 3, 1.0, 1.0)
    assert pvol == pytest.approx(4.5)
    assert press == pytest.approx(4.5)
    pvol, press = calculate_pressure_from_temp(particles, 3, 1.0, 1.0,
                                               dof=[1., 0., 0.])
    assert pvol == pytest.approx(4.583333333)
    assert press == pytest.approx(4.583333333)


def test_pressure_tensors():
    """Test calculation of the pressure tensor and scalar."""
    particles = empty_particles([1, 2, 3, 4])
    particles.vel = np.array([[1., 1., 1.], [2., 2., 2.],
                              [0.5, 0.5, 0.5], [1., 1., 1.]])
    press = calculate_pressure_tensor(particles, 1.0)
    for i in press.ravel():
        assert i == pytest.approx(13.75)
    presss = calculate_scalar_pressure(particles, 1.0, 3.0)
    assert presss == pytest.approx(13.75)


def test_calculate_thermo():
    """Test the calculate_thermo method."""
    system = System(units='reduced', box=create_box(cell=[1., 1., 1.]))
    particles = empty_particles([1, 2, 3, 4])
    particles.vel = np.array([[1., 1., 1.], [2., 2., 2.],
                              [0.5, 0.5, 0.5], [1., 1., 1.]])
    particles.vpot = 123.456
    system.particles = particles
    res = calculate_thermo(system)
    correct = {
        'etot': 36.020250000000004,
        'vpot': 30.864,
        'press': 13.75,
        'mom': np.array([10.5, 10.5, 10.5]),
        'temp': 4.583333333333333,
        'ekin': 5.15625,
        'press-tens': 13.75 * np.ones((3, 3)),
    }
    for key in ('etot', 'vpot', 'press', 'temp', 'ekin'):
        assert res[key] == pytest.approx(correct[key])
    for key in ('mom', 'press-tens'):
        assert np.allclose(correct[key], res[key])


def test_thermo_path():
    """Test thermo calculation for path-style output."""
    system = System(units='reduced', box=create_box(cell=[1., 1.5, 1.]))
    particles = empty_particles([1, 2])
    particles.vel = np.array([[1., 1., 1.], [2., 2., 2.]])
    particles.vpot = 314.21
    system.particles = particles
    res = calculate_thermo_path(system)
    correct = {
        'temp': 9.0,
        'vpot': 314.21,
        'ekin': 13.5,
        'etot': 327.70999999999998
    }
    for key, val in correct.items():
        assert res[key] == pytest.approx(val)


def test_reset_mom():
    """Test that we can reset momentum."""
    particles = empty_particles([1, 0.5, 2, 1])
    particles.vel = np.array([[1., 1., 1.], [2., 2., 2.],
                              [0.5, 0.5, 0.5], [1., 1., 1.]])
    vel = np.copy(particles.vel)
    reset_momentum(particles)
    mom = calculate_linear_momentum(particles)
    assert np.allclose(mom, np.zeros(3))
    particles.vel = vel
    reset_momentum(particles, dim=[False, True, False])
    mom = calculate_linear_momentum(particles)
    correct = [4., 0.0, 4.0]
    for i, j in zip(mom, correct):
        assert i == pytest.approx(j)
