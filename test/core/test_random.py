# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test functionality for the random generator classes."""
import logging
import pytest
import numpy as np
from numpy.random import RandomState
from pyretis.core.system import System
from pyretis.core.box import create_box
from pyretis.core.particles import Particles
from pyretis.core.particlefunctions import (
    calculate_linear_momentum,
    calculate_kinetic_temperature,
)
from pyretis.core.random_gen import (
    RandomGenerator,
    MockRandomGenerator,
    ReservoirSampler,
    RandomGeneratorBorg,
    MockRandomGeneratorBorg,
    create_random_generator,
)


def test_rand(caplog):
    """Test that we can draw random numbers in [0, 1)."""
    rgen = RandomGenerator(seed=0)
    for i in (0, 1, 10, 100, 10000):
        numbers = rgen.rand(shape=i)
        assert all([j >= 0 for j in numbers])
        assert all([j < 1 for j in numbers])
        assert len(numbers) == i

    # Without arguments
    number = rgen.rand()
    assert len(number) == 1
    # Test that it fails as we expect:
    with pytest.raises(TypeError):
        rgen.rand(1, 2)
    with pytest.raises(TypeError):
        rgen.rand((1, 2))


def test_random_integers(caplog):
    """Test generation for [low, high]."""
    rgen = RandomGenerator(seed=0)
    for i in (-5, 0, 10, 100):
        for j in (-5, 0, 10, 100):
            if i >= j + 1:
                with pytest.raises(ValueError):
                    rgen.random_integers(i, j)
            else:
                for _ in range(100):  # just repeat a bit
                    number = rgen.random_integers(i, j)
                    assert i <= number <= j
    # Just draw ones:
    numbers = [rgen.random_integers(1, 1) for _ in range(10)]
    assert all([i == 1 for i in numbers])
    # Just draw 1 or 2
    numbers = [rgen.random_integers(1, 2) for _ in range(100)]
    assert all([i in (1, 2) for i in numbers])


def test_random_normal(caplog):
    """Test generation of numbers from normal distribution."""
    rgen = RandomGenerator(seed=0)
    loc = 1.2345
    std = 0.2468
    numbers = rgen.normal(loc=loc, scale=std, size=10000)
    assert np.average(numbers) == pytest.approx(loc, abs=0.01)
    assert np.std(numbers) == pytest.approx(std, abs=0.01)


def test_random_normal_shape(caplog):
    """Test drawing of numbers from normal distribution with shape."""
    rgen = RandomGenerator(seed=0)
    shape = (10, 3)
    numbers = rgen.normal(loc=0.0, scale=2.0, size=shape)
    assert numbers.shape == shape

    # Pretend that we have 6 particles with different "mass"
    sigma = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0]
    tol = 0.1
    for dim in (1, 2, 3):
        # lets draw numbers a couple of times:
        pos = []
        for _ in range(500):
            numbers = rgen.normal(loc=0.0, scale=np.repeat(sigma, dim))
            numbers.shape = (len(sigma), dim)
            pos.append(numbers)
        pos = np.array(pos)
        # std over all drawn matrices:
        std = np.std(pos, axis=(0,))
        # compare for each dimension:
        for i in range(dim):
            std_diff = np.abs(std[:, i] - sigma) / sigma
            assert all([j < tol for j in std_diff])


def test_multivariate_normal(caplog):
    """Just test that we can draw from the multivariate distribution."""
    rgen = RandomGenerator(seed=0)
    mean = np.array([[1.0, 0.0], [0.0, 1.0]])
    cov = np.array([[1.0, 0.0], [0.0, 1.0]])
    numbers = rgen.multivariate_normal(mean, cov)
    assert numbers.shape == (1, 2, 2)
    numbers = rgen.multivariate_normal(mean, cov, size=2)
    assert numbers.shape == (2, 2, 2)


def test_draw_maxwellian_velocities(caplog):
    """Test that we can draw with the system object as input."""
    temperature = 2.0
    mass = np.array([1.0, 2.0, 4.0, 16.0, 256.0, 65536.0])
    sigv = np.sqrt(temperature / mass)
    tol = 0.1
    rgen = RandomGenerator(seed=0)
    for dim in (1, 2, 3):
        system = System(
            temperature=temperature,
            units='reduced',
            box=create_box(periodic=[False]*dim))
        system.particles = Particles(dim=dim)
        for i in mass:
            system.add_particle(np.zeros(dim), mass=i, name='Ar', ptype=0)
        all_vel = []
        for _ in range(1000):
            veli, _ = rgen.draw_maxwellian_velocities(system)
            all_vel.append(veli)
        vel = np.array(all_vel)
        # std over all drawn matrices:
        std = np.std(vel, axis=(0,))
        # compare for each dimension:
        for i in range(dim):
            std_diff = np.abs(std[:, i] - sigv) / sigv
            assert all([j < tol for j in std_diff])


def test_generate_maxwellian(caplog):
    """Test the generate_maxwellian_velocities method."""
    particles = Particles(dim=3)
    dof = [1., 1., 1.]
    heavy = []
    light = []
    for i in range(10):
        if i % 2 == 0:
            mass = 1
            light.append(i)
        else:
            mass = 2
            heavy.append(i)
        particles.add_particle(np.zeros(3), np.zeros(3),
                               np.zeros(3), mass=mass)
    rgen = RandomGenerator(seed=0)
    for mom in (True, False):
        rgen.generate_maxwellian_velocities(particles, 1.0, 2.0,
                                            dof=dof, momentum=mom)
        _, temp, _ = calculate_kinetic_temperature(particles, 1.0,
                                                   dof=dof)
        assert temp == pytest.approx(2.0)
        if mom:
            assert np.allclose(np.zeros(3),
                               calculate_linear_momentum(particles))

    rgen.generate_maxwellian_velocities(particles, 1.0, 2.0,
                                        dof=dof, momentum=True,
                                        selection=heavy)
    rgen.generate_maxwellian_velocities(particles, 1.0, 2.0,
                                        dof=None, momentum=False,
                                        selection=light)
    _, temp, _ = calculate_kinetic_temperature(particles, 1.0,
                                               dof=dof)
    assert temp == pytest.approx(2.0)
    assert np.allclose(
        np.zeros(3),
        calculate_linear_momentum(particles, selection=heavy)
    )
    assert not np.allclose(
        np.zeros(3),
        calculate_linear_momentum(particles, selection=light)
    )


def test_state(caplog):
    """Test that we can set and get the state of the generator."""
    rgen = RandomGenerator(seed=123)
    for _ in range(5):
        rgen.random_integers(1, 1000)
    state = rgen.get_state()
    numbers1 = [rgen.random_integers(1, 1000) for _ in range(10)]
    rgen.set_state(state)
    numbers2 = [rgen.random_integers(1, 1000) for _ in range(10)]
    assert numbers1 == numbers2


def test_create(caplog):
    """Test that we create random generators from settings."""
    settings = {}
    rgen = create_random_generator(settings)
    assert rgen.seed == 0
    assert isinstance(rgen, RandomGenerator)

    settings = {'seed': 100}
    rgen = create_random_generator(settings)
    assert rgen.seed == 100
    assert isinstance(rgen, RandomGenerator)

    settings = {'rgen': 'mock', 'seed': 101}
    rgen = create_random_generator(settings)
    assert rgen.seed == 101
    assert isinstance(rgen, MockRandomGenerator)


def test_reservoir_sampler_init(caplog):
    """Test the initialisation."""
    rgen = ReservoirSampler()
    assert isinstance(rgen.rgen, RandomState)
    rgen = ReservoirSampler(length=10, rgen=RandomGenerator(seed=1))
    assert isinstance(rgen.rgen, RandomGenerator)


def test_reservoir_sampler_add_get(caplog):
    """Test that we can add to the reservoir."""
    # Note: MockRandomGenerator used here.
    rgen = ReservoirSampler(length=10, rgen=MockRandomGenerator(seed=1))
    correct = [21, 0, 0, 0, 2, 2, 18, 29, 1, 0]
    for i in range(30):
        rgen.append(i)
        assert len(rgen.reservoir) == 10
    for i, j in zip(correct, rgen.reservoir):
        assert i == j
    for i in range(10):
        j = rgen.get_item()
        assert correct[i] == j

    with caplog.at_level(logging.CRITICAL):
        j = rgen.get_item()
    assert correct[0] == j
    assert "Out of bounds in the reservoir sampler!" in caplog.text


def test_mock_random_generator_state(caplog):
    """Test that we can get and set the state of the generator."""
    rgen = MockRandomGenerator(seed=987)
    for _ in range(5):
        rgen.random_integers(1, 100)
    state = rgen.get_state()
    numbers1 = [rgen.random_integers(1, 100) for _ in range(10)]
    rgen.set_state(state)
    numbers2 = [rgen.random_integers(1, 100) for _ in range(10)]
    assert numbers1 == numbers2


def test_mock_random_generator_rand(caplog):
    """Test that we can draw fake random numbers in [0, 1)."""
    rgen = MockRandomGenerator(seed=0)
    numbers = rgen.rand(shape=5)
    assert np.allclose(numbers, rgen.rgen[0:5])


def test_mock_random_generator_random_integers(caplog):
    """Test that we can draw fake random integers."""
    rgen = MockRandomGenerator(seed=0)
    correct = [14, 10, 14, 15, 13]
    for i in correct:
        j = rgen.random_integers(10, 15)
        assert i == j
    rgen = MockRandomGenerator(seed=0)
    for _ in range(len(rgen.rgen)):
        j = rgen.random_integers(4, 9)
        assert 4 <= j <= 9


def test_mock_random_generator_normal(caplog):
    """Test that we can draw fake normal numbers."""
    rgen = MockRandomGenerator(seed=0)
    numbers = rgen.normal()
    assert numbers[0] == pytest.approx(rgen.rgen[0])
    numbers = rgen.normal(loc=1.0, scale=10, size=5)
    for i, j in zip(numbers, rgen.rgen[1:6]):
        assert i == pytest.approx(j)


def test_mock_random_generator_multivariate_normal(caplog):
    """Test that we can draw fake normal numbers."""
    rgen = MockRandomGenerator(seed=0)
    correct = np.array([[0.0178008, 0.01044599]])
    numbers = rgen.multivariate_normal(1.0, None)
    assert np.allclose(correct, numbers)
    numbers = rgen.multivariate_normal([1.0, 1.0], None, size=3)
    correct = np.array([[0.01765968, 0.01976767],
                        [0.01537996, 0.01986571],
                        [0.01363436, 0.01553565]])
    assert np.allclose(correct, numbers)


def test_mock_random_generator_choice(caplog):
    a = [1, 2, 3, 4]
    rgen = MockRandomGenerator(seed=0)
    correct = [4, 1, 4, 4, 3, 4, 2]
    numbers = rgen.choice(a, 7)
    assert np.allclose(numbers, correct)


def test_mock_random_generator_choice_no_size(caplog):
    a = [1, 2, 3, 4]
    rgen = MockRandomGenerator(seed=0)
    correct = [4]
    numbers = rgen.choice(a)
    assert np.allclose(numbers, correct)


def test_mock_random_generator_choice_no_replace(caplog):
    a = [1, 2, 3, 4]
    rgen = MockRandomGenerator(seed=0)
    correct = [4, 1, 3, 2]
    numbers = rgen.choice(a, size=4, replace=False)
    assert np.allclose(numbers, correct)


def test_mock_random_generator_choice_altered_p(caplog):
    a = [1, 2, 3, 4]
    rgen = MockRandomGenerator(seed=0)
    correct = [1, 1, 1, 1]
    numbers = rgen.choice(a, size=4, p=[1, 0, 0, 0])
    assert np.allclose(numbers, correct)


def test_borg_mock(caplog):
    """Test that we share the Mock state."""
    rgens = [
        MockRandomGeneratorBorg(
            seed=i,
            norm_shift=(i == 0),
        ) for i in range(5)
    ]
    state0 = rgens[0].get_state()['state']
    # Are all states the same?
    assert all([i.get_state()['state'] == state0 for i in rgens])
    # Did the norm_shift get set to the first value?
    assert all([i.norm_shift for i in rgens])
    # Change the state of one member by requesting a number:
    _ = rgens[0].rand()
    state1 = rgens[0].get_state()['state']
    # Check that the state did change:
    assert state0 != state1
    # Are all states still the same?
    assert all([i.get_state()['state'] == state1 for i in rgens])
    # Set state manually
    rgens[-1].set_state({'state': 3})
    assert all([i.get_state()['state'] == 3 for i in rgens])
    MockRandomGeneratorBorg.reset_state()


def assert_equal_npstates(state1, state2):
    """Compare two random states coming from numpy.random.get_state()."""
    assert len(state1) == len(state2)
    assert state1[0] == 'MT19937'
    assert state1[0] == state2[0]
    assert np.allclose(state1[1], state2[1])
    assert state1[2] == state2[2]
    assert state1[3] == state2[3]
    assert np.allclose(state1[4], state2[4])


def test_borg_randomgen(caplog):
    """Test that we share the Mock state."""
    rgens = [
        RandomGeneratorBorg(seed=i) for i in range(5)
    ]
    # Check that rgen objects are identical:
    obj = rgens[0].rgen
    for i in rgens:
        assert obj is i.rgen
    # Check that the states are the same:
    state0 = rgens[0].get_state()['state']
    for i in rgens:
        assert_equal_npstates(state0, i.get_state()['state'])
    RandomGeneratorBorg.reset_state()


def test_new_swarm(caplog):
    """Test that we make a new swarm without altering the old one."""
    rgenA = RandomGeneratorBorg(seed=0)
    rgb2 = RandomGeneratorBorg.make_new_swarm()
    rgenB = rgb2(seed=0)
    rgenC = RandomGeneratorBorg(seed=0)
    rgenD = rgb2(seed=0)

    assert rgenA.rgen is not rgenB.rgen  # Test new swarm
    assert rgenA.rgen is rgenC.rgen  # Test that the old swarm still works
    assert rgenB.rgen is rgenD.rgen  # Test that the new swarm works


@pytest.mark.parametrize(
    ('settings', 'borg_class', 'expected'),
    (
        ({'rgen': 'rgen-borg', 'seed': 0}, RandomGeneratorBorg,
         'rgen-borg'),
        ({'rgen': 'mock-borg', 'seed': 0}, MockRandomGeneratorBorg,
         'mock-borg'),
    )
)
def test_borg_roundtrip_state(caplog, settings, borg_class, expected):
    """Test that Borg generators keep their type through a state roundtrip."""
    borg_class.reset_state()
    rgen = create_random_generator(settings)
    state = rgen.get_state()
    numbers1 = rgen.rand(shape=5)

    assert state['rgen'] == expected

    borg_class.reset_state()
    restored = create_random_generator(state)
    numbers2 = restored.rand(shape=5)

    assert isinstance(restored, borg_class)
    assert restored.get_state()['rgen'] == expected
    assert np.allclose(numbers1, numbers2)
    borg_class.reset_state()
