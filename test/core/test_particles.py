# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the classes and methods from pyretis.core.particles"""
import pytest
import numpy as np
from pyretis.core.particles import (
    Particles,
    ParticlesExt,
    get_particle_type,
    particles_from_restart,
)
from scipy.special import comb


def test_creation():
    """Test the creation a particle list."""
    for dim in range(1, 4):
        particles = Particles(dim=dim)
        particles.add_particle(np.zeros(dim), np.zeros(dim), np.zeros(dim))
        particles.add_particle(np.zeros(dim), np.zeros(dim), np.zeros(dim))
        assert particles.dim == dim
        assert particles.npart == 2
        assert particles.virial is None


def test_empty_list():
    """Test that we can empty the list."""
    particles = Particles(dim=3)
    for _ in range(10):
        particles.add_particle(np.zeros(3), np.zeros(3), np.zeros(3))
    assert particles.npart == 10
    particles.empty_list()
    assert particles.npart == 0
    particles.add_particle(np.ones(3), np.ones(3), np.ones(3))
    assert particles.npart == 1


def test_pos_setget():
    """Test that we can set and get positions."""
    particles = Particles(dim=3)
    particles.add_particle(np.ones(3), np.ones(3), np.ones(3))
    particles.add_particle(np.ones(3), np.ones(3), np.ones(3))
    pos = particles.get_pos()
    assert pos is not particles.pos
    pos[0][1] = 100
    for i in particles.pos:
        assert np.allclose(i, np.ones(3))
    assert not np.allclose(pos[0], np.ones(3))
    particles.set_pos(pos)
    assert pos is not particles.pos
    for i, j in zip(pos, particles.pos):
        assert np.allclose(i, j)


def test_vel_setget():
    """Test that we can set and get velocities."""
    particles = Particles(dim=3)
    particles.add_particle(np.ones(3), np.ones(3), np.ones(3))
    particles.add_particle(np.ones(3), np.ones(3), np.ones(3))
    vel = particles.get_vel()
    assert vel is not particles.vel
    vel[0][1] = 100
    for i in particles.vel:
        assert np.allclose(i, np.ones(3))
    assert not np.allclose(vel[0], np.ones(3))
    particles.set_vel(vel)
    assert vel is not particles.vel
    for i, j in zip(vel, particles.vel):
        assert np.allclose(i, j)


def test_force_setget():
    """Test that we can set and get forces."""
    particles = Particles(dim=3)
    particles.add_particle(np.ones(3), np.ones(3), np.ones(3))
    particles.add_particle(np.ones(3), np.ones(3), np.ones(3))
    force = particles.get_force()
    assert force is not particles.force
    force[0][1] = 100
    for i in particles.force:
        assert np.allclose(i, np.ones(3))
    assert not np.allclose(force[0], np.ones(3))
    particles.set_force(force)
    assert force is not particles.force
    for i, j in zip(force, particles.force):
        assert np.allclose(i, j)


def test_get_selection():
    """Test that we can get a selection of properties."""
    particles = Particles(dim=3)
    for i in range(10):
        particles.add_particle(np.ones(3) * i, np.ones(3) * i,
                               np.ones(3) * i, name=f'H{i}')
    force = particles.get_selection(['force', 'name'], selection=None)[0]
    assert np.allclose(force, particles.force)
    force[0] = np.ones(3) * 101.
    assert np.allclose(force[0], np.array([101., 101., 101.]))
    out = particles.get_selection(['force', 'name'], selection=[1, 3])
    assert np.allclose(out[0][0], particles.force[1])
    assert np.allclose(out[0][1], particles.force[3])
    assert out[1][0] == 'H1'
    assert out[1][1] == 'H3'


def test_iterate():
    """Test that we can iterate over the particles."""
    particles = Particles(dim=3)
    for i in range(10):
        particles.add_particle(np.ones(3) * i, np.ones(3) * i,
                               np.ones(3) * i, name=f'A{i}')
    for i, part in enumerate(particles):
        assert np.allclose(part['pos'], particles.pos[i])


def test_pairs():
    """Test that we can iterate over pairs."""
    particles = Particles(dim=3)
    npart = 21
    for i in range(npart):
        particles.add_particle(np.ones(3) * i, np.ones(3) * i,
                               np.ones(3) * i, name=f'A{i}')
    pairs = set()
    for pair in particles.pairs():
        pairs.add((pair[0], pair[1]))
    assert len(pairs) == int(comb(npart, 2))


def test_generate_restart_info():
    """Test that we can generate restart info."""
    particles = Particles(dim=3)
    particles.add_particle(np.ones(3), np.ones(3), np.ones(3))
    particles.add_particle(np.ones(3) * 10, np.ones(3), np.ones(3))
    particles.vpot = 101
    particles.ekin = 102
    restart = particles.restart_info()
    assert restart['type'] == 'internal'
    particles2 = Particles(dim=3)
    particles2.load_restart_info(restart)
    assert np.allclose(particles.pos, particles2.pos)


def test_generate_restart_missing():
    """Test that we can generate restart info if attribs. are missing."""
    particles = Particles(dim=3)
    particles.add_particle(np.ones(3), np.ones(3), np.ones(3))
    particles.add_particle(np.ones(3) * 10, np.ones(3), np.ones(3))
    particles.vpot = 101
    if hasattr(particles, 'ekin'):
        del particles.ekin
    restart = particles.restart_info()
    particles2 = Particles(dim=3)
    particles2.load_restart_info(restart)
    assert np.allclose(particles.pos, particles2.pos)
    assert hasattr(particles2, 'ekin')
    assert not hasattr(particles, 'ekin')
    # Test load and create:
    simulation_restart = {'particles': restart}
    particles3 = particles_from_restart(simulation_restart)
    assert particles != particles3
    particles.ekin = None
    assert particles == particles3
    particles3 = particles_from_restart({})
    assert particles3 is None


def test_copy():
    """Test that we can copy particles."""
    particles = Particles(dim=3)
    particles.add_particle(np.ones(3), np.ones(3), np.ones(3))
    particles.add_particle(np.ones(3), np.ones(3), np.ones(3))
    particles.ekin = 102
    particles2 = particles.copy()
    assert particles is not particles2
    assert particles == particles2
    del particles.ekin
    particles2 = particles.copy()
    assert particles is not particles2
    assert particles != particles2
    # Test that we can change one without changing the other:
    particles2 = particles.copy()
    for attr in ('pos', 'vel', 'force'):
        val1 = getattr(particles2, attr)
        setattr(particles2, attr, np.zeros_like(val1))
        val1 = getattr(particles2, attr)
        val2 = getattr(particles, attr)
        assert not np.allclose(val1, val2)


def test_get_particle_type():
    """Test that we get correct classes."""
    cls1 = get_particle_type('internal')
    assert cls1 is Particles
    cls2 = get_particle_type('external')
    assert cls2 is ParticlesExt
    with pytest.raises(ValueError):
        get_particle_type('missing someone')


def test_particles_ext_creation():
    """Test the creation a particle list."""
    particles = ParticlesExt(dim=3)
    assert hasattr(particles, 'pos')
    assert hasattr(particles, 'vel')


def test_particles_ext_empty_list():
    """Test that we can empty the list."""
    particles = ParticlesExt(dim=3)
    particles.set_pos(('test', 100))
    assert particles.config == ('test', 100)
    particles.empty_list()
    assert particles.config == (None, None)


def test_particles_ext_setget_pos():
    """Test that we can set/get positions."""
    particles = ParticlesExt(dim=3)
    pos = ('some file', 101)
    assert particles.get_pos() == (None, None)
    particles.set_pos(pos)
    for i, j in zip(particles.get_pos(), pos):
        assert i == j
    assert particles.get_pos() is not pos
    pos2 = particles.get_pos()
    assert pos2 is particles.config


def test_particles_ext_set_vel():
    """Test that we can set the time direction for positions."""
    particles = ParticlesExt(dim=3)
    particles.set_vel(True)
    assert particles.vel_rev
    particles.set_vel(False)
    assert not particles.vel_rev


def test_particles_ext_restart_info():
    """Test restart methods."""
    particles = ParticlesExt(dim=3)
    pos = ('skinny love', 2)
    particles.set_pos(pos)
    particles.set_vel(True)
    restart = particles.restart_info()
    assert restart['vel_rev']
    for i, j in zip(particles.get_pos(), pos):
        assert i == j
    pos2 = ('My my my', 4)
    restart['config'] = pos2
    particles2 = ParticlesExt(dim=3)
    particles2.load_restart_info(restart)
    for i, j in zip(particles2.get_pos(), pos2):
        assert i == j
    del restart['config']
    with pytest.raises(ValueError):
        particles2.load_restart_info(restart)


def test_particles_ext_copy():
    """Test that we can copy particles."""
    particles = ParticlesExt(dim=3)
    particles.set_pos(('filename.ext', 987))
    particles.ekin = 123
    particles.vpot = 456
    particles2 = particles.copy()
    assert particles is not particles2
    assert particles == particles2
