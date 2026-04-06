# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the engine classes."""
import logging
import pytest
import copy
from io import StringIO
from unittest.mock import patch
import numpy as np
from pyretis.engines import (
    MDEngine,
    Langevin,
    VelocityVerlet,
    Verlet,
)
from pyretis.core import System, create_box, Particles
from pyretis.core.common import big_fat_comparer
from pyretis.core.particlefunctions import (
    calculate_thermo,
    calculate_thermo_path,
)
from pyretis.core.path import Path
from pyretis.core.random_gen import MockRandomGenerator
from pyretis.orderparameter import PositionVelocity
from pyretis.forcefield import ForceField
from pyretis.forcefield.potentials.potentials import DoubleWell
logging.disable(logging.CRITICAL)


TISMOL_POS = np.array([-1.00000000, -0.99799432, -0.98665264, -0.97520352,
                       -0.96807892, -0.96159279, -0.95014721, -0.94447852,
                       -0.94017809, -0.93804918, -0.92699119, -0.92490484,
                       -0.91348463, -0.90195902, -0.89475989, -0.88820117,
                       -0.87668488, -0.87094729, -0.86657973, -0.86438545,
                       -0.85326381, -0.85111548, -0.83963493, -0.82805057,
                       -0.82079420, -0.81417970, -0.80260905, -0.79681846,
                       -0.79239923, -0.79015456, -0.77898382, -0.77678760,
                       -0.76526037, -0.75363047, -0.74632964, -0.73967169,
                       -0.72805861, -0.72222650, -0.71776665, -0.71548223,
                       -0.70427262, -0.70203832, -0.69047381, -0.67880735,
                       -0.67147062, -0.66477741, -0.65312967, -0.64726344,
                       -0.64276997, -0.64045243, -0.62921016])

TISMOL_VEL = np.array([0.78008020, 0.78725597, 0.79204310, 0.79490675,
                       0.79431301, 0.80064533, 0.80501899, 0.80970383,
                       0.81149089, 0.81460505, 0.82094536, 0.82705648,
                       0.83080101, 0.83265977, 0.83109156, 0.83647149,
                       0.83992168, 0.84371050, 0.84461736, 0.84686145,
                       0.85235247, 0.85763502, 0.86057204, 0.86165886,
                       0.85934749, 0.86400519, 0.86676039, 0.86988007,
                       0.87013281, 0.87173251, 0.87659871, 0.88127588,
                       0.88362726, 0.88416154, 0.88132438, 0.88547579,
                       0.88774997, 0.89041257, 0.89022221, 0.89138792,
                       0.89583821, 0.90011745, 0.90208908, 0.90227395,
                       0.89911189, 0.90295628, 0.90494650, 0.90734699,
                       0.90690729, 0.90783201, 0.91205777])

TISMOL_POS_HF = np.array([-1.00000000, -0.95540084, -0.19165497, 0.78018908,
                          1.32632619, 2.28604749, 2.39185991, 2.64430617,
                          2.25348211, 2.49825935, 2.88324951, 3.30550209,
                          2.98182643, 2.76033520, 2.86756839, 2.58347298,
                          2.50423709, 2.20295643, 2.42537267, 3.05261057,
                          3.15554594, 2.44639287, 2.88716400, 3.29914672,
                          2.96754315, 3.33636578, 2.79841879, 2.84221288,
                          2.33746733, 2.54907741, 2.90952943, 3.31484556,
                          2.98322875, 2.76077697, 2.86775262, 2.58354092,
                          2.50427057, 2.20297400, 2.42538389, 3.05261681,
                          3.15554769, 2.44639327, 2.88716422, 3.29914680,
                          2.96754316, 3.33636579, 2.79841879, 2.84221288,
                          2.33746733, 2.54907741, 2.90952943])

TISMOL_VEL_HF = np.array([0.78008020, 0.04459916, 0.76596773, 0.97676712,
                          0.53799599, 0.98657113, 0.36343554, 0.55356508,
                          0.03172585, 0.48984683, 0.73416686, 0.98453450,
                          0.55129904, 0.40598753, 0.59448391, 0.26823255,
                          0.31168371, 0.05072849, 0.44876367, 0.94301707,
                          0.78008020, 0.04459916, 0.76596773, 0.97676712,
                          0.53799599, 0.98657113, 0.36343554, 0.55356508,
                          0.03172585, 0.48984683, 0.73416686, 0.98453450,
                          0.55129904, 0.40598753, 0.59448391, 0.26823255,
                          0.31168371, 0.05072849, 0.44876367, 0.94301707,
                          0.78008020, 0.04459916, 0.76596773, 0.97676712,
                          0.53799599, 0.98657113, 0.36343554, 0.55356508,
                          0.03172585, 0.48984683, 0.73416686])

TISMOL_POS_VV = np.array([-1.00000000, -0.99843984, -0.99687973, -0.99531972,
                          -0.99375986, -0.99220019, -0.99064077, -0.98908165,
                          -0.98752287, -0.98596448, -0.98440654, -0.98284908,
                          -0.98129215, -0.97973581, -0.97818009, -0.97662505,
                          -0.97507074, -0.97351719, -0.97196445, -0.97041258,
                          -0.96886161])

TISMOL_VEL_VV = np.array([0.78008020, 0.78006775, 0.78003045, 0.77996842,
                          0.77988179, 0.77977066, 0.77963517, 0.77947542,
                          0.77929154, 0.77908365, 0.77885188, 0.77859634,
                          0.77831715, 0.77801444, 0.77768833, 0.77733895,
                          0.77696642, 0.77657086, 0.77615240, 0.77571116,
                          0.77524727])


def prepare_test_system():
    """Prepare a small system we can integrate."""
    box = create_box(periodic=[False])
    system = System(temperature=1.0, units='lj', box=box)
    system.particles = Particles(system.get_dim())
    pos = np.array([-1.0])
    vel = np.array([0.78008018])
    system.add_particle(name='Ar', pos=pos, vel=vel, mass=1.0, ptype=0)
    forcefield = ForceField('Double well force field for testing')
    param = {'a': 1.0, 'b': 2.0, 'c': 0.0}
    pot = DoubleWell()
    forcefield.add_potential(pot, parameters=param)
    system.forcefield = forcefield
    return system


class TestEngine:
    """Run the tests for the different Engines."""

    def test_restart_internal_engine_langevin(self):
        """Run test to check generic restart info."""
        system = [prepare_test_system(),
                  prepare_test_system(),
                  prepare_test_system()]

        eng = [
            Langevin(0.002, 0.3, rgen='mock', seed=1, high_friction=False),
            Langevin(0.002, 0.3, rgen='mock', seed=16, high_friction=False),
            Langevin(0.002, 0.3, rgen='mock', seed=1, high_friction=False)
            ]

        for _ in range(5):
            eng[0].integration_step(system[0])
            eng[1].integration_step(system[1])
            assert not np.allclose(system[0].particles.pos,
                                   system[1].particles.pos)

        for _ in range(16):
            eng[1].integration_step(system[1])
            eng[2].integration_step(system[2])
            assert not np.allclose(system[1].particles.pos,
                                   system[2].particles.pos)

        sav_en = copy.deepcopy(eng[0].restart_info())
        eng[1].load_restart_info(sav_en)

        sav_sys = copy.deepcopy(system[0].restart_info())
        system[1].load_restart_info(sav_sys)

        for _ in range(11):
            eng[1].integration_step(system[1])

        assert np.allclose(system[1].particles.pos,
                           system[2].particles.pos)
        assert np.allclose(system[1].particles.vel,
                           system[2].particles.vel)
        assert not np.allclose(system[0].particles.pos,
                               system[2].particles.pos)
        assert not np.allclose(system[0].particles.vel,
                               system[2].particles.vel)

    def test_restart_internal_engine_velocity_verlet(self):
        """Run test to check generic restart info."""

        system = [prepare_test_system(),
                  prepare_test_system(),
                  prepare_test_system()]

        eng = [VelocityVerlet(0.02),
               VelocityVerlet(0.02),
               VelocityVerlet(0.02)]

        for _ in range(5):
            eng[0].integration_step(system[0])
            eng[1].integration_step(system[1])

        for _ in range(16):
            eng[1].integration_step(system[1])
            eng[2].integration_step(system[2])
            assert not np.allclose(system[1].particles.pos,
                                   system[2].particles.pos)

        sav_sys = copy.deepcopy(system[0].restart_info())
        system[1].load_restart_info(sav_sys)

        for _ in range(11):
            eng[1].integration_step(system[1])

        assert np.allclose(system[1].particles.pos,
                           system[2].particles.pos)
        assert np.allclose(system[1].particles.vel,
                           system[2].particles.vel)
        assert not np.allclose(system[0].particles.pos,
                               system[2].particles.pos)
        assert not np.allclose(system[0].particles.vel,
                               system[2].particles.vel)

    def test_retarts_internal(self):
        """Check all internal engine restarts."""
        engs = [Verlet(timestep=0.02),
                VelocityVerlet(timestep=0.02),
                Langevin(timestep=0.02, gamma=0.1)]

        engs_base = [Verlet(timestep=1),
                     VelocityVerlet(timestep=1),
                     Langevin(timestep=1, gamma=1)]

        for engine, reference in zip(engs, engs_base):
            restart = engine.restart_info()
            reference.load_restart_info(copy.deepcopy(restart))
            if engine != reference:
                ergen = reference.rgen.get_state()
                rrgen = engine.rgen.get_state()
                # todo __eq__ method for rgen should be improved?
                assert big_fat_comparer(ergen, rrgen, hard=True)
                reference.rgen = engine.rgen
            assert engine == reference

    def test_restart_internal_engine_verlet(self):
        """Run test to check verlet error in restart."""
        system = [prepare_test_system(),
                  prepare_test_system(),
                  prepare_test_system()]

        eng = [Verlet(timestep=0.02),
               Verlet(timestep=0.02),
               Verlet(timestep=0.02)]
        eng[0].set_initial_positions(system[0].particles)
        eng[1].set_initial_positions(system[0].particles)
        eng[2].set_initial_positions(system[0].particles)

        for _ in range(6):
            eng[0].integration_step(system[0])
            eng[1].integration_step(system[1])

        for _ in range(16):
            eng[1].integration_step(system[1])
            eng[2].integration_step(system[2])
            assert not np.allclose(system[1].particles.pos,
                                   system[2].particles.pos)

        sav_sys = copy.deepcopy(system[0].restart_info())
        system[1].load_restart_info(sav_sys)
        eng[1].set_initial_positions(system[1].particles)

        for _ in range(10):
            eng[1].integration_step(system[1])

        assert not np.allclose(system[1].particles.pos,
                               system[2].particles.pos)
        assert not np.allclose(system[1].particles.vel,
                               system[2].particles.vel)
        assert not np.allclose(system[0].particles.pos,
                               system[2].particles.pos)
        assert not np.allclose(system[0].particles.vel,
                               system[2].particles.vel)

    def test_invert_dt(self):
        """Test that we can invert time."""
        eng = MDEngine(1, description='Mock engine for testing')
        invert = eng.invert_dt()
        assert not invert
        assert eng.timestep == pytest.approx(-1)
        invert = eng.invert_dt()
        assert invert
        assert eng.timestep == pytest.approx(1)

    def test_thermo(self):
        """Test that we can select the thermo function."""
        eng = MDEngine(1, description='Mock engine for testing')
        thermo = eng.select_thermo_function(thermo='full')
        assert thermo.__code__.co_code == calculate_thermo.__code__.co_code
        thermo = eng.select_thermo_function(thermo='path')
        assert (thermo.__code__.co_code ==
                calculate_thermo_path.__code__.co_code)
        thermo = eng.select_thermo_function(thermo='gibberish')
        assert (thermo.__code__.co_code ==
                calculate_thermo_path.__code__.co_code)
        thermo = eng.select_thermo_function(thermo=None)
        assert thermo is None

    def test_calculate_order(self):
        """Test that we can calculate the order parameter."""
        system = prepare_test_system()
        _box = create_box(periodic=[False], low=[0.0])
        system.box = _box
        eng = MDEngine(1, description='Mock engine for testing')
        order_function = PositionVelocity(0, dim='x', periodic=False)
        ensemble = {'system': system, 'order_function': order_function}
        order = eng.calculate_order(ensemble)
        assert order[0] == pytest.approx(system.particles.pos[0])
        assert order[1] == pytest.approx(system.particles.vel[0])
        pos = np.array([[1.0], ])
        vel = np.array([[-1.0], ])
        box = [100.0]
        order = eng.calculate_order(ensemble, xyz=pos, vel=vel, box=box)
        assert np.allclose(order, [1.0, -1.0])
        assert np.allclose(system.box.cell, box)

    def test_modify_velocities(self):
        """Test that the engine can modify velocities."""
        system = prepare_test_system()
        eng = MDEngine(1, description='Mock engine for testing')
        rgen = MockRandomGenerator(seed=1)
        ensemble = {'system': system, 'rgen': rgen}
        eng.modify_velocities(ensemble, {'aimless': True,
                                         'momentum': False})
        assert system.particles.vel[0][0] == pytest.approx(rgen.rgen[1])
        deltav = np.ones_like(system.particles.vel)
        eng.modify_velocities(ensemble, {'sigma_v': deltav,
                                         'aimless': False,
                                         'momentum': False})
        assert system.particles.vel[0][0] == rgen.rgen[1] + rgen.rgen[2]
        eng.modify_velocities(ensemble, {'aimless': True,
                                         'momentum': True})
        assert system.particles.vel[0][0] == pytest.approx(0.0)

        system.particles.pos = np.ones_like(system.particles.pos)
        system.particles.vpot = system.potential()
        dek, kin_new = eng.modify_velocities(ensemble,
                                             {'aimless': True,
                                              'momentum': False,
                                              'rescale': 6})
        assert 6 == pytest.approx(kin_new + system.particles.vpot)
        kin_old = kin_new
        dek, kin_new2 = eng.modify_velocities(ensemble,
                                              {'aimless': True,
                                               'momentum': False,
                                               'rescale': -1})
        assert kin_old == pytest.approx(kin_new2)
        assert dek == pytest.approx(0.0)

    def test_langevin_inertia(self):
        """Test the Langevin engine.

        In this test we compare the trajectory obtained by PyRETIS
        with a trajectory obtained independently by the TISMOL program.
        Since we are using a Langevin engine, we are using a mock
        random generator. This random generator is implemented identically
        in TISMOL and PyRETIS.
        """
        system = prepare_test_system()
        eng = Langevin(0.002, 0.3, rgen='mock', seed=1, high_friction=False)

        for i in range(51):
            assert np.allclose(system.particles.pos, TISMOL_POS[i])
            assert np.allclose(system.particles.vel, TISMOL_VEL[i])
            eng.integration_step(system)

    def test_langevin_negative_gamma(self):
        """Test that we fail for a negative gamma value."""
        system = prepare_test_system()
        eng = Langevin(0.002, -0.3, rgen='mock', seed=1, high_friction=False)
        with pytest.raises(ValueError):
            eng._init_parameters(system)  # pylint: disable=protected-access

    def test_langevin_hf(self):
        """Test the high friction variant of the Langevin engine."""
        system = prepare_test_system()
        eng = Langevin(0.002, 0.3, rgen='mock', seed=1, high_friction=True)
        for i in range(51):
            assert np.allclose(system.particles.pos, TISMOL_POS_HF[i])
            assert np.allclose(system.particles.vel, TISMOL_VEL_HF[i])
            eng.integration_step(system)

    def test_velocityverlet(self):
        """Test the Velocity Verlet engine."""
        system = prepare_test_system()
        eng = VelocityVerlet(0.002)
        for i in range(21):
            assert np.allclose(system.particles.pos, TISMOL_POS_VV[i])
            assert np.allclose(system.particles.vel, TISMOL_VEL_VV[i])
            eng.integration_step(system)

    def test_velocityverlet2(self):
        """Test that the engine is callable."""
        system = prepare_test_system()
        eng = VelocityVerlet(0.002)
        eng(system)
        assert np.allclose(system.particles.pos, TISMOL_POS_VV[1])
        assert np.allclose(system.particles.vel, TISMOL_VEL_VV[1])

    def test_verlet(self):
        """Test the Verlet engine."""
        system = prepare_test_system()
        eng = Verlet(0.002)
        eng.set_initial_positions(system.particles)
        assert eng.previous_pos[0][0] == pytest.approx(-1.00156016)
        correct_pos = np.array([-0.99843984, -0.99687973, -0.99531972,
                                -0.99375986, -0.99220019])
        correct_vel = np.array([0.78008018, 0.78006773, 0.78003043,
                                0.77996841, 0.77988177])
        for i in range(5):
            eng.integration_step(system)
            assert system.particles.pos[0][0] == pytest.approx(correct_pos[i])
            assert system.particles.vel[0][0] == pytest.approx(correct_vel[i])

    def test_add_to_path(self):
        """Test that we can add to paths from the engine."""
        eng = MDEngine(1.0, 'Just testing')

        path = Path(None, maxlen=10)
        system = System()
        system.particles = Particles(system.get_dim())
        left = -1.0
        right = 8.0
        for i in range(8):
            snapshot = {'order': [1.0 * i], 'pos': np.ones(3) * i,
                        'vel': np.ones(3) * i, 'vpot': i, 'ekin': i}
            phasepoint = eng.snapshot_to_system(system, snapshot)
            status, success, stop, add = eng.add_to_path(path, phasepoint,
                                                         left, right)
            assert status == 'Running propagate...'
            assert not success
            assert not stop
            assert add
        snapshot = {'order': [10.], 'pos': np.ones(3),
                    'vel': np.ones(3), 'vpot': 0.0, 'ekin': 0.0}
        phasepoint = eng.snapshot_to_system(system, snapshot)
        status, success, stop, add = eng.add_to_path(path, phasepoint,
                                                     left, right)
        assert status == 'Crossed right interface!'
        assert success
        assert stop
        assert add
        snapshot = {'order': [9.], 'pos': np.ones(3),
                    'vel': np.ones(3), 'vpot': 0.0, 'ekin': 0.0}
        phasepoint = eng.snapshot_to_system(system, snapshot)
        status, success, stop, add = eng.add_to_path(path, phasepoint,
                                                     left, right)
        assert status == 'Max. path length exceeded!'
        assert not success
        assert stop
        assert add
        phasepoint = {'order': [9.], 'pos': np.ones(3),
                      'vel': np.ones(3), 'vpot': 0.0, 'ekin': 0.0}
        status, success, stop, add = eng.add_to_path(path, phasepoint,
                                                     left, right)
        assert status == 'Max. path length exceeded!'
        assert not success
        assert stop
        assert not add

    def test_eq_and_ne(self):
        """Test that engines can be compared."""
        eng = MDEngine(1.0, 'Just testing')
        eng1 = MDEngine(1.0, 'Just Kidding')
        eng2 = VelocityVerlet(0.002)

        with patch('sys.stdout', new=StringIO()):
            assert eng != eng1
            assert eng != eng2

        eng = Langevin(0.002, 0.3, rgen='mock', seed=1, high_friction=False)
        eng1 = Langevin(0.002, 0.3, rgen='rgen', seed=1, high_friction=False)
        eng2 = Langevin(0.002, 0.3)
        with patch('sys.stdout', new=StringIO()):
            assert eng != eng1
            assert eng != eng2

        eng = [VelocityVerlet(0.02),
               VelocityVerlet(0.02),
               VelocityVerlet(0.02),
               VelocityVerlet(0.01)]

        assert eng[0] == eng[1]
        assert eng[0] == eng[2]
        with patch('sys.stdout', new=StringIO()):
            assert eng[0] != eng[3]

        eng = [Langevin(0.002, 0.3, rgen='mock', seed=1, high_friction=False),
               Langevin(0.002, 0.3, rgen='mock', seed=16, high_friction=False),
               Langevin(0.002, 0.3, rgen='mock', seed=1, high_friction=False)]

        assert eng[0] == eng[2]
        with patch('sys.stdout', new=StringIO()):
            assert eng[0] != eng[1]
            assert eng[0] != 'gino'

        eng[0].just = 'kidding'
        assert eng[0] != eng[2]
