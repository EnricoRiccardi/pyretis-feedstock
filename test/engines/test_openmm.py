# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the OpenMMEngine class."""
import unittest
import pytest
import os
import numpy as np
import pyretis
from pyretis.engines import OpenMMEngine
from pyretis.engines.openmm import HAS_OPENMM, make_pyretis_system
from pyretis.core.path import Path
from pyretis.orderparameter import Distance
from pyretis.core.random_gen import RandomGenerator
from pyretis.setup.createsystem import create_system
from .test_engines import (prepare_test_system,
                           MockRandomGenerator)


if HAS_OPENMM:
    from .test_helpers.test_helpers import (create_openmm_simulation,
                                            write_test_pdb, FakeOp,
                                            write_openmm_simulation)
else:
    def create_openmm_simulation(*args, **kwargs):
        """Dummy function for when OpenMM is missing."""
        raise RuntimeError("OpenMM is not installed")

    def write_test_pdb(*args, **kwargs):
        """Dummy function for when OpenMM is missing."""
        pass

    def write_openmm_simulation(*args, **kwargs):
        """Dummy function for when OpenMM is missing."""
        pass

    class FakeOp:
        """Dummy class for when OpenMM is missing."""
        def __init__(self, *args, **kwargs):
            pass


@pytest.mark.skipif(not HAS_OPENMM, reason="OpenMM is not installed")
class TestOpenMMEngine:
    """Test the OpenMMEngine."""
    pdb = 'test_openmm_pyretis'
    f_name = pdb+'.pdb'

    openmm_simulation = 'simulation'
    openmm_module = 'openmm_module.py'

    @classmethod
    def setup_class(cls):
        if os.path.isfile(cls.f_name):
            i = 1
            while os.path.isfile(cls.pdb+str(i)+'.pdb'):
                i += 1
            cls.f_name = cls.pdb + str(i) + '.pdb'
        write_test_pdb(cls.f_name)
        write_openmm_simulation(cls.openmm_module, cls.f_name)

    @classmethod
    def teardown_class(cls):
        files = [cls.f_name, cls.openmm_module]
        for f in files:
            if os.path.isfile(f):
                os.remove(f)

    def setup_method(self):
        self.sim = create_openmm_simulation(self.f_name)
        # Distance between the two water oxygens periodic= false, because
        # triclininc boxes are hard:
        self.order_parameter = Distance((0, 3), periodic=False)
        self.settings = {'engine': {
            'openmm_module': self.openmm_module,
            'openmm_simulation': self.openmm_simulation}
                         }

        self.system = make_pyretis_system(self.settings)

    def test_init(self):
        """Test initiation of the OpenMMEngine."""
        eng = OpenMMEngine(openmm_simulation=self.sim, subcycles=4)
        assert eng.simulation is self.sim
        assert eng.subcycles == 4

    def test_init_module(self):
        eng = OpenMMEngine(openmm_simulation=self.openmm_simulation,
                           openmm_module=self.openmm_module)
        # Test if it runs
        eng.integration_step(self.system)

    def test_init_string_but_no_module(self):
        with pytest.raises(RuntimeError):
            OpenMMEngine(openmm_simulation=self.openmm_simulation)

    def test_installation_error(self):
        """Test failure when the OpenMMEngine is missing."""
        with pytest.raises(RuntimeError):
            pyretis.engines.openmm.HAS_OPENMM = False
            _ = OpenMMEngine(openmm_simulation=self.sim)
        # Reset change:
        pyretis.engines.openmm.HAS_OPENMM = True

    def test_create_system(self):
        settings = self.settings.copy()
        settings['engine']['type'] = 'OpenMM'
        create_system(settings)

    def test_integration_step(self):
        """Test the integration step of the OpenMMEngine."""
        eng = OpenMMEngine(openmm_simulation=self.sim)
        eng.integration_step(self.system)
        # Test if resetting works:
        eng.integration_step(self.system)

    def test_integration_step_no_box(self):
        """Test the integration step of the OpenMMEngine with no box."""
        eng = OpenMMEngine(openmm_simulation=self.sim)
        self.system.box = None
        eng.integration_step(self.system)
        # Test if resetting works:
        eng.integration_step(self.system)
        assert self.system.box is None

    def test_calculate_order(self):
        """Test order parameter calculation with the OpenMMEngine."""
        eng = OpenMMEngine(openmm_simulation=self.sim)
        answer = eng.calculate_order({'order_function': self.order_parameter,
                                      'system': self.system})
        assert self.order_parameter.calculate(self.system) == answer

    def test_calculate_order_xyz(self):
        """Test order parameter calculation with the OpenMMEngine."""
        eng = OpenMMEngine(openmm_simulation=self.sim)
        xyz = np.array([[1, 0, 0], [0, 0, 0], [0, 0, 0],
                        [-1, 0, 0], [0, 0, 0], [0, 0, 0]])
        answer = eng.calculate_order({'order_function': self.order_parameter,
                                      'system': self.system},
                                     xyz=xyz,
                                     vel=self.system.particles.vel,
                                     box=self.system.box.length)
        assert answer == [2]

    def test_pass_functions(self):
        """Test "empty" OpenMMEngine methods."""
        eng = OpenMMEngine(openmm_simulation=self.sim)
        eng.clean_up()
        eng.dump_phasepoint(1)

    def test_kick_across_middle_from_smaller(self):
        """Test the intialisation method with the OpenMMEngine."""
        eng = OpenMMEngine(openmm_simulation=self.sim)
        tis_settings = {'zero_momentum': True,
                        'rescale_energy': False}
        out = eng.kick_across_middle({'system': self.system,
                                      'order_function': self.order_parameter,
                                      'rgen': RandomGenerator()},
                                     middle=3.6,
                                     tis_settings=tis_settings)
        xyz0 = out[0].particles.get_pos()
        xyz1 = out[1].particles.get_pos()
        answer0 = eng.calculate_order({'order_function': self.order_parameter,
                                       'system': self.system},
                                      xyz=xyz0,
                                      vel=self.system.particles.vel,
                                      box=self.system.box.length)[0]
        answer1 = eng.calculate_order({'order_function': self.order_parameter,
                                       'system': self.system},
                                      xyz=xyz1,
                                      vel=self.system.particles.vel,
                                      box=self.system.box.length)[0]
        assert answer0 < answer1
        assert answer0 < 3.6
        assert answer1 > 3.6

    def test_kick_across_middle_from_bigger(self):
        """Test the intialisation method with the OpenMMEngine."""
        eng = OpenMMEngine(openmm_simulation=self.sim)
        tis_settings = {'zero_momentum': True,
                        'rescale_energy': False}
        out = eng.kick_across_middle({'system': self.system,
                                      'order_function': self.order_parameter,
                                      'rgen': RandomGenerator()},
                                     middle=3.55,
                                     tis_settings=tis_settings)
        xyz0 = out[0].particles.get_pos()
        xyz1 = out[1].particles.get_pos()
        answer0 = eng.calculate_order({'order_function': self.order_parameter,
                                       'system': self.system},
                                      xyz=xyz0,
                                      vel=self.system.particles.vel,
                                      box=self.system.box.length)[0]
        answer1 = eng.calculate_order({'order_function': self.order_parameter,
                                       'system': self.system},
                                      xyz=xyz1,
                                      vel=self.system.particles.vel,
                                      box=self.system.box.length)[0]
        assert answer0 > answer1
        assert answer0 > 3.55
        assert answer1 < 3.55

    def test_kick_across_middle_if_on_interface(self):
        """Test the intialisation method with the OpenMMEngine."""
        eng = OpenMMEngine(openmm_simulation=self.sim)
        tis_settings = {'zero_momentum': True,
                        'rescale_energy': False}

        # Init frame is ON the interface:
        middle = eng.calculate_order({'order_function': self.order_parameter,
                                      'system': self.system},
                                     xyz=self.system.particles.pos,
                                     vel=self.system.particles.vel,
                                     box=self.system.box.length)[0]

        out = eng.kick_across_middle({'system': self.system,
                                      'order_function': self.order_parameter,
                                      'rgen': RandomGenerator()},
                                     middle=middle,
                                     tis_settings=tis_settings)
        xyz0 = out[0].particles.get_pos()
        xyz1 = out[1].particles.get_pos()
        answer0 = eng.calculate_order({'order_function': self.order_parameter,
                                       'system': self.system},
                                      xyz=xyz0,
                                      vel=self.system.particles.vel,
                                      box=self.system.box.length)[0]
        answer1 = eng.calculate_order({'order_function': self.order_parameter,
                                       'system': self.system},
                                      xyz=xyz1,
                                      vel=self.system.particles.vel,
                                      box=self.system.box.length)[0]
        # Order of fluctuation around OP is not set:
        if answer0 < answer1:
            assert answer0 > middle
            assert answer1 < middle
        else:
            assert answer1 < middle
            assert answer0 > middle

    def test_kick_across_middle_if_first_away(self):
        """Test the intialisation method with the OpenMMEngine."""
        eng = OpenMMEngine(openmm_simulation=self.sim)
        tis_settings = {'zero_momentum': True,
                        'rescale_energy': False}

        # 2 op calculations per frame force OP to run away from:
        order_parameter = FakeOp([1, 1, 2, 2, -1, -1])
        _ = eng.kick_across_middle({'system': self.system,
                                    'order_function': order_parameter,
                                    'rgen': RandomGenerator()},
                                   middle=0,
                                   tis_settings=tis_settings)
        # Should have done 5 calls:
        assert order_parameter.n == 5

    def test_modify_velocities(self):
        """Test that we can modify velocities with the OpenMMEngine.

        Copied from pyretis.test.engines.test_engine.py
        """
        system = prepare_test_system()
        eng = OpenMMEngine(1)
        rgen = MockRandomGenerator(seed=1)
        eng.modify_velocities({'system': system, 'rgen': rgen},
                              {'sigma_v': None, 'aimless': True,
                               'momentum': False, 'rescale': None})
        assert system.particles.vel[0][0] == pytest.approx(rgen.rgen[1])
        deltav = np.ones_like(system.particles.vel)
        eng.modify_velocities({'system': system, 'rgen': rgen},
                              {'sigma_v': deltav, 'aimless': False,
                               'momentum': False, 'rescale': None})
        assert (system.particles.vel[0][0] ==
                pytest.approx(rgen.rgen[1] + rgen.rgen[2]))
        eng.modify_velocities({'system': system, 'rgen': rgen},
                              {'sigma_v': None, 'aimless': True,
                               'momentum': True, 'rescale': None})
        assert system.particles.vel[0][0] == pytest.approx(0.0)

        system.particles.pos = np.ones_like(system.particles.pos)
        system.particles.vpot = system.potential()
        dek, kin_new = eng.modify_velocities({'system': system, 'rgen': rgen},
                                             {'sigma_v': None, 'aimless': True,
                                              'momentum': False, 'rescale': 6})
        assert 6 == pytest.approx(kin_new + system.particles.vpot)
        kin_old = kin_new
        dek, kin_new2 = eng.modify_velocities({'system': system, 'rgen': rgen},
                                              {'sigma_v': None,
                                               'aimless': True,
                                               'momentum': False,
                                               'rescale': -1})
        assert kin_old == pytest.approx(kin_new2)
        assert dek == pytest.approx(0.0)

    def propagate_test(self, reverse=False):
        """Test the propagation method with the OpenMMEngine."""
        path = Path(None, maxlen=10)
        ifaces = [1, 2, 3]
        order_parameter = FakeOp([1, 2, 2.5, 3, 4])
        eng = OpenMMEngine(openmm_simulation=self.sim)
        ensemble = {'system': self.system,
                    'order_function': order_parameter,
                    'interfaces': ifaces}

        eng.propagate(path=path, ensemble=ensemble, reverse=reverse)

        assert order_parameter.n == 5

    def test_propagate_forward(self):
        """Test the propagation method (forward) with the OpenMMEngine."""
        self.propagate_test(reverse=False)

    def test_propagate_backward(self):
        """Test the propagation method (backward) with the OpenMMEngine."""
        self.propagate_test(reverse=True)
