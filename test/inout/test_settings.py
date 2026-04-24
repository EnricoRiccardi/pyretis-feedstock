# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test parsing from a settings input file.

Here we test that we parse the input file correctly and also that
we fail in predictable ways.
"""
import logging
import numpy as np
import os
import tempfile
import pytest
from io import StringIO
from pyretis.core.common import big_fat_comparer
from pyretis.core.units import create_conversion_factors, CONVERT
from pyretis.engines import Verlet, VelocityVerlet, Langevin
from pyretis.forcefield import PotentialFunction
from pyretis.forcefield.potentials import (
    DoubleWell,
    DoubleWellWCA,
    PairLennardJonesCut,
    PairLennardJonesCutnp,
    RectangularWell
)
from pyretis.inout.settings import (
    copy_settings,
    look_for_input_files,
    parse_settings_file,
    settings_to_text,
    write_settings_file,
    add_default_settings,
    add_specific_default_settings,
    settings_from_restart,
    _apply_restart_overrides,
    _clean_settings,
    _parse_raw_section,
    _parse_all_raw_sections,
    _parse_sections,
    RESTART_OVERRIDE_KEYWORDS,
    SECTIONS
)
from pyretis.orderparameter import (
    Angle,
    Dihedral,
    Distance,
    Distancevel,
    DistanceVelocity,
    OrderParameter,
    Position,
    PositionVelocity,
    Velocity,
)
from pyretis.setup.common import (
    create_engine,
    create_orderparameter
)
from pyretis.setup.createforcefield import (
    create_potentials,
    create_force_field
)
from pyretis.setup.createsystem import create_initial_positions
from unittest.mock import patch
from .help import turn_on_logging

logging.disable(logging.CRITICAL)

LOCAL_DIR = os.path.abspath(os.path.dirname(__file__))


class _ScalarInternalOrderParameter(OrderParameter):
    """Internal order parameter returning a scalar for one frame."""

    def __init__(self):
        super().__init__(description='Internal scalar order parameter')

    def calculate(self, system):
        """Return a scalar value for one frame."""
        return 1.25


class _MultiFrameInternalOrderParameter(OrderParameter):
    """Internal order parameter incorrectly returning multiple frames."""

    def __init__(self):
        super().__init__(description='Bad internal multi-frame output')

    def calculate(self, system):
        """Return multiple frames instead of one frame."""
        return [[1.25], [2.5]]


def _read_raw_settings(filename):
    """Read raw settings from a local file.

    Parameters
    ----------
    filename : string
        The file we are going to read.

    Return
    ------
    out : string
        The data read from the file.
    """
    inputfile = os.path.join(LOCAL_DIR, filename)
    data = None
    with open(inputfile, 'r', encoding='utf-8') as indata:
        data = indata.read()
    return data


def _test_correct_parsing(test, data, correct):
    """Helper method to test that we correctly parse settings.

    Parameters
    ----------
    test : object like unittest.TestCase
        The test that should pass or fail.
    data : string
        Raw data to be parsed.
    correct : dict
        The correct data.

    Returns
    -------
    out : dict
        The parsed settings.
    """
    raw = _parse_sections(data.split('\n'))
    settings = _parse_all_raw_sections(raw)
    for key in settings:
        assert settings[key] == correct[key]
    return settings


class TestKeywordParsing:
    """Test the parsing of input settings."""

    def test_parse_file(self):
        """Test that we can parse an input file."""
        inputfile = os.path.join(LOCAL_DIR, 'settings.rst')
        settings = parse_settings_file(inputfile)
        correct = {
            'heading': {'text': ('Molecular dynamics example settings\n'
                                 '===================================\n'
                                 'this is a simple test system to check '
                                 'the input statements are correct')},
            'system': {'units': 'lj',
                       'dimensions': 3,  # Added by default.
                       'temperature': 2.0},
            'simulation': {'task': 'md-nve',
                           'steps': 100,
                           'rgen': 'rgen',
                           'priority_shooting': False},
            'engine': {'class': 'velocityverlet',
                       'timestep': 0.002,
                       'rgen': 'rgen',
                       'type': 'internal'},
            'particles': {'position': {'input_file': 'initial.gro'},
                          'velocity': {'generate': 'maxwell',
                                       'temperature': 2.0,
                                       'momentum': True,
                                       'seed': 0},
                          'mass': {'Ar': 1.0},
                          'type': 'internal'},
            'forcefield': {'description': 'Lennard Jones test'},
            'potential': [{'class': 'PairLennardJonesCutnp',
                           'shift': True}],
        }
        del settings['engine']['exe_path']
        del settings['simulation']['exe_path']

        for key in correct:
            assert key in settings
            assert correct[key] == settings[key]

    def test_parse_file_without_units_uses_engine_default(self):
        """Test that units are filled from the selected engine."""
        settings_txt = """
Simulation settings
-------------------
task = md
steps = 10

Engine settings
---------------
class = velocityverlet
timestep = 0.002
"""
        with tempfile.NamedTemporaryFile('w', suffix='.rst',
                                         encoding='utf-8',
                                         delete=False) as handle:
            handle.write(settings_txt.strip())
            filename = handle.name
        try:
            settings = parse_settings_file(filename)
        finally:
            os.unlink(filename)
        assert settings['system']['units'] == 'lj'

    def test_parse_file_without_units_uses_openmm_default(self):
        """Test that OpenMM supplies its own default units."""
        settings_txt = """
Simulation settings
-------------------
task = md
steps = 10

Engine settings
---------------
class = openmm
openmm_simulation = sim
subcycles = 1
"""
        with tempfile.NamedTemporaryFile('w', suffix='.rst',
                                         encoding='utf-8',
                                         delete=False) as handle:
            handle.write(settings_txt.strip())
            filename = handle.name
        try:
            settings = parse_settings_file(filename)
        finally:
            os.unlink(filename)
        assert settings['system']['units'] == 'openmm'

    def test_keyword_format(self):
        """Test different forms of some simple keywords."""
        test_data = [
            (["units = 'lj'"], {'units': 'lj'}),  # Normal input
            (["unITS = 'lj'"], {'units': 'lj'}),  # Test case-sensitive?
            (["units = 'lj' # comment"], {'units': 'lj'}),  # Add a comment.
            (["units = 'lj' # comment units='a'"], {'units': 'lj'}),
            (["unITS='lj'"], {'units': 'lj'}),  # Test spacing.
            (["unITS= 'lj'"], {'units': 'lj'}),  # Test spacing.
            (['units = "lj"'], {'units': 'lj'}),  # Test " vs '.
            (['units =          "lj"'], {'units': 'lj'}),  # Test spacing.
            (['units        = "lj"'], {'units': 'lj'}),  # Test spacing.
            (['units = lj'], {'units': 'lj'}),  # Test quotations.
            (['units = "lj'], {'units': '"lj'})  # Test quotations.
        ]

        for data in test_data:
            setting = _parse_raw_section(data[0], 'system')
            assert setting == data[1]

    def test_keyword_dict(self):
        """Test some cases when reading dicts."""
        teststr = []
        correct = []
        # First == some simple input:
        teststr.append(['Engine settings',
                        '---------------',
                        'class = velocityverlet', 'timestep = 0.002'])
        correct.append({'timestep': 0.002, 'class': 'velocityverlet'})
        # Test input with comment(s):
        teststr.append(['Engine settings',
                        '---------------',
                        'class = velocityverlet', 'timestep = 0.002 # fs'])
        correct.append({'timestep': 0.002, 'class': 'velocityverlet'})
        # Test with white space, etc.:
        teststr.append(['Engine settings',
                        '---------------', 'junk',
                        '    ', 'junk = 10',
                        'class =    velocityverlet', 'timestep=0.002'])
        correct.append({'timestep': 0.002, 'junk': 10,
                        'class': 'velocityverlet'})
        teststr.append(['Engine settings',
                        '---------------', 'junk',
                        '    ', 'junk = 10',
                        'class =    velocityverlet', 'timestep=0.002'])
        correct.append({'timestep': 0.002, 'class': 'velocityverlet',
                        'junk': 10})
        teststr.append(['Engine settings',
                        '---------------',
                        'class = langevin', 'timestep = 0.002',
                        'gamma = 0.3', 'high_friction = False',
                        'seed = 0'])
        correct.append({'timestep': 0.002, 'class': 'langevin',
                        'gamma': 0.3, 'high_friction': False,
                        'seed': 0})

        for tst, corr in zip(teststr, correct):
            raw = _parse_sections(tst)
            setting = _parse_raw_section(raw['engine'], 'engine')
            assert setting == corr

    def test_write_and_read(self):
        """Test that we can parse some data == write it and read it."""
        data = _read_raw_settings('settings-rw-test.rst')
        correct = {'engine': {'timestep': 0.002,
                              'class': 'velocityverlet'},
                   'system': {'dimensions': 2, 'temperature': 1.0},
                   'simulation': {'task': 'md-nve', 'steps': 100}}
        raw = _parse_sections(data.split('\n'))
        settings = _parse_all_raw_sections(raw)
        settings = _test_correct_parsing(self, data, correct)
        with tempfile.NamedTemporaryFile() as temp:
            txt = settings_to_text(settings)
            temp.write(txt.encode('utf-8'))
            temp.flush()
            settings_read = parse_settings_file(temp.name, add_default=False)
        assert settings_read == correct

    def test_ignore_section(self):
        """Test that we ignore unknown sections."""
        data = _read_raw_settings('settings-unknown-section.rst')
        correct = {'engine': {'timestep': 0.002, 'class': 'velocityverlet'}}
        settings = _test_correct_parsing(self, data, correct)
        assert 'junk' not in settings


class TestKeywordEngine:
    """Test the parsing of input settings for engines."""

    def test_load_external_engine(self):
        """Test that we can load external python modules for engines."""
        data = _read_raw_settings('external-engine.rst')
        correct = {'engine': {'class': 'FooEngine', 'module': 'fooengine.py',
                              'timestep': 0.5,
                              'extra': 100}}
        settings = _test_correct_parsing(self, data, correct)
        # Here we add the exe_path key to the settings to tell
        # PyRETIS where we are executing from. This is to locate the
        # script we want to run.
        settings['simulation'] = {'exe_path': LOCAL_DIR}
        fooengine = create_engine(settings)
        assert fooengine.timestep == correct['engine']['timestep']
        # pylint: disable=no-member
        assert fooengine.extra == correct['engine']['extra']

    def test_fail_external_engine(self):
        """Test that external loads fail in a predicable way."""
        # First test: an engine that forgot to define a required method.
        test_data, correct = [], []
        test_data.append('---------------\n'
                         'Engine settings\n'
                         '---------------\n'
                         'class = BarEngine\n'
                         'module = fooengine.py')
        correct.append({'engine': {'class': 'BarEngine',
                                   'module': 'fooengine.py'}})
        test_data.append('Engine settings\n'
                         '---------------\n'
                         'class = BazEngine\n'
                         'module = fooengine.py')
        correct.append({'engine': {'class': 'BazEngine',
                                   'module': 'fooengine.py'}})
        test_data.append('Engine\n'
                         '------\n'
                         'module =  dummy')
        correct.append({'engine': {'module': 'dummy'}})
        test_data.append('Engine\n'
                         '------\n'
                         'module = dummy\n'
                         'class = dummy')
        correct.append({'engine': {'module': 'dummy', 'class': 'dummy'}})

        for data, corr in zip(test_data, correct):
            settings = _test_correct_parsing(self, data, corr)
            settings['simulation'] = {'exe_path': LOCAL_DIR}
            with pytest.raises(ValueError):
                create_engine(settings)

    def test_internal_engine(self):
        """Test that we can load all internal engines."""
        klass, test_data, correct = [], [], []
        klass.append(Verlet)
        test_data.append('Engine\n'
                         '------\n'
                         'class = Verlet\n'
                         'timestep = 0.5')
        correct.append({'engine': {'class': 'Verlet', 'timestep': 0.5}})
        klass.append(VelocityVerlet)
        test_data.append('Engine\n'
                         '------\n'
                         'class = VelocityVerlet\n'
                         'timestep = 0.314\n'
                         'desc = Test VV engine')
        correct.append({'engine': {'class': 'VelocityVerlet',
                                   'timestep': 0.314,
                                   'desc': 'Test VV engine'}})
        klass.append(Langevin)
        test_data.append('Engine\n'
                         '------\n'
                         'class = Langevin\n'
                         'timestep = 0.1\n'
                         'gamma = 2.718281828\n'
                         'seed = 101\n'
                         'high_friction = True')
        correct.append({'engine': {'class': 'Langevin',
                                   'timestep': 0.1,
                                   'gamma': 2.718281828,
                                   'seed': 101,
                                   'high_friction': True}})
        klass.append(Langevin)
        test_data.append('Engine\n'
                         '------\n'
                         'class = Langevin\n'
                         'timestep = 0.25\n'
                         'gamma = 2.718281828\n'
                         'seed = 11\n'
                         'high_friction = False')
        correct.append({'engine': {'class': 'Langevin',
                                   'timestep': 0.25,
                                   'gamma': 2.718281828,
                                   'seed': 11,
                                   'high_friction': False}})
        for data, corr, cls in zip(test_data, correct, klass):
            settings = _test_correct_parsing(None, data, corr)
            engine = create_engine(settings)
            assert isinstance(engine, cls)
            assert engine.timestep == \
                pytest.approx(corr['engine']['timestep'])
            for key, val in corr['engine'].items():
                if hasattr(engine, key):
                    assert getattr(engine, key) == val


class TestKeywordOrderParameter:
    """Test creation of order parameters."""

    @staticmethod
    def _internal_orderparameter_settings():
        """Return a minimal valid internal order-parameter section."""
        data = (
            'Orderparameter\n'
            '--------------\n'
            'class = Position\n'
            'index = 0\n'
            'dim = x\n'
            'periodic = False'
        )
        return _test_correct_parsing(
            None,
            data,
            {
                'orderparameter': {
                    'class': 'Position',
                    'index': 0,
                    'dim': 'x',
                    'periodic': False,
                }
            },
        )

    def test_load_orderparameter(self):
        """Test loading of external order parameter."""
        data = [
            'Orderparameter',
            '--------------',
            'class = FooOrderParameter',
            'module = fooorderparameter.py',
            'name = Dummy'
        ]
        correct = {
            'orderparameter': {
                'class': 'FooOrderParameter',
                'module': 'fooorderparameter.py',
                'name': 'Dummy'
            }
        }
        settings = _test_correct_parsing(self, '\n'.join(data), correct)
        # Here we add the exe-path key to the settings to tell
        # PyRETIS where we are executing from. This is to locate the
        # script we want to run.
        settings['simulation'] = {'exe_path': LOCAL_DIR}
        orderp = create_orderparameter(settings)
        assert correct['orderparameter']['class'] == orderp.__class__.__name__
        assert orderp.calculate(None) == [0.0]

    def test_external_orderparameter_scalar_is_normalized(self):
        """External scalar output is wrapped into a single-frame list."""
        data = (
            'Orderparameter\n'
            '--------------\n'
            'class = ScalarOrderParameter\n'
            'module = fooorderparameter.py'
        )
        settings = _test_correct_parsing(
            self,
            data,
            {
                'orderparameter': {
                    'class': 'ScalarOrderParameter',
                    'module': 'fooorderparameter.py',
                }
            },
        )
        settings['simulation'] = {'exe_path': LOCAL_DIR}
        orderp = create_orderparameter(settings)
        assert orderp.calculate(None) == [1.25]

    def test_internal_orderparameter_scalar_is_normalized(self):
        """Internal scalar output is wrapped into a single-frame list."""
        settings = self._internal_orderparameter_settings()
        with patch(
            'pyretis.setup.common.order_factory',
            return_value=_ScalarInternalOrderParameter(),
        ):
            orderp = create_orderparameter(settings)
        assert orderp.calculate(None) == [1.25]

    def test_external_orderparameter_array_is_normalized(self):
        """External array output is converted into a plain list."""
        data = (
            'Orderparameter\n'
            '--------------\n'
            'class = ArrayOrderParameter\n'
            'module = fooorderparameter.py'
        )
        settings = _test_correct_parsing(
            self,
            data,
            {
                'orderparameter': {
                    'class': 'ArrayOrderParameter',
                    'module': 'fooorderparameter.py',
                }
            },
        )
        settings['simulation'] = {'exe_path': LOCAL_DIR}
        orderp = create_orderparameter(settings)
        assert orderp.calculate(None) == [1.25, 2.5]

    def test_external_orderparameter_rejects_multi_frame_output(self):
        """External order parameters must return one frame at a time."""
        data = (
            'Orderparameter\n'
            '--------------\n'
            'class = MultiFrameOrderParameter\n'
            'module = fooorderparameter.py'
        )
        settings = _test_correct_parsing(
            self,
            data,
            {
                'orderparameter': {
                    'class': 'MultiFrameOrderParameter',
                    'module': 'fooorderparameter.py',
                }
            },
        )
        settings['simulation'] = {'exe_path': LOCAL_DIR}
        orderp = create_orderparameter(settings)
        with pytest.raises(ValueError, match='single frame'):
            orderp.calculate(None)

    def test_internal_orderparameter_rejects_multi_frame_output(self):
        """Internal order parameters must return one frame at a time."""
        settings = self._internal_orderparameter_settings()
        with patch(
            'pyretis.setup.common.order_factory',
            return_value=_MultiFrameInternalOrderParameter(),
        ):
            orderp = create_orderparameter(settings)
        with pytest.raises(ValueError, match='single frame'):
            orderp.calculate(None)

    def test_fail_orderparameter(self):
        """Test that loading external order parameters fails."""
        test_data, correct = [], []
        test_data.append('Orderparameter\n'
                         '--------------\n'
                         'class = BarOrderParameter\n'
                         'module = fooorderparameter.py')
        correct.append(
            {'orderparameter': {'class': 'BarOrderParameter',
                                'module': 'fooorderparameter.py'}})

        for data, corr in zip(test_data, correct):
            settings = _test_correct_parsing(None, data, corr)
            settings['simulation'] = {'exe_path': LOCAL_DIR}
            with pytest.raises(ValueError):
                create_orderparameter(settings)

    def test_create_orderparameter(self):
        """Test that we can create internal order parameters."""
        test_data = _read_raw_settings('internal-order.rst')
        klass = [
            OrderParameter,
            Position,
            Velocity,
            PositionVelocity,
            Distance,
            Distancevel,
            DistanceVelocity,
            Angle,
            Dihedral,
        ]
        correct = {
            'orderparameter': {
                'class': 'OrderParameter',
                'name': 'test'
            },
            'collective-variable':
                [
                    {
                        'class': 'Position',
                        'name': 'Position',
                        'index': 0, 'dim': 'x',
                        'periodic': False
                    },
                    {
                        'class': 'Velocity',
                        'name': 'Velocity',
                        'index': 0, 'dim': 'x',
                    },
                    {
                        'class': 'PositionVelocity',
                        'name': 'PositionVelocity',
                        'index': 0, 'dim': 'x',
                        'periodic': True
                    },
                    {
                        'class': 'Distance',
                        'name': 'My distance',
                        'index': (100, 101), 'periodic': False
                    },
                    {
                        'class': 'Distancevel',
                        'name': 'My distance Velocity',
                        'index': (102, 103), 'periodic': True
                    },
                    {
                        'class': 'DistanceVelocity',
                        'name': 'My distance Position Velocity',
                        'index': (1001, 1002), 'periodic': True
                    },
                    {
                        'class': 'Angle',
                        'name': 'Angle',
                        'index': (1001, 1002, 1003), 'periodic': True
                    },
                    {
                        'class': 'Dihedral',
                        'name': 'Dihedral',
                        'index': (1001, 1002, 1003, 1004), 'periodic': True
                    },
                ],
        }
        raw = _parse_sections(test_data.split('\n'))
        settings = _parse_all_raw_sections(raw)
        # Compare for order parameter:
        for key, val in settings['orderparameter'].items():
            assert val == correct['orderparameter'][key]
        for setting, corr in zip(settings['collective-variable'],
                                 correct['collective-variable']):
            for key, val in setting.items():
                assert val == corr[key]
        order = create_orderparameter(settings)
        for i, j in zip(order.order_parameters, klass):
            assert isinstance(i, j)


class TestKeywordParticles:
    """Test initialisation of particles."""

    def test_lattice(self):
        """Test initialisation on a lattice."""
        data = _read_raw_settings('settings-lattice.rst')
        correct = {'particles': {'position': {'generate': 'fcc',
                                              'repeat': [6, 6, 6],
                                              'lcon': 1.0}},
                   'system': {'units': 'lj'}}
        settings = _test_correct_parsing(self, data, correct)
        create_conversion_factors(settings['system']['units'])
        particles, size, _ = create_initial_positions(settings)
        correct_size = {'low': np.array([0., 0., 0.]),
                        'high': np.array([6., 6., 6.])}
        for key, val in correct_size.items():
            assert np.allclose(size[key], val)
        assert particles.npart == 4 * 6 * 6 * 6
        assert particles.mass[0][0] == pytest.approx(1.0)
        assert particles.imass[0][0] == pytest.approx(1.0)
        for i in range(particles.npart):
            assert particles.name[i] == 'Ar'

    def test_lattice_type(self):
        """Test initialisation on a lattice with types."""
        data = _read_raw_settings('settings-lattice-type.rst')
        correct = {'particles': {
            'position': {'generate': 'fcc', 'repeat': [3, 3, 3], 'lcon': 1.0},
            'ptype': [0, 1]},
            'system': {'units': 'lj'}}
        settings = _test_correct_parsing(self, data, correct)
        particles, _, _ = create_initial_positions(settings)
        for i in range(particles.npart):
            assert particles.name[i] == 'Ar'
        for i in range(particles.npart):
            if i == 0:
                assert particles.ptype[i] == 0
            else:
                assert particles.ptype[i] == 1

    def test_lattice_dens(self):
        """Test initialisation on a lattice with density set."""
        data = _read_raw_settings('settings-lattice-dens.rst')
        correct = {
            'particles': {
                'position': {'generate': 'fcc', 'repeat': [3, 3, 3],
                             'density': 0.9}},
            'system': {'units': 'lj'}}
        settings = _test_correct_parsing(self, data, correct)
        particles, size, _ = create_initial_positions(settings)
        correct_size = []
        lcon = 3.0 * (4.0 / 0.9)**(1.0 / 3.0)
        for _ in settings['particles']['position']['repeat']:
            correct_size.append([0.0, lcon])
        correct_size = np.array(correct_size)
        assert np.allclose(size['low'], correct_size[:, 0])
        assert np.allclose(size['high'], correct_size[:, 1])
        for i in range(particles.npart):
            assert particles.name[i] == 'Ar'
            assert particles.ptype[i] == 0

    def test_lattice_dens_lcon(self):
        """Test initialisation on a lattice with density and lcon set."""
        data = _read_raw_settings('settings-lattice-dens-lcon.rst')
        correct = {'particles': {'position': {'generate': 'fcc',
                                              'repeat': [3, 3, 3],
                                              'density': 0.9,
                                              'lcon': 1000.}},
                   'system': {'units': 'lj'}}
        settings = _test_correct_parsing(self, data, correct)
        _, size, _ = create_initial_positions(settings)
        correct_size = []
        # `lcon` should be replaced by density:
        lcon = 3.0 * (4.0 / 0.9)**(1.0 / 3.0)
        for _ in settings['particles']['position']['repeat']:
            correct_size.append([0.0, lcon])
        correct_size = np.array(correct_size)
        assert np.allclose(size['low'], correct_size[:, 0])
        assert np.allclose(size['high'], correct_size[:, 1])

    def test_lattice_and_mass(self):
        """Test initialisation on a lattice and setting of masses/types."""
        data = _read_raw_settings('settings-lattice-mass.rst')
        correct = {'particles': {'position': {'generate': 'fcc',
                                              'repeat': [3, 3, 3],
                                              'lcon': 1.},
                                 'ptype': [0, 1],
                                 'name': ['Ar', 'Kr'],
                                 'mass': {'Ar': 1.0}},
                   'system': {'units': 'lj'}}
        settings = _test_correct_parsing(self, data, correct)
        particles, _, _ = create_initial_positions(settings)
        for i in range(particles.npart):
            if i == 0:
                assert particles.ptype[i] == 0
                assert particles.name[i] == 'Ar'
                assert particles.mass[i][0] == pytest.approx(1.0)
            else:
                assert particles.ptype[i] == 1
                assert particles.name[i] == 'Kr'
                assert particles.mass[i][0] == pytest.approx(2.09767698)

    def test_inconsistent_dimlattice(self):
        """Test initialisation on a lattice with inconsistent dimensions."""
        data = _read_raw_settings('settings-lattice-dims.rst')
        correct = {'particles': {'position': {'generate': 'sq',
                                              'repeat': [6, 6],
                                              'lcon': 1.},
                                 'name': ['Ar'],
                                 'mass': {'Ar': 1.0}},
                   'system': {'units': 'lj', 'dimensions': 3}}
        settings = _test_correct_parsing(self, data, correct)
        args = [settings]
        with pytest.raises(ValueError):
            create_initial_positions(*args)

    def test_file_xyz(self):
        """Test initialisation from a XYZ file."""
        data = _read_raw_settings('initial-xyz2.rst')
        correct = {'particles': {'position': {'input_file': 'config.xyz'}},
                   'system': {'units': 'lj'}}
        settings = _test_correct_parsing(self, data, correct)
        units = settings['system']['units']
        create_conversion_factors(units)
        # Add path to the file for this test:
        settings['simulation'] = {'exe_path': LOCAL_DIR}
        particles, size, vel_read = create_initial_positions(settings)
        assert not vel_read
        assert size is None
        pos = particles.pos * CONVERT['length'][units, 'A']
        correct_pos = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5],
                                [0.5, 0.5, 0.0], [0.5, 0.0, 0.5],
                                [0.0, 0.5, 0.5]])
        assert np.allclose(pos, correct_pos)
        assert all([i == j for i, j in zip(particles.ptype,
                                           [0, 1, 2, 2, 2])])
        assert all([i == j for i, j in zip(particles.name,
                                           ['Ba', 'Hf', 'O',
                                            'O', 'O'])])
        masses = []
        for i in particles.mass:
            masses.append(i[0] * CONVERT['mass'][units, 'g/mol'])
        assert np.allclose(masses, [137.327, 178.49, 15.9994,
                                    15.9994, 15.9994])

    def test_file_gro(self):
        """Test initialisation from a GRO file."""
        data = _read_raw_settings('initial-gro.rst')
        correct = {'particles': {'position': {'input_file': 'config.gro'}},
                   'system': {'units': 'gromacs'}}
        settings = _test_correct_parsing(self, data, correct)
        # Add path to the file for this test:
        create_conversion_factors(settings['system']['units'])
        settings['simulation'] = {'exe_path': LOCAL_DIR}
        particles, size, vel_read = create_initial_positions(settings)
        assert vel_read
        assert np.allclose(size['cell'], [2., 2., 2.])
        correct_pos = np.array([[0., 0., 0.], [0.05, 0.05, 0.05],
                                [0.05, 0.05, 0.], [0.05, 0., 0.05],
                                [0., 0.05, 0.05]])
        assert np.allclose(particles.pos, correct_pos)
        assert all([i == j for i, j in zip(particles.ptype, [0, 1, 2, 2, 2])])
        assert all([i == j for i, j in zip(particles.name,
                                           ['Ba', 'Hf', 'O', 'O', 'O'])])
        assert np.allclose(particles.mass.T, [137.327, 178.49, 15.9994,
                                              15.9994, 15.9994])
        correct_vel = np.array([[1.0, 1.0, 1.0], [-1.0, -1.0, -1.0],
                                [2.0, 0.0, -2.0], [-2.0, 1.0, 2.0],
                                [0.0, -1.0, 0.0]])
        assert np.allclose(particles.vel, correct_vel)

    def test_file_xyztab(self):
        """Test initialisation from a XYZ file with mass dict."""
        data = _read_raw_settings('initial-xyz.rst')
        correct = {'particles': {'position': {'input_file': 'configtag.xyz'},
                                 'ptype': [0, 0, 0, 1, 1],
                                 'mass': {'Ar': 1., 'Kr': 2.09767698,
                                          'Kr2': 2.09767698}},
                   'system': {'units': 'lj'}}
        settings = _test_correct_parsing(self, data, correct)
        units = settings['system']['units']
        create_conversion_factors(units)
        # Add path to the file for this test:
        settings['simulation'] = {'exe_path': LOCAL_DIR}
        particles, size, vel_read = create_initial_positions(settings)
        assert not vel_read
        assert size is None
        assert all([i == j for i, j in zip(particles.ptype,
                                           [0, 0, 0, 1, 1])])
        assert all([i == j for i, j in zip(particles.name,
                                           ['Ar', 'Ar', 'Ar', 'Kr', 'Kr2'])])
        masses = []
        for i in particles.mass:
            masses.append(i[0] * CONVERT['mass'][units, 'g/mol'])
        assert np.allclose(masses, [39.948, 39.948, 39.948,
                                    83.798, 83.798])


class TestKeywordForcefield:
    """Test initialisation of force fields."""

    def test_forcefield(self):
        """Test initialisation of a simple force field."""
        data = _read_raw_settings('forcefield-simple.rst')
        correct = {'forcefield': {'description': 'My first force field'},
                   'potential': [{'class': 'PairLennardJonesCutnp',
                                  'shift': True,
                                  'parameter': {0: {'sigma': 1.0,
                                                    'epsilon': 1.0,
                                                    'rcut': 2.5}}}]}
        settings = _test_correct_parsing(self, data, correct)
        forcefield = create_force_field(settings)
        assert isinstance(forcefield.potential[0], PairLennardJonesCutnp)
        assert forcefield.potential[0].shift == \
            correct['potential'][0]['shift']

    def test_potential_parse(self):
        """Test creation of potentials while parsing input."""
        data = _read_raw_settings('settings-potential1.rst')
        correct = {'potential': [{'class': 'PairLennardJonesCut',
                                  'shift': False, 'mixing': 'geometric',
                                  'parameter': {0: {'sigma': 1.0,
                                                    'epsilon': 1.0,
                                                    'rcut': 2.5},
                                                1: {'sigma': 2.0,
                                                    'epsilon': 2.0,
                                                    'rcut': 2.5}}}]}
        settings = _test_correct_parsing(self, data, correct)
        potentials, pot_par = create_potentials(settings)
        assert isinstance(potentials[0], PairLennardJonesCut)
        assert potentials[0].shift == correct['potential'][0]['shift']
        # Test that we can assign parameters:
        for pot, params in zip(potentials, pot_par):
            pot.set_parameters(params)
        potparam = potentials[0].params
        assert potparam[(0, 0)]['epsilon'] == pytest.approx(1.0)
        assert potparam[(0, 0)]['sigma'] == pytest.approx(1.0)
        assert potparam[(0, 0)]['rcut'] == pytest.approx(2.5)
        assert potparam[(1, 1)]['epsilon'] == pytest.approx(2.0)
        assert potparam[(1, 1)]['sigma'] == pytest.approx(2.0)
        assert potparam[(1, 1)]['rcut'] == pytest.approx(2.5)
        assert potparam[(0, 1)]['epsilon'] == pytest.approx(1.414213562373095)
        assert potparam[(0, 1)]['sigma'] == pytest.approx(1.4142135623730951)
        assert potparam[(0, 1)]['rcut'] == pytest.approx(2.5)

    def test_empty_create_potential(self):
        """Test the create potential method."""
        correct = {'potential': [{'class': 'Fake'}]}
        potentials, pot_par = create_potentials(correct)
        for i, j in zip(potentials, pot_par):
            assert i is None
            assert j is None

    def test_potential_inconsitentdim(self):
        """Test creation of potentials with inconsistent dims."""
        data = _read_raw_settings('settings-potential2.rst')
        correct = {'system': {'dimensions': 2},
                   'potential': [{'class': 'PairLennardJonesCut',
                                  'shift': True,
                                  'parameter': {0: {'sigma': 1.0,
                                                    'epsilon': 1.0,
                                                    'rcut': 2.5}}}]}
        settings = _test_correct_parsing(self, data, correct)
        args = [settings]
        with pytest.raises(ValueError):
            create_potentials(*args)

    def test_potential_create(self):
        """Test that we can create all potentials."""
        data = _read_raw_settings('allpotentials.rst')
        all_potentials = [PairLennardJonesCut,
                          PairLennardJonesCutnp,
                          DoubleWellWCA,
                          DoubleWell,
                          RectangularWell]
        raw = _parse_sections(data.split('\n'))
        settings = _parse_all_raw_sections(raw)
        potentials, _ = create_potentials(settings)
        for pot, pot_input in zip(potentials, all_potentials):
            assert isinstance(pot, pot_input)

    def test_ext_potential(self):
        """Test creation of potentials while parsing input from externals."""
        data = _read_raw_settings('settings-ext-potential.rst')
        correct = {'potential': [{'class': 'FooPotential',
                                  'module': 'foopotential.py',
                                  'parameter': {'a': 2.0}}]}
        settings = _test_correct_parsing(self, data, correct)
        assert settings == correct
        # Add path for testing:
        settings['simulation'] = {'exe_path': LOCAL_DIR}
        potentials, pot_param = create_potentials(settings)
        assert isinstance(potentials[0], PotentialFunction)
        assert potentials[0].params['a'] == pytest.approx(0.0)
        for pot, pot_param in zip(potentials, pot_param):
            pot.set_parameters(pot_param)
        assert potentials[0].params['a'] == pytest.approx(2.0)

    def test_ext_potentialfail(self):
        """Test failure of external potential creation."""
        data = _read_raw_settings('settings-ext-potential-fail.rst')
        correct = {'potential': [{'class': 'BarPotential',
                                  'module': 'foopotential.py',
                                  'parameter': {'a': 2.0}}]}
        settings = _test_correct_parsing(self, data, correct)
        settings['simulation'] = {'exe_path': LOCAL_DIR}
        args = [settings]
        with pytest.raises(ValueError):
            create_potentials(*args)

    def test_complicated_input(self):
        """Test that we can read 'complex' force field input."""
        data = _read_raw_settings('forcefield.rst')
        correct = {'forcefield': {'description': 'My force field mix'},
                   'potential': [{'class': 'PairLennardJonesCutnp',
                                  'shift': True,
                                  'parameter': {0: {'sigma': 1.0,
                                                    'epsilon': 1.0,
                                                    'rcut': 2.5}}},
                                 {'class': 'DoubleWellWCA',
                                  'parameter': {'types': [(0, 0)],
                                                'rzero': 1. * (2.**(1./6.)),
                                                'height': 6.0,
                                                'width': 0.25}},
                                 {'class': 'FooPotential',
                                  'module': 'foopotential.py',
                                  'parameter': {'a': 10.0}}]}
        settings = _test_correct_parsing(self, data, correct)
        assert settings == correct
        settings['simulation'] = {'exe_path': LOCAL_DIR}
        forcefield = create_force_field(settings)
        assert len(forcefield.potential) == 3
        assert isinstance(forcefield.potential[0], PairLennardJonesCutnp)
        assert isinstance(forcefield.potential[1], DoubleWellWCA)
        assert isinstance(forcefield.potential[2], PotentialFunction)

    def test_too_many(self):
        """Test what happens when we add too many potentials."""
        settings = {
            'forcefield': {'description': 'My force field mix'},
            'potential': [],
        }
        for i in range(66):
            settings['potential'].append(
                {'class': 'DoubleWellWCA',
                 'parameter': {'types': [(i, i)],
                               'rzero': 1. * (2.**(1./6.))}}
            )
        txt = settings_to_text(settings)
        raw = _parse_sections(txt.split('\n'))
        settings2 = _parse_all_raw_sections(raw)
        for i, pot in enumerate(settings2['potential']):
            assert pot['parameter']['types'][0] == (i, i)

        # Looking for troubles:
        settings['ensemble'] = []
        for i in range(66):
            settings['ensemble'].append(
                {'potential': [{'class': 'DoubleWellWCA',
                                'parameter': {'types': [(i, i)],
                                              'rzero': 16}}]})

        txt = settings_to_text(settings)
        raw = _parse_sections(txt.split('\n'))
        settings3 = _parse_all_raw_sections(raw)
        hidden = settings3['ensemble'][10]['potential'][0]
        assert hidden['parameter']['rzero'] == 16

    def test_copy_settings(self):
        """Test that the copy method works as intended."""
        inputfile = os.path.join(LOCAL_DIR, 'settings-retis.rst')
        settings = parse_settings_file(inputfile)
        settings2 = copy_settings(settings)
        for key, val in settings.items():
            assert key in settings2
            assert val == settings2[key]

    def test_settings_to_txt(self):
        """Test the settings to text in more detail."""
        inputfile = os.path.join(LOCAL_DIR, 'settings-retis.rst')
        settings = parse_settings_file(inputfile)
        txt = settings_to_text(settings)
        to_find = ['RETIS EXAMPLE FOR TESTING',
                   'TIS', 'RETIS']
        found = [False for _ in to_find]
        for lines in txt.split('\n'):
            for i, key in enumerate(to_find):
                if found[i]:
                    continue
                if lines.startswith(key):
                    found[i] = True
        assert all(found)

    def test_write_settings(self):
        """Test that we can write settings to a file."""
        inputfile = os.path.join(LOCAL_DIR, 'settings-retis.rst')
        settings = parse_settings_file(inputfile)
        with tempfile.TemporaryDirectory() as tempdir:
            out_file = os.path.join(tempdir, '_pyretis_settings_temp.rst')
            # Check that we can write a file:
            write_settings_file(settings, out_file)
            assert os.path.isfile(out_file)
            # Check that we can write a file and backup:
            write_settings_file(settings, out_file, backup=True)
            out_file2 = f'{out_file}_000'
            with open(out_file, 'r', encoding='utf-8') as fileh:
                raw_sections = _parse_sections(fileh)
            settings2 = _parse_all_raw_sections(raw_sections)
            assert os.path.isfile(out_file2)
            assert big_fat_comparer(settings, settings2, hard=True)

    def test_look_for_files(self):
        """Test that we look and find the files."""
        input_files = look_for_input_files(input_path=LOCAL_DIR,
                                           required_files={},
                                           extra_files=[])
        assert input_files == {}

        with turn_on_logging(), patch('sys.stdout', new=StringIO()):
            # Check that ferrari.41 doesn't exist.
            look_for_input_files(input_path='.',
                                 required_files={},
                                 extra_files=['ferrari.f41'])

        # the index key should actually have an index.ndx value, but
        # I just use the available fake.mdp file instead...
        input_files = look_for_input_files(input_path=LOCAL_DIR,
                                           required_files={'mdp': 'gigi.mdp'},
                                           extra_files={'index': 'fake.mdp'})
        assert 'fake.mdp' in input_files['index']

        input_files = look_for_input_files(input_path=LOCAL_DIR,
                                           required_files={'mdp': 'gigi.mdp'},
                                           extra_files=[])
        assert 'fake.mdp' in input_files['mdp']
        input_files = look_for_input_files(input_path=LOCAL_DIR,
                                           required_files={'mdp': 'gigi.mdp'},
                                           extra_files=['traj.trr'])
        assert 'fake.mdp' in input_files['mdp']

        input_files = look_for_input_files(input_path=LOCAL_DIR,
                                           required_files={'mdp': 'fake.mdp'})
        assert 'fake.mdp' in input_files['mdp']

        input_files = look_for_input_files(input_path=LOCAL_DIR,
                                           required_files={'mdp': 'gigi.mdp',
                                                           'g96': 'dumped.g96',
                                                           'trr': 'traj.trr'})
        assert 'fake.mdp' in input_files['mdp']
        assert 'dumped.g96' in input_files['g96']
        assert 'traj.trr' in input_files['trr']

        look_for_input_files(input_path=LOCAL_DIR,
                             required_files={'cf': 'config.gro'},
                             extra_files=[])
        with pytest.raises(ValueError) as err:
            look_for_input_files(
                input_path=LOCAL_DIR,
                required_files={'lol': 'config.gRo'},
                extra_files=[])
        assert ('Missing input file "config.gRo" and multiple files '
                + 'have extension ".gro"') == str(err.value)

        with pytest.raises(ValueError) as err:
            input_files = look_for_input_files(input_path=LOCAL_DIR,
                                               required_files={'': 'ma.gro'})
        assert ('Missing input file "ma.gro" and multiple '
                + 'files have extension ".gro"') == str(err.value)

        with pytest.raises(ValueError) as err:
            input_files = look_for_input_files(input_path='banana',
                                               required_files={'': 'suka.gr'})
        assert 'Input path folder banana not existing' == str(err.value)

        with pytest.raises(ValueError) as err:
            input_files = look_for_input_files(
                input_path=LOCAL_DIR,
                required_files={'yyy': 'config.xxx'},
                extra_files=[])
        assert 'Missing input file "config.xxx" ' in str(err.value)

    def test_clean_settings(self):
        """Test that we can clean settings."""
        inputfile = os.path.join(LOCAL_DIR, 'settings-retis.rst')
        settings = parse_settings_file(inputfile)
        settings['junk'] = {'some': 'junk'}
        settings['system']['junk'] = 'this is junk'
        settings['box'] = []
        settings_c = _clean_settings(settings)
        assert 'junk' not in settings_c
        assert 'junk' not in settings_c['system']
        assert 'box' not in settings_c


class TestAddDefault:
    """Test adding default values."""

    def test_add_default_settings(self):
        """Test that we can add default settings."""
        settings_base = SECTIONS
        add_default_settings(settings_base)
        assert big_fat_comparer(settings_base, SECTIONS)

    def test_add_specific_default_settings(self, caplog):
        """Test that we can add default settings."""
        settings = {}
        add_default_settings(settings)

        settings['simulation']['task'] = 'tis'
        settings['engine']['class'] = 'velocityverlet'
        add_specific_default_settings(settings)
        assert settings['system']['temperature'] == 1
        assert settings['system']['units'] == 'lj'
        assert not settings['simulation']['flux']
        assert not settings['simulation']['zero_ensemble']

        settings['simulation']['task'] = 'retis'
        del settings['simulation']['flux']
        del settings['simulation']['zero_ensemble']
        add_specific_default_settings(settings)
        assert settings['simulation']['flux']
        assert settings['simulation']['zero_ensemble']

        assert settings['particles']['type'] == 'internal'
        assert settings['engine']['type'] == 'internal'

        settings['engine']['class'] = 'cp2k'
        del settings['system']['units']
        settings['engine']['exe_path'] += '/test/inout'
        add_specific_default_settings(settings)
        assert settings['particles']['type'] == 'external'
        assert settings['engine']['type'] == 'external'
        assert settings['system']['units'] == 'cp2k'
        input_path = settings['engine'].get('input_path', '.')
        add_specific_default_settings(settings)
        assert settings['engine']['input_files']['conf'][-11:] == 'initial.xyz'
        assert (settings['engine']['input_files']['template'][-9:] ==
                '/cp2k.inp')
        assert settings['engine']['input_files']['template'] == \
            os.path.join(os.path.abspath('.'), 'test/inout/',
                         input_path, 'cp2k.inp')
        assert settings['system']['temperature'] == 500
        del settings['system']['temperature']
        add_specific_default_settings(settings)
        assert settings['system']['temperature'] == 500
        no_temp_inp = os.path.join(os.path.abspath('.'),
                                   'test/inout/',
                                   input_path,
                                   'no_temp.inp')
        settings['engine']['template'] = no_temp_inp
        with turn_on_logging(), patch('sys.stdout', new=StringIO()):
            with caplog.at_level(logging.INFO):
                del settings['system']['temperature']
                add_specific_default_settings(settings)
                assert settings['system']['temperature'] == 300
                assert 'Temperature not set in CP2K' in caplog.text
                settings['system']['temperature'] = 500
                add_specific_default_settings(settings)
                assert settings['system']['temperature'] == 300
                assert 'And temperature in input' in caplog.text

        settings = {}
        add_default_settings(settings)
        settings['simulation']['task'] = 'md'
        settings['simulation']['exe_path'] = os.path.join(
            os.path.abspath('.'), 'test/engines'
        )
        settings['engine']['class'] = 'lammps'
        settings['engine']['input_path'] = 'lammps_input'
        add_specific_default_settings(settings)
        assert settings['system']['units'] == 'real'

        settings = {}
        settings['simulation'] = {'covid': 'kill_us_all'}
        settings['initial-path'] = {'method': 'restart'}
        add_default_settings(settings)
        settings['engine']['class'] = 'gromacs2'
        settings['engine']['topology'] = 'TrainToRide'
        add_specific_default_settings(settings)
        assert settings['particles']['type'] == 'external'
        assert settings['engine']['type'] == 'external'
        assert settings['system']['units'] == 'gromacs'
        input_path = settings['engine'].get('input_path', '.')
        assert settings['engine']['input_files']['topology'] == 'TrainToRide'
        assert settings['engine']['input_files']['conf'][-8:] == 'conf.gro'
        assert settings['engine']['input_files']['input_o'][-11:] == \
            '/grompp.mdp'
        assert settings['engine']['input_files']['index'] == \
            os.path.join(os.path.abspath('.'), input_path, 'index.ndx')
        assert settings['simulation']['restart'] == 'pyretis.restart'


class TestKeywordEnsemble:
    """Test ensemble and setting keywords."""

    def test_fill_up_tis_and_retis_settings(self):
        """Test if ensembles are correctly generated."""
        inputfile = os.path.join(LOCAL_DIR, 'settings-retis.rst')
        settings = parse_settings_file(inputfile)

        assert (settings['ensemble'][2]['collective-variable'][3]['something']
                ['unexpected'] == [1, 2, 3])
        assert (settings['ensemble'][6]['particles']
                ['velocity']['fantasy'] == 'game')
        assert (settings['ensemble'][6]['particles']['velocity']
                ['generate'] == 'Priapo')
        assert settings['ensemble'][5]['retis'] == settings['retis']
        assert settings['ensemble'][2]['box'] == settings['box']
        assert (settings['ensemble'][0]['tis']['ensemble_number'] !=
                settings['ensemble'][1]['tis']['ensemble_number'])
        assert (settings['ensemble'][2]
                ['collective-variable'][5]['name'] == 'Bugno')
        assert (settings['ensemble'][2]['collective-variable'][1]
                ['position']['nonsense'] == 'pineapple_on_pizza')
        assert (settings['ensemble'][1]['particles']
                ['velocity']['generate'] == 'maxwell')
        assert len(settings['ensemble']) == 8
        assert settings['simulation']['zero_left'] == -99


class TestRestartOverride:
    """Test the restart keyword override mechanism."""

    def _make_settings(self, steps=100, shooting_move='sh'):
        """Build a minimal settings dict for testing."""
        return {
            'simulation': {
                'steps': steps,
                'task': 'retis',
                'exe_path': '/tmp',
            },
            'tis': {
                'shooting_move': shooting_move,
                'freq': 0.5,
                'high_accept': True,
            },
            'retis': {
                'swapfreq': 0.05,
            },
        }

    def test_restart_override_keywords_structure(self):
        """RESTART_OVERRIDE_KEYWORDS must be a non-empty dict of lists."""
        assert isinstance(RESTART_OVERRIDE_KEYWORDS, dict)
        assert len(RESTART_OVERRIDE_KEYWORDS) > 0
        for section, keys in RESTART_OVERRIDE_KEYWORDS.items():
            assert isinstance(section, str)
            assert isinstance(keys, list)
            assert len(keys) > 0

    def test_restart_override_keywords_contains_wf_ss_wt(self):
        """tis section must include wf/ss/wt-specific keywords."""
        tis_keys = RESTART_OVERRIDE_KEYWORDS.get('tis', [])
        for kw in ('high_accept', 'n_jumps', 'mirror_freq',
                   'target_freq', 'target_indices'):
            assert kw in tis_keys, f'Missing wf/ss/wt keyword: {kw}'

    def test_apply_no_warning_when_values_match(self, caplog):
        """No warning when restarted and new values are identical."""
        new_set = self._make_settings(steps=100, shooting_move='sh')
        settings = self._make_settings(steps=100, shooting_move='sh')
        with turn_on_logging():
            with caplog.at_level(logging.WARNING):
                _apply_restart_overrides(new_set, settings)
        assert 'Restart override' not in caplog.text

    def test_apply_warning_on_mismatch(self, caplog):
        """A WARNING is logged when an allowed keyword value changes."""
        new_set = self._make_settings(steps=100, shooting_move='sh')
        settings = self._make_settings(steps=200, shooting_move='wf')
        with turn_on_logging():
            with caplog.at_level(logging.WARNING):
                _apply_restart_overrides(new_set, settings)
        assert 'Restart override' in caplog.text

    def test_apply_overrides_value(self):
        """The new value replaces the restarted one for allowed keywords."""
        new_set = self._make_settings(steps=100)
        settings = self._make_settings(steps=999)
        _apply_restart_overrides(new_set, settings)
        assert new_set['simulation']['steps'] == 999

    def test_apply_override_shooting_move(self):
        """shooting_move override works for wf/ss/wt keywords."""
        new_set = self._make_settings(shooting_move='sh')
        settings = self._make_settings(shooting_move='wf')
        _apply_restart_overrides(new_set, settings)
        assert new_set['tis']['shooting_move'] == 'wf'

    def test_apply_skips_missing_section(self):
        """Sections absent from either dict are silently skipped."""
        new_set = {'simulation': {'steps': 10}}
        settings = {'tis': {'freq': 0.5}}
        _apply_restart_overrides(new_set, settings)
        assert new_set == {'simulation': {'steps': 10}}

    def test_apply_skips_key_absent_from_new_input(self):
        """Keys absent from the new input are not touched in the restart."""
        new_set = {'simulation': {'steps': 10, 'interfaces': [0.1, 0.5]}}
        settings = {'simulation': {'steps': 20}}
        _apply_restart_overrides(new_set, settings)
        assert new_set['simulation']['steps'] == 20
        # interfaces is not in the new input, so the restart value is kept
        assert new_set['simulation']['interfaces'] == [0.1, 0.5]

    def test_apply_ensemble_override_warning(self, caplog):
        """Ensemble-level overrides emit a warning when values differ."""
        new_set = {
            'ensemble': [
                {'tis': {'shooting_move': 'sh', 'freq': 0.5}},
            ]
        }
        settings = {
            'ensemble': [
                {'tis': {'shooting_move': 'wf', 'freq': 0.5}},
            ]
        }
        with turn_on_logging():
            with caplog.at_level(logging.WARNING):
                _apply_restart_overrides(new_set, settings)
        assert 'Restart override (ensemble 0)' in caplog.text
        assert new_set['ensemble'][0]['tis']['shooting_move'] == 'wf'

    def test_apply_ensemble_override_no_warning_match(self, caplog):
        """No warning when ensemble values already match."""
        new_set = {
            'ensemble': [
                {'tis': {'shooting_move': 'sh'}},
            ]
        }
        settings = {
            'ensemble': [
                {'tis': {'shooting_move': 'sh'}},
            ]
        }
        with turn_on_logging():
            with caplog.at_level(logging.WARNING):
                _apply_restart_overrides(new_set, settings)
        assert 'Restart override' not in caplog.text

    def test_settings_from_restart_applies_overrides(self, tmp_path):
        """settings_from_restart honours new input values for allowed keys."""
        import pickle
        # Build a restart dict that looks like what pyretis stores.
        restarted_settings = {
            'simulation': {
                'steps': 50,
                'task': 'retis',
                'exe_path': str(tmp_path),
                'startcycle': 0,
                'restart': str(tmp_path / 'pyretis.restart'),
                'priority_shooting': False,
            },
            'engine': {
                'class': 'velocityverlet',
                'input_path': '.',
                'exe_path': str(tmp_path),
            },
            'tis': {'shooting_move': 'sh'},
        }
        restart = {
            'cycle': 50,
            'settings': restarted_settings,
        }
        restart_file = tmp_path / 'pyretis.restart'
        with open(restart_file, 'wb') as fh:
            pickle.dump(restart, fh)

        # New input requests more steps and a different shooting_move.
        new_settings = {
            'simulation': {
                'steps': 200,
                'task': 'retis',
                'exe_path': str(tmp_path),
                'restart': str(restart_file),
                'priority_shooting': False,
            },
            'engine': {
                'class': 'velocityverlet',
                'input_path': '.',
                'exe_path': str(tmp_path),
            },
            'tis': {'shooting_move': 'wf'},
        }
        result, _ = settings_from_restart(new_settings)
        assert result['simulation']['steps'] == 200
        assert result['tis']['shooting_move'] == 'wf'
