# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the methods in pyretis.inout.setup.createsystem."""
import os
import logging
import filecmp
import pytest
import tempfile
import numpy as np
from pyretis.core.system import System
from pyretis.inout.common import make_dirs
from pyretis.engines.external import ExternalMDEngine
from pyretis.engines import CP2KEngine
from pyretis.inout.formats.xyz import (
    read_xyz_file,
    convert_snapshot,
)
from pyretis.orderparameter.orderparameter import PositionVelocity
from pyretis.core.particles import ParticlesExt
from .test_helpers.test_helpers import remove_dir
logging.disable(logging.CRITICAL)


HERE = os.path.abspath(os.path.dirname(__file__))


class DummyExternal(ExternalMDEngine):
    """A dummy external engine. Only useful for testing."""

    def __init__(self, input_path, timestep, subcycles):
        """Initialise the dummy engine."""
        super().__init__('External engine for testing!', timestep,
                         subcycles)
        self.input_path = os.path.abspath(input_path)
        self.input_files = {
            'conf': os.path.join(self.input_path, 'dummy_config'),
            'template': os.path.join(self.input_path, 'dummy_template'),
        }

    def step(self, system, name):
        """Perform a single dummy step."""

    def _read_configuration(self, filename):
        """Read a xyz configuration."""
        xyz, vel, box, names = None, None, None, None
        for snapshot in read_xyz_file(filename):
            box, xyz, vel, names = convert_snapshot(snapshot)
            break
        return box, xyz, vel, names

    def _reverse_velocities(self, filename, outfile):
        """Reverse velocities with a dummy method."""

    def _extract_frame(self, traj_file, idx, out_file):
        """Extract a frame, dummy method."""

    def _propagate_from(self, name, path, ensemble, msg_file, reverse=False):
        """Propagate with the engine, dummy method."""

    def modify_velocities(self, ensemble, vel_settings):
        """Modify velocities, dummy method."""


class TestExternalEngine:
    """Run tests for the external engine."""

    def test_exe_dir(self, caplog):
        """Test exe_dir property."""
        engine = DummyExternal('test', 1.0, 10)
        assert engine.exe_dir is None
        logging.disable(logging.INFO)
        with caplog.at_level(logging.WARNING, logger='pyretis.engines.engine'):
            engine.exe_dir = 'non-existing-dir'
        logging.disable(logging.CRITICAL)
        assert 'non-existing-dir' == engine.exe_dir

    def test_integration_step(self, caplog):
        """Test that the integration step fails for the external integrator."""
        engine = DummyExternal('test', 1.0, 10)
        with pytest.raises(NotImplementedError):
            engine.integration_step(None)

    def test_read_input_settings(self, caplog):
        """Test that we can read input settings."""
        filename = os.path.join(HERE, 'input.txt')
        engine = DummyExternal('test', 1.0, 10)
        # pylint: disable=protected-access
        settings = engine._read_input_settings(filename)
        correct = {'b': '100', 'a': '1', 'another setting': 'text'}
        for key, val in correct.items():
            assert val == settings[key]

    def test_modify_input(self, caplog):
        """Test that we can modify input settings."""
        engine = DummyExternal('test', 1.0, 10)
        filename = os.path.join(HERE, 'input.txt')
        settings = {'another setting': 'hello', 'c': '101'}
        # pylint: disable=protected-access
        with tempfile.NamedTemporaryFile() as temp:
            engine._modify_input(filename, temp.name, settings)
            temp.flush()
            settings2 = engine._read_input_settings(temp.name)
            correct = {'b': '100', 'a': '1', 'another setting': 'hello',
                       'c': '101'}
            for key, val in correct.items():
                assert val == settings2[key]

    def test_move_file(self, caplog):
        """Test that we can move files."""
        engine = DummyExternal('test', 1.0, 10)
        filename = os.path.join(HERE, 'empty_file')
        with open(filename, 'w', encoding="utf8") as temp:
            temp.write('Hello!')
        outfile = os.path.join(HERE, 'empty_file2')
        # pylint: disable=protected-access
        engine._movefile(filename, outfile)
        assert os.path.isfile(outfile)
        assert not os.path.isfile(filename)
        engine._removefile(outfile)
        engine._removefile(filename)

    def test_copy_file(self, caplog):
        """Test that we can copy files."""
        engine = DummyExternal('test', 1.0, 10)
        filename = os.path.join(HERE, 'empty_file')
        with open(filename, 'w', encoding="utf8") as temp:
            temp.write('Hello!')
        outfile = os.path.join(HERE, 'empty_file_copy')
        # pylint: disable=protected-access
        engine._copyfile(filename, outfile)
        assert os.path.isfile(outfile)
        assert os.path.isfile(filename)
        compare = filecmp.cmp(filename, outfile)
        assert compare
        engine._removefile(outfile)
        engine._removefile(filename)

    def test_removefiles(self, caplog):
        """Test that we can remove several files."""
        engine = DummyExternal('test', 1.0, 10)
        files = []
        for i in range(3):
            basename = f'empty_file{i}'
            filename = os.path.join(HERE, basename)
            with open(filename, 'w', encoding="utf8") as temp:
                temp.write('Hello!')
            files.append(basename)
        # pylint: disable=protected-access
        engine._remove_files(HERE, files)
        for i in files:
            filename = os.path.join(HERE, i)
            assert not os.path.isfile(filename)

    def test_cleanup(self, caplog):
        """Test the cleanup method."""
        dirname = os.path.join(HERE, 'testdir')
        make_dirs(dirname)
        files = []
        for i in range(3):
            basename = f'empty_file{i}'
            filename = os.path.join(dirname, basename)
            with open(filename, 'w', encoding="utf8") as temp:
                temp.write('Hello!')
            files.append(filename)
        engine = DummyExternal('test', 1.0, 10)
        engine.exe_dir = dirname
        engine.clean_up()
        for i in files:
            assert not os.path.isfile(i)
        remove_dir(dirname)

    def test_read_configuration(self, caplog):
        """Test that we can read a configuration."""
        engine = DummyExternal('test', 1.0, 10)
        filename = os.path.join(HERE, 'config.xyz')
        # pylint: disable=protected-access
        box, xyz, vel, atoms = engine._read_configuration(filename)
        correct_box = np.array([1., 2., 3.])
        assert np.allclose(box, correct_box)
        correct_atoms = ['Ba', 'Hf', 'O', 'O', 'O']
        for i, j in zip(correct_atoms, atoms):
            assert i == j
        correct_xyz = np.array([[0., 0., 0.], [0.5, 0.5, 0.5],
                                [0.5, 0.5, 0.], [0.5, 0., 0.5],
                                [0., 0.5, 0.5]])
        assert np.allclose(xyz, correct_xyz)
        correct_vel = np.array([[1., 1., 1.], [2., 2., 2.], [3., 3., 3.],
                                [4., 4., 4.], [5., 5., 5.]])
        assert np.allclose(vel, correct_vel)

    def test_calculate_order(self, caplog):
        """Test calculation of the order parameter."""
        engine = DummyExternal('test', 1.0, 10)
        filename = os.path.join(HERE, 'config.xyz')
        order_function = PositionVelocity(0, dim='x', periodic=False)
        system = System()
        system.particles = ParticlesExt(dim=3)
        system.particles.config = (filename, 0)
        ensemble = {'order_function': order_function, 'system': system}
        order = engine.calculate_order(ensemble)
        assert order[0] == pytest.approx(0.0)
        assert order[1] == pytest.approx(1.0)
        system.particles.vel_rev = True
        order = engine.calculate_order(ensemble)
        assert order[1] == pytest.approx(-1.0)
        engine.ext = 'A_cool_estention'
        engine.exe_dir = '.'
        engine.dump_phasepoint(system)
        assert engine.ext in system.particles.config[0]

    def test_calculate_order_external(self, caplog):
        """Test calculation of the order parameter with external engine."""
        cmd = os.path.join(HERE, 'mockcp2k.py')
        dir_name = os.path.join(HERE, 'cp2k_input')
        extra_files = ['extra_file']
        engine = CP2KEngine(cmd,
                            input_path=dir_name,
                            timestep=0.002,
                            subcycles=10,
                            extra_files=extra_files)
        order_function = PositionVelocity(0, dim='x', periodic=False)
        system = System()
        system.particles = ParticlesExt(dim=3)
        system.particles.config = (engine.input_files['conf'], 0)
        ensemble = {'order_function': order_function, 'system': system,
                    'engine': engine}
        order = engine.calculate_order(ensemble)
        # Check that have assigned the box to correctly calculate the order.
        assert order[0] == 0.7731715808
        assert order[1] == 1.

    def test_execute_command(self, caplog):
        """Test what happens when the execution fails."""
        engine = DummyExternal('test', 1.0, 10)
        cmd = [os.path.join(HERE, 'aprogram.py')]
        engine.execute_command(cmd, cwd=HERE, inputs=b'')
        cmd.append('arg')
        with pytest.raises(RuntimeError):
            engine.execute_command(cmd, cwd=HERE, inputs=b'')
        # The outputs should be retained after the previous error:
        with open(os.path.join(HERE, 'engine.log'), 'r',
                  encoding="utf8") as stdout:
            content = stdout.read()
            assert 'This is a program for testing external' in content
        with open(os.path.join(HERE, 'engine.err'), 'r',
                  encoding="utf8") as stderr:
            content = stderr.read()
            assert 'ERROR: Program got arguments:' in content
            assert ' '.join(cmd) in content
        for fname in ('engine.log', 'engine.err'):
            # pylint: disable=protected-access
            engine._removefile(os.path.join(HERE, fname))
