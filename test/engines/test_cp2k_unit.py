# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Unit tests for the CP2KEngine class."""
import os
import tempfile
from unittest.mock import MagicMock, patch
import pytest
import numpy as np
from pyretis.engines.cp2k import (
    CP2KEngine, write_for_step_vel, write_for_integrate,
    write_for_continue
)


def test_write_functions():
    """Test the helper functions for writing CP2K input files."""
    with patch('pyretis.engines.cp2k.update_cp2k_input') as mock_update:
        # Test write_for_step_vel
        vel = np.array([[0.1, 0.2, 0.3]])
        write_for_step_vel('in', 'out', 0.5, 10, 'pos', vel)
        mock_update.assert_called()

        # Test write_for_integrate
        write_for_integrate('in', 'out', 0.5, 10, 'pos')
        mock_update.assert_called()

        # Test write_for_continue
        write_for_continue('in', 'out', 0.5, 10)
        mock_update.assert_called()


@patch('pyretis.engines.cp2k.look_for_input_files')
def test_cp2k_engine_init(mock_look):
    """Test CP2KEngine initialization."""
    mock_look.return_value = {'conf': 'init.xyz', 'template': 'cp2k.inp'}
    with tempfile.TemporaryDirectory() as tempdir:
        engine = CP2KEngine(cp2k='cp2k_exe', input_path=tempdir,
                            timestep=0.5, subcycles=10)
        assert engine.cp2k == ['cp2k_exe']
        assert engine.timestep == 0.5
        assert engine.subcycles == 10


class TestCP2KEngineMethods:
    """Test methods of the CP2KEngine class."""

    @pytest.fixture
    def engine_env(self):
        """Fixture to set up a full engine env with temp directories."""
        with tempfile.TemporaryDirectory() as tempdir:
            input_path = os.path.join(tempdir, 'input')
            os.makedirs(input_path)
            exe_path = os.path.join(tempdir, 'exe')
            os.makedirs(exe_path)

            # Create some dummy input files
            with open(os.path.join(input_path, 'cp2k.inp'), 'w') as f:
                f.write('dummy template')
            with open(os.path.join(input_path, 'initial.xyz'), 'w') as f:
                f.write('2\n# Step: 0 Box: 1 1 1\n'
                        'H 0 0 0 0 0 0\nH 1 1 1 0 0 0')

            engine = CP2KEngine(cp2k='cp2k', input_path=input_path,
                                timestep=0.1, subcycles=1)
            engine.exe_dir = exe_path
            yield engine, tempdir

    def test_run_cp2k(self, engine_env):
        """Test run_cp2k method."""
        engine, _ = engine_env
        with patch.object(engine, 'execute_command') as mock_exec:
            out = engine.run_cp2k('input.inp', 'proj')
            mock_exec.assert_called()
            assert out['energy'].endswith('proj-1.ener')

    @patch('pyretis.engines.cp2k.read_xyz_file')
    @patch('pyretis.engines.cp2k.convert_snapshot')
    @patch('pyretis.engines.cp2k.write_xyz_trajectory')
    def test_extract_frame(self, mock_write, mock_conv, mock_read, engine_env):
        """Test _extract_frame method."""
        engine, _ = engine_env
        mock_read.return_value = [MagicMock()]
        mock_conv.return_value = ('box', 'xyz', 'vel', 'names')
        engine._extract_frame('traj.xyz', 0, 'out.xyz')
        mock_write.assert_called()

    @patch('pyretis.engines.cp2k.read_cp2k_restart')
    @patch('pyretis.engines.cp2k.read_cp2k_energy')
    @patch('pyretis.engines.cp2k.write_xyz_trajectory')
    @patch('pyretis.engines.cp2k.update_cp2k_input')
    def test_step(self, mock_upd, mock_write, mock_energy, mock_restart,
                  engine_env):
        """Test step method."""
        engine, _ = engine_env
        mock_restart.return_value = ('atoms', 'xyz', 'vel', 'box', None)
        mock_energy.return_value = {'ekin': [1.0], 'vpot': [2.0]}
        mock_system = MagicMock()

        vel = np.array([[0.1, 0.2, 0.3]])
        # CP2KEngine.step calls dump_frame, should return a string filename
        with patch.object(engine, 'dump_frame', return_value='conf.xyz'):
            with patch.object(engine, '_read_configuration',
                              return_value=('box', 'xyz', vel, 'atoms')):
                with patch.object(engine, 'run_cp2k',
                                  return_value={'restart': 'r',
                                                'energy': 'e'}):
                    res = engine.step(mock_system, 'name')
                    assert isinstance(res, str)
                    assert 'name.xyz' in res
                    mock_system.particles.set_pos.assert_called()

    @patch('pyretis.engines.cp2k.read_cp2k_restart')
    @patch('pyretis.engines.cp2k.read_cp2k_energy')
    @patch('pyretis.engines.cp2k.write_xyz_trajectory')
    @patch('pyretis.engines.cp2k.update_cp2k_input')
    def test_integrate(self, mock_upd, mock_write, mock_energy, mock_restart,
                       engine_env):
        """Test integrate method."""
        engine, _ = engine_env
        mock_restart.return_value = ('atoms', 'xyz', np.zeros((1, 3)),
                                     'box', None)
        mock_energy.return_value = {'ekin': [1.0, 1.1], 'vpot': [2.0, 2.1]}
        mock_ensemble = {'system': MagicMock(), 'order_function': MagicMock()}

        with patch.object(engine, 'dump_frame', return_value='frame.xyz'):
            with patch.object(engine, '_read_configuration',
                              return_value=('box', 'xyz', 'vel', 'atoms')):
                with patch.object(engine, 'run_cp2k',
                                  return_value={'restart': 'r', 'energy': 'e',
                                                'wfn': 'w'}):
                    with patch.object(engine, 'calculate_order',
                                      return_value=[0.1]):
                        with patch.object(engine, '_movefile'):
                            with patch.object(engine, '_remove_files'):
                                # The engine.integrate is a generator
                                res = list(engine.integrate(mock_ensemble, 2))
                                assert len(res) == 2
                                assert 'thermo' in res[0]

    @patch('pyretis.engines.cp2k.read_cp2k_restart')
    @patch('pyretis.engines.cp2k.read_cp2k_energy')
    @patch('pyretis.engines.cp2k.write_xyz_trajectory')
    @patch('pyretis.engines.cp2k.update_cp2k_input')
    def test_propagate_from(self, mock_upd, mock_write, mock_energy,
                            mock_restart, engine_env):
        """Test _propagate_from method."""
        engine, _ = engine_env
        mock_restart.return_value = ('atoms', 'xyz', 'vel', 'box', None)
        # return enough energy values for multiple steps
        mock_energy.return_value = {'ekin': [1.0]*20, 'vpot': [2.0]*20}

        mock_path = MagicMock()
        mock_path.maxlen = 2
        mock_system = MagicMock()
        # system.particles.get_pos() should return a tuple (filename, index)
        mock_system.particles.get_pos.return_value = ('initial.xyz', 0)
        mock_ensemble = {
            'system': mock_system,
            'interfaces': [0.1, 0.5, 0.9]
        }
        mock_msg = MagicMock()

        vel = np.array([[0.1, 0.2, 0.3]])
        with patch.object(engine, '_read_configuration',
                          return_value=('box', 'xyz', vel, 'atoms')):
            with patch.object(engine, 'run_cp2k',
                              return_value={'restart': 'r', 'energy': 'e',
                                            'wfn': 'w'}):
                with patch.object(engine, 'calculate_order',
                                  return_value=[0.2]):
                    with patch.object(engine, 'snapshot_to_system'):
                        with patch.object(engine, 'add_to_path',
                                          return_value=('status', True,
                                                        True, True)):
                            success, status = engine._propagate_from(
                                'name', mock_path, mock_ensemble, mock_msg
                            )
                            assert success is True
