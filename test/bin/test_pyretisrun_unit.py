# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Unit tests for the pyretisrun binary script."""
import os
from unittest.mock import MagicMock, patch, ANY
import pytest
from pyretis.bin.pyretisrun import (
    hello_world, bye_bye_world, use_tqdm,
    run_md_flux_simulation, run_md_simulation,
    explore_simulation, run_path_simulation,
    make_tis_files, run_generic_simulation,
    set_up_simulation, store_simulation_settings,
    main, soft_exit_ignore
)


@patch('pyretis.bin.pyretisrun.print_to_screen')
@patch('pyretis.bin.pyretisrun.logger')
def test_hello_world(mock_logger, mock_print):
    """Test hello_world function."""
    hello_world('input.rst', 'rundir', 'logfile.log')
    mock_print.assert_called()
    mock_logger.info.assert_called()


@patch('pyretis.bin.pyretisrun.print_to_screen')
@patch('pyretis.bin.pyretisrun.logger')
def test_bye_bye_world(mock_logger, mock_print):
    """Test bye_bye_world function."""
    bye_bye_world()
    mock_print.assert_called()
    mock_logger.info.assert_called()


def test_use_tqdm():
    """Test use_tqdm function."""
    # Test True
    pbar = use_tqdm(True)
    import tqdm
    assert pbar == tqdm.tqdm

    # Test False
    pbar = use_tqdm(False)
    assert pbar != tqdm.tqdm
    # Call the returned dummy
    res = pbar([1, 2, 3])
    assert list(res) == [1, 2, 3]
    res2 = pbar(iterable=[4, 5])
    assert list(res2) == [4, 5]


@patch('pyretis.bin.pyretisrun.use_tqdm')
def test_run_simulation_runners(mock_tqdm):
    """Test various simulation runners."""
    mock_sim = MagicMock()
    mock_sim.run.return_value = [1, 2, 3]
    mock_sim.cycle = {'startcycle': 0, 'endcycle': 10}
    mock_sim_settings = {
        'simulation': {'exe_path': 'some_path', 'task': 'md'},
        'ensemble': []
    }

    # Mock tqdm to just return the iterator
    mock_tqdm.return_value = lambda x, **kwargs: x

    # Test run_md_flux_simulation
    assert run_md_flux_simulation(mock_sim, mock_sim_settings)
    mock_sim.set_up_output.assert_called()
    mock_sim.write_restart.assert_called()

    # Test run_md_simulation
    mock_sim.reset_mock()
    assert run_md_simulation(mock_sim, mock_sim_settings)
    mock_sim.set_up_output.assert_called()
    mock_sim.write_restart.assert_called()

    # Test run_generic_simulation
    mock_sim.reset_mock()
    assert run_generic_simulation(mock_sim, mock_sim_settings)
    mock_sim.set_up_output.assert_called()
    mock_sim.write_restart.assert_called()


@patch('pyretis.bin.pyretisrun.use_tqdm')
def test_path_simulation_runners(mock_tqdm):
    """Test explore and path simulation runners."""
    mock_sim = MagicMock()
    mock_sim.run.return_value = [1, 2, 3]
    mock_sim.initiate.return_value = True
    mock_sim.cycle = {'startcycle': 0, 'endcycle': 10}
    mock_sim_settings = {
        'simulation': {'task': 'tis', 'exe_path': 'some_path'},
        'ensemble': [{'tis': {'freq': 0.5, 'allowmaxlength': False}}]
    }
    mock_tqdm.return_value = lambda x, **kwargs: x

    # Test explore_simulation
    assert explore_simulation(mock_sim, mock_sim_settings)
    assert mock_sim_settings['ensemble'][0]['tis']['freq'] == 0
    assert mock_sim_settings['ensemble'][0]['tis']['allowmaxlength'] is True

    # Test run_path_simulation
    mock_sim.reset_mock()
    assert run_path_simulation(mock_sim, mock_sim_settings)
    mock_sim.initiate.assert_called()

    # Test initiate failure
    mock_sim.initiate.return_value = False
    assert not run_path_simulation(mock_sim, mock_sim_settings)


@patch('pyretis.bin.pyretisrun.write_settings_file')
def test_make_tis_files(mock_write):
    """Test make_tis_files function."""
    settings = {
        'ensemble': [
            {'simulation': {'zero_ensemble': True},
             'engine': {'exe_path': 'p1'}},
            {'simulation': {'zero_ensemble': False},
             'engine': {'exe_path': 'p2'}}
        ],
        'simulation': {'zero_ensemble': False}
    }
    assert make_tis_files(None, settings)
    assert mock_write.call_count == 2


@patch('pyretis.bin.pyretisrun.parse_settings_file')
@patch('pyretis.bin.pyretisrun.create_simulation')
def test_set_up_simulation(mock_create, mock_parse):
    """Test set_up_simulation function."""
    mock_parse.return_value = {
        'simulation': {'task': 'md'},
        'engine': {},
        'ensemble': [{'simulation': {}, 'engine': {}}]
    }
    # Create a dummy file
    with patch('os.path.isfile', return_value=True):
        runner, sim, settings = set_up_simulation('input.rst', 'runpath')
        assert runner == run_md_simulation
        mock_create.assert_called()
        assert settings['simulation']['exe_path'] == 'runpath'
        assert settings['ensemble'][0]['simulation']['exe_path'] == 'runpath'


@patch('pyretis.bin.pyretisrun.write_settings_file')
def test_store_simulation_settings(mock_write):
    """Test store_simulation_settings function."""
    store_simulation_settings({}, 'indir', True)
    mock_write.assert_called_with({}, os.path.join('indir', 'out.rst'),
                                  backup=True)


@patch('pyretis.bin.pyretisrun.set_up_simulation')
@patch('pyretis.bin.pyretisrun.store_simulation_settings')
@patch('pyretis.bin.pyretisrun.soft_exit_ignore')
@patch('pyretis.bin.pyretisrun.bye_bye_world')
def test_main_func(mock_bye, mock_soft, mock_store, mock_setup):
    """Test the main function logic."""
    mock_run = MagicMock()
    mock_sim = MagicMock()
    mock_sim.cycle = {'step': 5}
    mock_setup.return_value = (mock_run, mock_sim, {'simulation': {}})

    with patch('os.path.isfile', return_value=False):  # No EXIT file
        main('infile', 'indir', 'exe_dir', False, 20)
        mock_setup.assert_called_with('infile', 'exe_dir')
        mock_run.assert_called()
        mock_bye.assert_called()


@patch('pyretis.bin.pyretisrun.print_to_screen')
@patch('pyretis.bin.pyretisrun.bye_bye_world')
def test_main_exit_file(mock_bye, mock_print):
    """Test main function when EXIT file exists."""
    with patch('os.path.isfile', return_value=True):
        main('infile', 'indir', 'exe_dir', False, 20)
        # Verify it prints error and calls bye_bye_world
        mock_print.assert_any_call(ANY, level='error')
        mock_bye.assert_called()


@patch('pyretis.bin.pyretisrun.set_up_simulation')
@patch('pyretis.bin.pyretisrun.store_simulation_settings')
@patch('pyretis.bin.pyretisrun.print_to_screen')
@patch('pyretis.bin.pyretisrun.logger')
def test_main_exception(mock_logger, mock_print, mock_store, mock_setup):
    """Test main function exception handling."""
    mock_setup.side_effect = Exception("Test Error")
    with patch('os.path.isfile', return_value=False):
        # We test with log_level = 20 (INFO) so it doesn't re-raise
        main('infile', 'indir', 'exe_dir', False, 20)
        mock_print.assert_any_call('ERROR - execution stopped.', level='error')

    # Test with DEBUG level to see it re-raises
    with patch('os.path.isfile', return_value=False):
        with pytest.raises(Exception, match="Test Error"):
            main('infile', 'indir', 'exe_dir', False, 10)  # 10 is DEBUG


def test_set_up_simulation_errors():
    """Test set_up_simulation error cases."""
    with patch('os.path.isfile', return_value=False):
        with pytest.raises(ValueError, match='NOT found'):
            set_up_simulation('missing.rst', 'path')


@patch('pyretis.bin.pyretisrun.use_tqdm')
def test_explore_simulation_init_fail(mock_tqdm):
    """Test explore_simulation initiation failure."""
    mock_sim = MagicMock()
    mock_sim.initiate.return_value = False
    assert not explore_simulation(mock_sim, {'ensemble': []})


def test_soft_exit_ignore():
    """Test soft_exit_ignore function."""
    with patch('signal.signal') as mock_signal:
        soft_exit_ignore(True, 'some_dir')
        mock_signal.assert_called()
        soft_exit_ignore(False, 'some_dir')
        mock_signal.assert_called()
