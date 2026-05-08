# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the path simulation classes."""
from io import StringIO
import copy
import logging
import os
import tempfile
import pytest
import pathlib

from pyretis.inout.settings import (add_default_settings,
                                    fill_up_tis_and_retis_settings,)
from pyretis.setup.createsimulation import create_ensembles
from pyretis.simulation.path_simulation import (PathSimulation,
                                                SimulationTIS,
                                                SimulationRETIS,)
from contextlib import contextmanager
from unittest.mock import patch
from .help import TEST_SETTINGS

logging.disable(logging.CRITICAL)


@contextmanager
def capture_log_output():
    """Capture all logger output to a StringIO buffer."""
    stream = StringIO()
    handler = logging.StreamHandler(stream)
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter('%(message)s'))
    root_logger = logging.getLogger('')
    root_logger.addHandler(handler)
    prev_disable = logging.root.manager.disable
    prev_level = root_logger.level
    logging.disable(logging.NOTSET)
    root_logger.setLevel(logging.DEBUG)
    try:
        yield stream
    finally:
        root_logger.removeHandler(handler)
        root_logger.setLevel(prev_level)
        logging.disable(prev_disable)


HERE = os.path.abspath(os.path.dirname(__file__))

SETTINGS1 = {
    'simulation': {
        'task': 'tis',
        'steps': 10,
        'interfaces': [-0.9, -0.9, 1.0],
        'seed': 1,
        'rgen': 'rgen',
    },
    'system': {'units': 'lj', 'temperature': 0.1, 'dimensions': 1},
    'engine': {
        'class': 'Langevin',
        'gamma': 0.3,
        'high_friction': False,
        'seed': 321,
        'timestep': 0.2
    },
    'potential': [{'a': 1.0, 'b': 2.0, 'c': 0.0, 'class': 'DoubleWell'}],
    'orderparameter': {
        'class': 'Position',
        'dim': 'x',
        'index': 0,
        'name': 'Position',
        'periodic': False
    },
    'initial-path': {'method': 'kick'},
    'particles': {'position': {'input_file': 'initial.xyz'}},
    'retis': {
        'swapfreq': 0.5,
        'nullmoves': True,
        'swapsimul': True,
    },
    'tis': {
        'freq': 0.5,
        'maxlength': 2000,
        'aimless': True,
        'allowmaxlength': False,
        'sigma_v': -1,
        'seed': 1,
        'rgen': 'rgen',
        'zero_momentum': False,
        'rescale_energy': False,
    },
}


def create_variables(exe_path=None):
    """Create system, engine, order parameter etc."""
    settings = {key: copy.deepcopy(val) for key, val in TEST_SETTINGS.items()}
    settings['engine'] = {'class': 'Langevin', 'gamma': 0.3,
                          'high_friction': False,
                          'timestep': 0.2}
    settings['orderparameter'] = {'class': 'Position',
                                  'index': 0}

    if exe_path is not None:
        settings['simulation'] = {'exe_path': exe_path}
    add_default_settings(settings)
    settings['particles']['type'] = 'internal'
    settings['simulation']['interfaces'] = [-1.0, 0.0, 1.0, 2.0]
    settings['simulation']['task'] = 'tis'
    settings['simulation']['zero_ensemble'] = True
    settings['simulation']['zero_left'] = -16.0
    settings['simulation']['flux'] = True
    settings['retis'] = {'rgen': 'rgen'}
    with patch('sys.stdout', new=StringIO()):
        fill_up_tis_and_retis_settings(settings)
        ensembles = create_ensembles(settings)
    return settings, ensembles


def create_variables1(exe_path=None):
    """Create system, engine, order parameter etc."""
    settings = {key: copy.deepcopy(val) for key, val in SETTINGS1.items()}
    add_default_settings(settings)
    settings['particles']['position']['input_file'] =\
        os.path.join(HERE, 'initial.xyz')
    settings['particles']['type'] = 'internal'
    settings['simulation']['task'] = 'retis'
    settings['simulation']['zero_ensemble'] = True
    settings['simulation']['flux'] = True
    settings['simulation']['exe_path'] = exe_path

    with patch('sys.stdout', new=StringIO()):
        fill_up_tis_and_retis_settings(settings)
        ensembles = create_ensembles(settings)
    return settings, ensembles


class TestPathSimulation:
    """Run the tests for the PathSimulation class."""

    def test_init(self):
        """Test that we can create the simulation object."""
        settings, ensembles = create_variables()
        simulation = PathSimulation(ensembles, settings, controls={})
        assert hasattr(simulation, 'rgen')
        settings['engine']['sigma_v'] = 0.1
        simulation = PathSimulation(ensembles, settings, controls={})
        assert len(simulation.ensembles) == len(ensembles)
        for i, j in zip(simulation.ensembles, ensembles):
            assert i is j
        assert simulation.settings['tis']['aimless']
        settings['simulation']['task'] = 'make-tis-files'
        del settings['ensemble']
        fill_up_tis_and_retis_settings(settings)
        for ens_set in settings['ensemble']:
            assert 'tis' == ens_set['simulation']['task']
        del settings['tis']
        with pytest.raises(ValueError):
            PathSimulation(ensembles, settings, controls={})

    def test_init_sigma_v(self):
        """Test that we can create the simulation object."""
        settings, ensembles = create_variables()
        settings['ensemble'][1]['tis']['sigma_v'] = 0.1
        simulation = PathSimulation(ensembles, settings, controls={})
        assert not simulation.settings['ensemble'][1]['tis']['aimless']

    def test_restart_info(self):
        """Test generation of restart info."""
        settings, ensembles = create_variables()
        simulation = PathSimulation(ensembles,
                                    settings,
                                    controls={'steps': 16})
        info = simulation.restart_info()
        for key in ('cycle', 'rgen'):
            assert key in info['simulation']
        assert info['simulation']['rgen']['state'][2] == 624
        assert 'system' in info
        simulation = PathSimulation(ensembles, settings, controls={})
        simulation.load_restart_info(info)
        assert simulation.rgen.get_state()['state'][2] == 624
        assert simulation.cycle['startcycle'] == 16

    def test_write_restart_file(self):
        """Test that we can create the restart files."""
        with tempfile.TemporaryDirectory() as tempdir:
            settings, ensembles = create_variables(exe_path=tempdir)
            simulation = PathSimulation(ensembles, settings, controls={})
            simulation.set_up_output(settings)
            simulation.write_restart(now=True)
            dirs = [i.name for i in os.scandir(tempdir) if i.is_dir()]
            assert len(dirs) == 4
            for i in ensembles:
                assert i['path_ensemble'].ensemble_name_simple in dirs
                for j in os.scandir(
                        i['path_ensemble'].directory['path_ensemble']):
                    if j.is_file():
                        assert 'ensemble.restart' == j.name
            files = [i.name for i in os.scandir(tempdir) if i.is_file()]
            assert len(files) == 1
            assert 'pyretis.restart' in files

    def test_initiation(self):
        """Test the initiation method."""
        with tempfile.TemporaryDirectory() as tempdir:
            # Generate a input file in the temp directory
            with open(os.path.join(tempdir, 'initial.xyz'), 'w',
                      encoding='utf-8') as scrivi:
                scrivi.write('1\n \n Ar  -1.0  -0.0  -0.0 ')

            settings, ensembles = create_variables1(exe_path=tempdir)
            simulation = PathSimulation(ensembles, settings, controls={})
            simulation.set_up_output(settings)
            with patch('sys.stdout', new=StringIO()):
                assert simulation.initiate(settings)
            # Check the soft exit option:
            open(os.path.join(tempdir, 'EXIT'), 'w', encoding='utf-8').close()
            with patch('sys.stdout', new=StringIO()):
                assert not simulation.initiate(settings)


class TestSimulationRETIS:
    """Run the tests for the SimulationRETIS class."""

    def test_step(self):
        """Test the initiation method."""
        moves = {'move-0': 'tr', 'move-1': 'tr', 'move-2': 'tr'}

        with tempfile.TemporaryDirectory() as tempdir:
            settings, ensembles = create_variables1(exe_path=tempdir)
            simulation = SimulationRETIS(ensembles, settings, controls={})
            simulation.set_up_output(settings)
            with patch('sys.stdout', new=StringIO()):
                assert simulation.initiate(settings)
                result = simulation.step()
                for key, move in moves.items():
                    assert key in result
                    assert move == result[key]

    def test_soft_exit(self):
        """Test the soft exit method."""
        with tempfile.TemporaryDirectory() as tempdir:
            settings, ensembles = create_variables1(exe_path=tempdir)
            simulation = SimulationRETIS(ensembles, settings, controls={})
            settings['simulation']['exe_path'] = tempdir
            simulation.set_up_output(settings)
            simulation.cycle['endcycle'] = 1
            exit_file = os.path.join(tempdir, 'EXIT')
            with capture_log_output() as log_out:
                assert simulation.initiate(settings)
                pathlib.Path(exit_file).touch()
                for _ in enumerate(simulation.run()):
                    pass
                assert 'soft exit' in log_out.getvalue().strip()
                restart_file = os.path.join(tempdir, 'pyretis.restart')
                assert os.path.exists(restart_file)
                assert simulation.cycle['step'] == 1

    def test_priority(self):
        """Test the priority shooting feature."""
        with tempfile.TemporaryDirectory() as tempdir:
            settings, ensembles = create_variables1(exe_path=tempdir)
            fill_up_tis_and_retis_settings(settings)
            simulation = SimulationRETIS(ensembles, settings, controls={})
            simulation.set_up_output(settings)
            with patch('sys.stdout', new=StringIO()):
                assert simulation.initiate(settings)
                simulation.cycle['step'] = 2
                settings['retis']['swapfreq'] = 0.0
                settings['simulation']['priority_shooting'] = True
                simulation.ensembles[0]['path_ensemble'].nstats['npath'] = 2
                simulation.ensembles[1]['path_ensemble'].nstats['npath'] = 1
                simulation.ensembles[2]['path_ensemble'].nstats['npath'] = 1

                # Get nstats['npath'] = 2 for all path_ensembles:
                result = simulation.step()
                assert 'move-0' not in result
                assert 'move-1' in result
                assert 'move-2' in result
                assert [2, 2, 2] == [i['path_ensemble'].nstats['npath']
                                     for i in simulation.ensembles]

                # test case when we 100% do swap with unequal nstats['npath']:
                simulation.ensembles[2]['path_ensemble'].nstats['npath'] = 0
                settings['retis']['swapfreq'] = 1.0

                # 1. We do TIS move in ensembles[2] only:
                result = simulation.step()
                assert 'move-0' not in result
                assert 'move-2' in result

                # We do swap move in when priority_shooting = False:
                settings['simulation']['priority_shooting'] = False
                result = simulation.step()
                for i in range(3):
                    move = 'swap' if i else 'nullmove'
                    assert move == result.get('move-' + str(i))

    def test_priority_tis(self):
        """Test the priority shooting feature for tis."""
        with tempfile.TemporaryDirectory() as tempdir:
            settings, _ = create_variables1(exe_path=tempdir)
            settings['simulation']['task'] = 'tis'
            settings['simulation']['interfaces'] = [-0.9, -0.8, -0.7, 1.0]
            settings['simulation']['flux'] = False
            del settings['ensemble']
            fill_up_tis_and_retis_settings(settings)
            with patch('sys.stdout', new=StringIO()):
                ensembles = create_ensembles(settings)
                simulation = SimulationTIS(ensembles, settings, controls={})
                simulation.set_up_output(settings)
                assert simulation.initiate(settings)
                simulation.cycle['endcycle'] = 4
                simulation.cycle['step'] = 2
                settings['simulation']['priority_shooting'] = True
                simulation.ensembles[0]['path_ensemble'].nstats['npath'] = 2
                simulation.ensembles[1]['path_ensemble'].nstats['npath'] = 1
                simulation.ensembles[2]['path_ensemble'].nstats['npath'] = 1

                cnt = 0
                while not simulation.is_finished():
                    simulation.step()
                    cnt += 1
                assert [4, 4, 4] == [i['path_ensemble'].nstats['npath']
                                     for i in simulation.ensembles]
                assert cnt == 3
                assert simulation.is_finished()
