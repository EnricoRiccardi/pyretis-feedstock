# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the simulationio module."""
from io import StringIO
import logging
import os
import pytest
from unittest.mock import patch
from pyretis.inout.fileio import FileIO
from pyretis.inout.archive import PathStorage
from pyretis.inout.formats import (
    PathIntFormatter,
    ThermoTableFormatter,
)
from pyretis.inout.screen import ScreenOutput
from pyretis.inout.simulationio import (
    get_task_type,
    get_file_mode,
    OutputTask,
    OUTPUT_TASKS,
    Task,
    task_from_settings,
)
from pyretis.engines.external import ExternalMDEngine
from pyretis.engines.internal import Verlet
from pyretis.inout.settings import add_default_settings
from pyretis.simulation import (
    Simulation,
    UmbrellaWindowSimulation,
    SimulationMD,
    SimulationNVE,
    SimulationMDFlux,
    SimulationTIS,
    SimulationRETIS,
)
logging.disable(logging.CRITICAL)


class TestTask:
    """Run the tests for :py:class`.Task`"""

    def test_when_change(self):
        """Test that we can change the "when" property."""
        task = Task(None)
        assert task.when is None
        task.when = None
        assert task.when is None
        task.when = {'every': 1}
        assert task.when == {'every': 1}
        task.when = None
        assert task.when is None
        task.when = {'every': 1, 'all-the-time': 100}
        assert task.when == {'every': 1}

    def test_task_dict(self):
        """Test that the dictionary returned from a task is ok."""
        task = Task(None)
        task.when = {'every': 2}
        task_dict = task.task_dict()
        assert task_dict['when'] == {'every': 2}

    def test_execute_now(self):
        """Test the execute now method."""
        step = {'step': 10, 'start': 2, 'stepno': 8}
        task = Task(None)

        task.when = None
        exe = task.execute_now(step)
        assert exe

        task.when = {'every': 2}
        exe = task.execute_now(step)
        assert exe

        task.when = {'every': 3}
        exe = task.execute_now(step)
        assert not exe

        task.when = {'at': 10}
        exe = task.execute_now(step)
        assert exe

        task.when = {'at': (9, 10)}
        exe = task.execute_now(step)
        assert exe


class DummyExternal(ExternalMDEngine):
    """A dummy external engine. Only useful for testing."""

    def __init__(self, input_path, timestep, subcycles):
        """Initialise the dummy engine."""
        super().__init__('External engine for testing!', timestep,
                         subcycles)
        self.input_path = input_path

    def step(self, system, name):
        """Perform a single dummy step."""

    def _read_configuration(self, filename):
        """Read a xyz configuration."""

    def _reverse_velocities(self, filename, outfile):
        """Reverse velocities with a dummy method."""

    def _extract_frame(self, traj_file, idx, out_file):
        """Extract a frame, dummy method."""

    def _propagate_from(self, name, path, system, order_function, interfaces,
                        msg_file, reverse=False):
        """Propagate with the engine, dummy method."""
        # pylint: disable=too-many-arguments

    def modify_velocities(self, system, rgen, sigma_v=None, aimless=True,
                          momentum=False, rescale=None):
        """Modify velocities, dummy method."""
        # pylint: disable=too-many-arguments


class TestMethods:
    """Run the tests for methods defined in the module."""

    def test_get_task_type(self):
        """Test the get_task_type method."""
        task = {'type': 'pathensemble'}
        task_type = get_task_type(task, None)
        assert task_type == task['type']
        task = {'type': 'path-traj-{}'}
        task_type = get_task_type(task, None)
        assert task_type == 'path-traj-int'
        engine = DummyExternal('test', 0.0, 1)
        task_type = get_task_type(task, engine)
        assert task_type == 'path-traj-ext'

    def test_get_file_mode(self):
        """Test the get_file_mode method."""
        settings = {'output': {'backup': 'overwrite'}}
        mode = get_file_mode(settings)
        assert mode == 'w'

        settings = {'output': {'backup': 'append'}}
        mode = get_file_mode(settings)
        assert mode == 'a'

        settings = {'output': {'backup': 10}}
        mode = get_file_mode(settings)
        assert mode == 'w'
        assert settings['output']['backup'] == 'backup'

    def test_task_from_settings(self):
        """Test that we can create tasks from settings."""
        for engine in (None, DummyExternal('test', 0.1, 1)):
            for key, val in OUTPUT_TASKS.items():
                task_dict = {'type': key, 'name': f'test-{key}'}
                settings = {}
                add_default_settings(settings)
                directory = os.path.join('this', 'is', 'a', 'fake', 'path')
                task = task_from_settings(task_dict, settings, directory,
                                          engine)
                assert task.result == val['result']
                assert task.target == val['target']
                assert task.when['every'] == settings['output'][val['when']]
                if val['target'] == 'file':
                    assert isinstance(task.writer, FileIO)
                elif val['target'] == 'screen':
                    assert isinstance(task.writer, ScreenOutput)
                elif val['target'] == 'file-archive':
                    assert isinstance(task.writer, PathStorage)
                if 'formatter' in val:
                    assert isinstance(task.writer.formatter, val['formatter'])

        task_dict = {'type': 'energy', 'name': 'test-energy'}
        settings = {}
        add_default_settings(settings)
        settings['output']['energy-file'] = 0
        task = task_from_settings(task_dict, settings, '.', None)
        assert task is None

        task_dict = {'type': 'path-traj-{}', 'name': 'test-path-traj'}
        settings = {}
        add_default_settings(settings)
        task = task_from_settings(task_dict, settings, '.', Verlet(0.1))
        assert isinstance(task.writer.formatter, PathIntFormatter)

        task = task_from_settings(task_dict, settings, '.',
                                  DummyExternal('test', 0.1, 1))
        assert isinstance(task.writer, PathStorage)

        task_dict = {'type': 'does-not-exist', 'name': 'test-energy'}
        settings = {}
        add_default_settings(settings)
        with pytest.raises(KeyError):
            task_from_settings(task_dict, settings, '.', None)

        OUTPUT_TASKS['thermo-web'] = {
            'target': 'web',
            'filename': 'thermo.txt',
            'result': ('thermo',),
            'when': 'energy-file',
            'formatter': ThermoTableFormatter,
        }
        task_dict = {'type': 'thermo-web', 'name': 'test-web'}
        settings = {}
        add_default_settings(settings)
        task = task_from_settings(task_dict, settings, '.', None)
        assert task is None
        task_dict = {'type': 'thermo-screen', 'name': 'thermo-screen'}
        task = task_from_settings(task_dict, settings, '.', None,
                                  progress=False)
        assert task is not None
        task = task_from_settings(task_dict, settings, '.', None,
                                  progress=True)
        assert task is None

    def test_create_tasks(self):
        """Test that we can create tasks from settings for simulations."""
        settings = {}
        add_default_settings(settings)
        for simulation in (UmbrellaWindowSimulation, SimulationMD,
                           SimulationNVE, SimulationMDFlux,
                           SimulationTIS, SimulationRETIS):
            sim = Simulation({}, {})
            sim.simulation_output = simulation.simulation_output
            sim.create_output_tasks(settings, progress=False)
            for task in sim.output_tasks:
                assert isinstance(task, OutputTask)


class TestOutputTask:
    """Test the functionality of the output task class."""

    def test_task_dict(self):
        """Test the task dict output."""
        formatter = ThermoTableFormatter()
        writer = ScreenOutput(formatter)
        task = OutputTask('thermo-test', ('thermo', ), writer, {'every': 1})
        task_dict = task.task_dict()
        assert task_dict['writer'] == ScreenOutput
        assert task_dict['formatter'] == ThermoTableFormatter

    def test_output(self):
        """Test that we can output from the task."""
        formatter = ThermoTableFormatter()
        writer = ScreenOutput(formatter)
        task = OutputTask('thermo-test', ('thermo', ), writer, {'every': 1})
        result = {
            'cycle': {'step': 110, 'stepno': 10},
            'thermo': {'vpot': 1.0, 'ekin': 2.0,
                       'etot': 3.0, 'temp': 4.0, 'press': 5.0},
        }
        correct = [
            ('#     Step          Temp           Pot           '
             'Kin           Tot         Press'),
            ('       110             4             1             2'
             '             3             5'),
            ('      1000             4             1             2'
             '             3             5'),
        ]
        idx = 0
        with patch('sys.stdout', new=StringIO()) as stdout:
            task.output(result)
            result['cycle']['step'] = 1000
            task.output(result)
            for linei in stdout.getvalue().strip().split('\n'):
                if linei:
                    assert linei == correct[idx]
                    idx += 1
        task = OutputTask('thermo-test', ('thermo', ), writer, {'every': 99})
        out = task.output(result)
        assert not out
        result['cycle'] = {'step': 99, 'stepno': 99}
        del result['thermo']
        out = task.output(result)
        assert not out
