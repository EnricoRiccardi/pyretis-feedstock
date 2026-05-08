# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the Simulation class."""
from io import StringIO
import os
import logging
import pytest
import tempfile
import pathlib
import numpy as np
from pyretis.core.particles import Particles
from pyretis.core.random_gen import RandomGenerator
from pyretis.core.units import units_from_settings
from pyretis.engines.internal import Langevin
from pyretis.inout.screen import ScreenOutput
from pyretis.inout.settings import add_default_settings
from pyretis.inout.simulationio import OutputTask
from pyretis.setup.createsystem import create_system
from pyretis.simulation.simulation import Simulation
from contextlib import contextmanager
from unittest.mock import patch
from .help import turn_on_logging, TEST_SETTINGS, TxtWriter
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


class TestSimulation:
    """Run the tests for Simulation class."""

    def test_add_task(self):
        """Test that we can add task(s) to the simulation."""
        simulation = Simulation({}, {})

        def task1():  # pylint: disable=missing-docstring
            return 'Hello there!'

        def task0():  # pylint: disable=missing-docstring
            return 'I should be first!'

        task = {'func': task1, 'result': 'hello'}
        add = simulation.add_task(task)
        assert add

        task = {'func': task0, 'result': 'first'}
        add = simulation.add_task(task, position=0)
        assert add

        task = {'func': 100, 'result': 'should-not-be-added'}
        add = simulation.add_task(task, position=0)
        assert not add

    def test_run_extend(self):
        """Test that we can run and extend a simulation."""
        simulation = Simulation({}, {'steps': 10})

        def task1(simulation):
            """Dummy task for the simulation."""
            return f"Hello there {simulation.cycle['step']:03d}"

        task = {'func': task1, 'result': 'hello', 'args': (simulation,)}
        add = simulation.add_task(task)
        assert add

        writer = ScreenOutput(
            TxtWriter(
                'test',
                header={'labels': ['Step', 'Message'], 'width': [9, 20]}
            )
        )
        simulation.output_tasks = [
            OutputTask('Hello-test', ['hello'], writer, {'every': 1}),
        ]

        lines = []
        with patch('sys.stdout', new=StringIO()) as stdout:
            for i, step in enumerate(simulation.run()):
                if i == 0:
                    assert 'hello' not in step
                else:
                    assert 'hello' in step
            simulation.extend_cycles(5)
            for i, step in enumerate(simulation.run()):
                assert i + 1 + 10 == step['cycle']['step']
            for linei in stdout.getvalue().strip().split('\n'):
                if linei:
                    lines.append(linei)
        for i, linei in enumerate(lines):
            if i == 0:
                assert '#    Step              Message' == linei
            else:
                msg = f'Hello there {i:03d}'
                correct = f'{i:>10d} {msg:>20s}'
                assert correct == linei

    def test_restart_simple(self):
        """Test restart methods."""
        simulation = Simulation(settings={}, controls={'steps': 10})
        for _ in simulation.run():
            pass
        restart = simulation.restart_info()
        simulation2 = Simulation(settings={}, controls={'steps': 0})
        assert simulation2.cycle['steps'] == 0
        simulation2.load_restart_info(restart)
        assert simulation2.cycle['stepno'] == 10
        sr1 = simulation.restart_info()
        sr2 = simulation2.restart_info()
        # Steps and end are not copied
        assert sr1['simulation']['cycle'] != sr2['simulation']['cycle']
        assert sr1['simulation']['type'] == sr2['simulation']['type']
        # Add random generator and try again:
        rgen = RandomGenerator(seed=101)
        simulation2.load_restart_info(restart)
        # Try adding an engine when we don't expect it:
        restart['engine'] = {'rgen': rgen.get_state()}
        simulation2.load_restart_info(restart)

        # Try setting fail system:
        simulation2.system = 'bumpy'
        with pytest.raises(AttributeError):
            simulation2.restart_info()

        # Add and try again:
        simulation2 = Simulation(settings={}, controls={'steps': 0})
        simulation2.engine = Langevin(1.0, 2.0, rgen=None)
        simulation2.load_restart_info(restart)

    def test_other_restart_features(self):
        """Test restart methods."""
        settings = TEST_SETTINGS.copy()
        add_default_settings(settings)
        simulation = Simulation(settings, controls={'steps': 0})
        simulation2 = Simulation(settings, controls={'steps': 0})
        units_from_settings(settings)
        system = create_system(settings)
        simulation.system = system
        simulation2.system = system
        particle = Particles(dim=16)
        particle.vpot = 101
        simulation.system.particles = particle

        rs1 = simulation.restart_info()
        simulation2.load_restart_info(rs1)

        assert np.allclose(simulation.system.particles.vpot,
                           simulation2.system.particles.vpot)

    def test_set_up_output(self, caplog):
        """Test the set_up_output things."""
        settings = {}
        add_default_settings(settings)
        simulation = Simulation(settings, {'steps': 100})
        with tempfile.TemporaryDirectory() as tempdir:
            settings['simulation']['exe_path'] = tempdir
            simulation.set_up_output(settings)
            assert (simulation.restart_freq ==
                    settings['output']['restart-file'])
            settings['output']['restart-file'] = -1
            with turn_on_logging():
                with caplog.at_level(logging.WARNING,
                                     logger='pyretis.simulation.simulation'):
                    simulation.set_up_output(settings)
                    assert simulation.restart_freq is None

    def test_soft_exit(self):
        """Test the soft exit method."""
        settings = TEST_SETTINGS.copy()
        with tempfile.TemporaryDirectory() as tempdir:
            simulation = Simulation({}, {'steps': 100})
            add_default_settings(settings)
            units_from_settings(settings)
            system = create_system(settings)
            simulation.system = system
            settings['simulation']['exe_path'] = tempdir
            simulation.set_up_output(settings)
            assert tempdir == simulation.exe_dir
            for i, _ in enumerate(simulation.run()):
                if i > simulation.restart_freq:
                    break
            files = [i.name for i in os.scandir(tempdir)]
            assert 1 == len(files)
            assert 'pyretis.restart' in files
            exit_file = os.path.join(tempdir, 'EXIT')
            pathlib.Path(exit_file).touch()
            with capture_log_output() as log_out:
                for _ in enumerate(simulation.run()):
                    pass
                assert 'soft exit' in log_out.getvalue().strip()
            assert (simulation.cycle['step'] ==
                    simulation.restart_freq + 2)
