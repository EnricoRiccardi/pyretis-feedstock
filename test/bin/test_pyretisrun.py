# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the bin."""
import logging
import os
import pytest
import shutil
import subprocess
import types
import tempfile
import tqdm  # For a progress bar
from contextlib import contextmanager
from io import StringIO
from unittest.mock import patch
from collections.abc import Iterable
from pyretis.bin.pyretisrun import (hello_world,
                                    use_tqdm,
                                    bye_bye_world,
                                    set_up_simulation,
                                    run_generic_simulation,
                                    main,
                                    _RUNNERS)
from pyretis.setup import create_simulation
from pyretis.inout.settings import parse_settings_file

logging.disable(logging.CRITICAL)
HERE = os.path.abspath(os.path.dirname(__file__))


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


class TestPyretisMainRunner:
    """Test Main."""

    def test_pyretisrun(self):
        """Test pyretisrun."""
        with patch('sys.stdout', new=StringIO()) as stdout:
            with subprocess.Popen(['pyretisrun', '-V'],
                                  stdout=subprocess.PIPE) as ex:
                asd = ex.stdout.read().split()
            assert b'PyRETIS' in asd

            with subprocess.Popen(['pyretisrun'],
                                  stderr=subprocess.PIPE) as ex:
                asd = ex.stderr.read().split()
            assert b'pyretisrun' in asd

    def test_pyretisrun2(self):
        """Test pyretisrun by actually running a retis sim with 0 steps."""
        input_file = 'dummy_input.rst'
        with patch('sys.stdout', new=StringIO()):
            with tempfile.TemporaryDirectory() as tempdir:
                shutil.copyfile(os.path.join(HERE, input_file),
                                os.path.join(tempdir, input_file))
                with subprocess.Popen(['pyretisrun', '-i', input_file, '-p',
                                       '-l', 'DEBUG'],
                                      cwd=tempdir,
                                      stderr=subprocess.DEVNULL,
                                      stdout=subprocess.PIPE) as ex:
                    asd = ex.stdout.read().split()
        assert b'Riccardi' in asd
        assert b'"kick"' in asd

    def test_exit(self):
        """Test pyretisrun exit by creating an EXIT file."""
        infile = 'dummy_input.rst'
        with tempfile.TemporaryDirectory() as tempdir:
            input_file = os.path.join(tempdir, infile)
            shutil.copyfile(os.path.join(HERE, infile), input_file)
            with capture_log_output() as log_out:
                open(tempdir + '/EXIT', 'w').close()
                main(input_file, tempdir, tempdir,
                     progress=False, log_level=20)
            assert 'EXIT file found' in log_out.getvalue().strip()

    def test_main(self):
        """Test the main() function for a tis task."""
        infile = 'dummy_input.rst'
        with tempfile.TemporaryDirectory() as tempdir:
            input_file = os.path.join(tempdir, infile)
            shutil.copyfile(os.path.join(HERE, infile), input_file)
            with capture_log_output() as log_out:
                main(input_file, tempdir, tempdir, progress=False, log_level=9)
            assert 'Execution ended' in log_out.getvalue().strip()
            with capture_log_output() as log_out:
                input_file = os.path.join(tempdir, 'does_not_exist.rst')
                pytest.raises(ValueError, main, input_file, tempdir,
                              tempdir, progress=False,
                              log_level=9)
            assert 'ERROR - execution stopped.' in log_out.getvalue().strip()

    def test_run_simulation(self):
        """Test all simulations functions."""
        infile = 'dummy_input.rst'
        runners = _RUNNERS
        runners['umbrellawindow'] = run_generic_simulation
        for key, runner in runners.items():
            with patch('sys.stdout', new=StringIO()):
                with tempfile.TemporaryDirectory() as tempdir:
                    if key in ['repptis', 'retis']:
                        continue
                    input_file = os.path.join(tempdir, infile)
                    shutil.copyfile(os.path.join(HERE, infile), input_file)
                    sim_settings = parse_settings_file(input_file)
                    sim_settings['simulation']['exe_path'] = tempdir
                    sim_settings['engine']['exe_path'] = tempdir
                    for ens in sim_settings['ensemble']:
                        ens['simulation']['exe_path'] = tempdir
                        ens['engine']['exe_path'] = tempdir
                    sim_settings['simulation']['task'] = key
                    sim = create_simulation(sim_settings)
                    assert runner(sim, sim_settings)

    def test_simulation_exit(self):
        """Test exit feature of explore_simulation and run_path_simulation."""
        infile = 'dummy_input.rst'
        runners = _RUNNERS
        for key, runner in runners.items():
            with patch('sys.stdout', new=StringIO()):
                with tempfile.TemporaryDirectory() as tempdir:
                    if key not in ('tis', 'explore'):
                        continue
                    input_file = os.path.join(tempdir, infile)
                    shutil.copyfile(os.path.join(HERE, infile), input_file)
                    sim_settings = parse_settings_file(input_file)
                    sim_settings['simulation']['exe_path'] = tempdir
                    sim_settings['engine']['exe_path'] = tempdir
                    for ens in sim_settings['ensemble']:
                        ens['simulation']['exe_path'] = tempdir
                        ens['engine']['exe_path'] = tempdir
                    sim_settings['simulation']['task'] = key
                    sim = create_simulation(sim_settings)
                    open(tempdir + '/EXIT', 'w').close()
                    assert not runner(sim, sim_settings)

    def test_use_tqdm(self):
        """Test tqdm."""
        bar = use_tqdm(progress=True)
        assert isinstance(bar, type(tqdm.tqdm))
        assert bar is tqdm.tqdm
        assert not isinstance(bar, Iterable)

        bar = use_tqdm(progress=False)
        assert bar is not tqdm.tqdm
        assert isinstance(bar, types.FunctionType)
        value = bar(iterable='something1', dummy='something2')
        assert value == 'something1'

    def test_hello_world(self):
        """Test that we are polite."""
        with capture_log_output() as log_out:
            hello_world(infile='I_can_read_your_mind.rst',
                        rundir=HERE,
                        logfile='nothing_to_declare')
        assert 'Start' in log_out.getvalue().strip()

    def test_bye_world(self):
        """Test that we can die."""
        with capture_log_output() as log_out:
            bye_bye_world()
        assert 'reference' in log_out.getvalue().strip()

    def test_set_up_simulation(self):
        """Test that we know how to set up a simulation."""
        inputfile = 'a_non_existent_input.rst'
        with pytest.raises(ValueError) as err:
            set_up_simulation(inputfile, HERE)
        assert f'Input file "{inputfile}" NOT found!' == str(err.value)

        inputfile = os.path.join(HERE, 'dummy_input.rst')
        with capture_log_output() as log_out:
            set_up_simulation(inputfile, HERE)
        assert 'Reading' in log_out.getvalue().strip()
