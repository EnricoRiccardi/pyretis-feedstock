# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the GromacsEngine class."""
import logging
import os
import tempfile
import shutil
import struct
import pytest
from pyretis.core.path import Path
from pyretis.engines.gromacs2 import (GromacsEngine2,
                                      GromacsRunner,
                                      reopen_file,
                                      read_remaining_trr)
from pyretis.inout.common import make_dirs
from pyretis.orderparameter.orderparameter import Position
from .test_helpers.test_helpers import make_test_system


logging.disable(logging.CRITICAL)
HERE = os.path.abspath(os.path.dirname(__file__))
GMX_DIR = os.path.join(HERE, 'gmx_input2')
GMX = os.path.join(HERE, 'mockgmx.py')
MDRUN = os.path.join(HERE, 'mockmdrun2.py')


class TestGromacsEngine2:
    """Run the tests for the GromacsEngine."""

    def test_init(self):
        """Test the initiation."""
        eng = GromacsEngine2(gmx='echo',
                             mdrun='echo',
                             input_path=GMX_DIR,
                             timestep=0.002,
                             subcycles=10,
                             maxwarn=10,
                             gmx_format='g96',
                             write_vel=True,
                             write_force=False)
        eng.exe_dir = GMX_DIR
        with pytest.raises(ValueError):
            GromacsEngine2(gmx='echo',
                           mdrun='echo',
                           input_path='gmx_input',
                           timestep=0.002,
                           subcycles=10,
                           gmx_format='not-a-format')
        with pytest.raises(FileNotFoundError):
            GromacsEngine2(gmx='echo',
                           mdrun='echo',
                           input_path='missing-files',
                           timestep=0.002,
                           subcycles=10)

    def test_propagate_forward(self):
        """Test the propagate method, forward direction.

        Here, the .trr file will be created before starting the
        engine. The engine will then be lagging behind and this
        should trigger the `read_remaining_trr` method.
        """
        with tempfile.TemporaryDirectory() as tempdir:
            eng = GromacsEngine2(GMX, MDRUN, GMX_DIR, 0.002, 7,
                                 maxwarn=1, gmx_format='g96',
                                 write_vel=True,
                                 write_force=False)
            rundir = os.path.join(tempdir, 'generate2gmxf')
            # Create the directory for running:
            make_dirs(rundir)
            eng.exe_dir = rundir
            # Create the system:
            system = make_test_system((eng.input_files['conf'], 0))
            # Propagate:
            orderp = Position(0, dim='x', periodic=False)
            path = Path(None, maxlen=8)
            ensemble = {'system': system, 'order_function': orderp,
                        'interfaces': [-0.45, 10.0, 14.0]}
            success, _ = eng.propagate(path, ensemble, reverse=False)
            assert success
            assert path.length == 4
            initial_x = -0.422
            for i, point in enumerate(path.phasepoints):
                assert point.particles.ekin == pytest.approx(
                    eng.subcycles * i)
                assert point.particles.vpot == pytest.approx(
                    eng.subcycles * -1.0 * i)
                assert point.order[0] == pytest.approx(
                    i * eng.subcycles + initial_x, abs=1e-3)
                assert not point.particles.get_vel()

            eng.clean_up()

    def test_propagate_backward(self):
        """Test the propagate method, backward direction.

        Here, the .trr file will be created before starting the
        engine. The engine will then be lagging behind and this
        should trigger the `read_remaining_trr` method.
        """
        with tempfile.TemporaryDirectory() as tempdir:
            eng = GromacsEngine2(GMX, MDRUN, GMX_DIR, 0.002, 7,
                                 maxwarn=1, gmx_format='g96',
                                 write_vel=True,
                                 write_force=False)
            rundir = os.path.join(tempdir, 'generate2gmxb')
            # Create the directory for running:
            make_dirs(rundir)
            eng.exe_dir = rundir
            # Create the system:
            system = make_test_system((eng.input_files['conf'], 0))
            # Propagate:
            orderp = Position(0, dim='x', periodic=False)
            path = Path(None, maxlen=8)
            ensemble = {'system': system, 'order_function': orderp,
                        'interfaces': [-20., 10.0, 14.0]}
            success, _ = eng.propagate(path, ensemble, reverse=True)
            assert success
            assert path.length == 4
            initial_x = -0.422
            for i, point in enumerate(path.phasepoints):
                assert point.particles.ekin == pytest.approx(eng.subcycles * i)
                assert (point.particles.vpot ==
                        pytest.approx(eng.subcycles * -1.0 * i))
                # (point.order[0] was already broken in previous attempt)
                assert (point.order[0] ==
                        pytest.approx(initial_x - i * eng.subcycles, abs=1e-3))
                assert point.particles.get_vel()
            eng.clean_up()

    def test_propagate_crash(self):
        """Test the propagate method when engine crashes."""
        with tempfile.TemporaryDirectory() as tempdir:
            mdrun = f'{MDRUN} -crash'
            eng = GromacsEngine2(GMX, mdrun, GMX_DIR,
                                 timestep=0.002,
                                 subcycles=7,
                                 maxwarn=1, gmx_format='g96',
                                 write_vel=True,
                                 write_force=False)
            rundir = os.path.join(tempdir, 'generate3gmxf')
            # Create the directory for running:
            make_dirs(rundir)
            eng.exe_dir = rundir
            # Create the system:
            system = make_test_system((eng.input_files['conf'], 0))
            # Propagate:
            order_function = Position(0, dim='x', periodic=False)
            path = Path(None, maxlen=8)
            with pytest.raises(RuntimeError):
                ensemble = {'system': system, 'order_function': order_function,
                            'interfaces': [-0.45, 10.0, 14.0]}
                eng.propagate(path, ensemble, reverse=False)
            # Check the error - output:
            with open(os.path.join(rundir, 'engine.err'), 'r',
                      encoding="utf8") as infile:
                data = infile.read()
                assert 'Crash error for testing.' in data
            eng.clean_up()

    def test_propagate_sleep(self):
        """Test the propagate method.

        Here, we try to write the trr file a bit slower so that
        the class running the GROMACS simulation will have to wait
        for the data to be written.
        """
        with tempfile.TemporaryDirectory() as tempdir:
            mdrun = f'{MDRUN} -sleep'
            eng = GromacsEngine2(GMX, mdrun, GMX_DIR, 0.002, 7,
                                 maxwarn=1, gmx_format='g96',
                                 write_vel=True,
                                 write_force=False)
            rundir = os.path.join(tempdir, 'generate4gmxf')
            # Create the directory for running:
            make_dirs(rundir)
            eng.exe_dir = rundir
            # Create the system:
            system = make_test_system((eng.input_files['conf'], 0))
            # Propagate:
            order_function = Position(0, dim='x', periodic=False)
            path = Path(None, maxlen=8)
            ensemble = {'system': system, 'order_function': order_function,
                        'interfaces': [-0.45, 10.0, 14.0]}
            success, _ = eng.propagate(path, ensemble, reverse=False)
            assert success
            assert path.length == 4
            initial_x = -0.422
            for i, point in enumerate(path.phasepoints):
                assert point.particles.ekin == pytest.approx(
                    eng.subcycles * i)
                assert point.particles.vpot == pytest.approx(
                    eng.subcycles * -1.0 * i)
                assert point.order[0] == pytest.approx(
                    i * eng.subcycles + initial_x, abs=1e-3)
            eng.clean_up()

    def test_integrate(self):
        """Test the integrate method."""
        with tempfile.TemporaryDirectory() as tempdir:
            eng = GromacsEngine2(GMX, MDRUN, GMX_DIR, 0.002, 3,
                                 maxwarn=1, gmx_format='g96',
                                 write_vel=True,
                                 write_force=False)
            rundir = os.path.join(tempdir, 'gmxintegrate2')
            # Create the directory for running:
            make_dirs(rundir)
            eng.exe_dir = rundir
            # Create the system:
            system = make_test_system((eng.input_files['conf'], 0))
            # Propagate:
            order_function = Position(0, dim='x', periodic=False)
            i = 0
            initial_x = -0.422
            ensemble = {'system': system, 'order_function': order_function}
            for step in eng.integrate(ensemble, 5):
                j = max(i - 1, 0)
                assert (step['order'][0] ==
                        pytest.approx(j * eng.subcycles + initial_x, abs=1e-3))
                i += 1
            eng.clean_up()
            ensemble = {'system': system}
            correct = [11.578, 14.478]
            for i_s, step in enumerate(eng.integrate(ensemble, 1)):
                j = max(i - 1, 0)
                assert (correct[i_s] ==
                        pytest.approx(j * eng.subcycles + initial_x, abs=1e-3))
                i += 1
            eng.clean_up()


class TestGromacsRunner2:
    """Test the Runner."""
    def test_init(self):
        """Test that we can initiate the runner."""
        eng = GromacsRunner(cmd='echo',
                            trr_file='No thank you',
                            edr_file='Please dont insist',
                            exe_dir='Get lost already')
        with pytest.raises(RuntimeError) as err:
            eng.check_poll()
        assert 'GROMACS is not running.' in str(err.value)
        del eng

    def test_read_remaining_trr(self):
        """Test read_remaining_trr work with corrupted trr."""
        filepath = os.path.dirname(os.path.abspath(__file__))
        ff_trr = os.path.join(filepath, '../inout/traj.trr')

        with tempfile.TemporaryDirectory() as tempdir:
            new_file = os.path.join(tempdir, 'life.trr')
            # Empty file first
            with open(new_file, 'w', encoding="utf8"):
                pass
            with open(new_file, 'rb') as now:
                read_remaining_trr(new_file, now, 0)
            shutil.copy(ff_trr, new_file)

            # Inconsistent file after
            with open(new_file, 'ab') as now:
                now.write(b'\x00\x00\x00\x00\x00\x00\x00\x00')
            with open(new_file, 'rb') as now:
                with pytest.raises(struct.error) as err:
                    for _ in read_remaining_trr(new_file, now, 0):
                        read_remaining_trr(new_file, now, 0)
            assert 'unpack requires a buffer of 8 bytes' in str(err.value)

    def test_reopen_file(self):
        """Reopen a file if the inode has changed test."""
        filepath = os.path.dirname(os.path.abspath(__file__))
        ff_temp = os.path.join(filepath, 'config.xyz')
        f_temp = tempfile.TemporaryFile(mode='w+')
        inode = os.stat(ff_temp).st_ino
        out = reopen_file(filename=ff_temp, fileh=None,
                          inode=inode, bytes_read=None)

        assert out == (None, None)

        out = reopen_file(filename=ff_temp, fileh=f_temp,
                          inode=16, bytes_read=0)
        out[0].close()
