# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the GROMACS Engine class."""
import logging
import numpy as np
import os
import shutil
import tempfile
import pytest
from pyretis.core.path import Path
from pyretis.core.common import counter
from pyretis.engines import GromacsEngine
from pyretis.inout.common import make_dirs
from pyretis.inout.formats.gromacs import (read_gromacs_gro_file,
                                           read_gromacs_generic)
from pyretis.orderparameter.orderparameter import Position
from .test_helpers.test_helpers import make_test_system

logging.disable(logging.CRITICAL)
HERE = os.path.abspath(os.path.dirname(__file__))
GMX = os.path.join(HERE, 'mockgmx.py')
MDRUN = os.path.join(HERE, 'mockmdrun.py')
GMX_DIR = os.path.join(HERE, 'gmx_input')
GMX_DIR2 = os.path.join(HERE, 'gmx_input2')
GMX_DIR3 = os.path.join(HERE, 'gmx_input3')
GMX_LOAD = os.path.join(HERE, '../initiation/gromacs/gmx_input')


def copy_gmx_input(tempdir, source=GMX_DIR, name='gmx_input'):
    """Copy GROMACS input files to a private test directory."""
    inputdir = os.path.join(tempdir, name)
    shutil.copytree(source, inputdir)
    return inputdir


class TestGromacsEngine:
    """Run the tests for the GROMACS Engine."""

    def test_init(self):
        """Test that we can initiate the engine."""
        with tempfile.TemporaryDirectory() as tempdir:
            inputdir = copy_gmx_input(tempdir)
            eng = GromacsEngine(gmx='echo',
                                mdrun='echo',
                                input_path=inputdir,
                                timestep=0.002,
                                subcycles=10,
                                maxwarn=10,
                                gmx_format='gro',
                                write_vel=True,
                                write_force=False)
            eng.exe_dir = inputdir
            inputdir3 = copy_gmx_input(tempdir, source=GMX_DIR3,
                                       name='gmx_input3')
            eng2 = GromacsEngine(gmx='echo',
                                 mdrun='echo',
                                 input_path=inputdir3,
                                 timestep=0.002,
                                 subcycles=10,
                                 maxwarn=10,
                                 gmx_format='gro',
                                 write_vel=True,
                                 write_force=False)
            eng2.exe_dir = inputdir3
        with pytest.raises(ValueError):
            GromacsEngine(gmx='echo',
                          mdrun='echo',
                          input_path='gmx_input',
                          timestep=0.002,
                          subcycles=10,
                          gmx_format='not-a-format')

    def test_read_generic(self):
        """Test a single read."""
        for i, _ in enumerate(read_gromacs_generic('no_file.g96')):
            max_i = i
        assert max_i == 0
        any_g96 = os.path.join(HERE, '..', 'inout', 'config.g96')
        for i, _ in enumerate(read_gromacs_generic(any_g96)):
            max_i = i
        assert max_i == 0
        any_gro = os.path.join(HERE, 'gmx_input', 'conf.gro')
        for i, _ in enumerate(read_gromacs_generic(any_gro)):
            max_i = i
        assert max_i == 0
        any_trr = os.path.join(HERE, '..', 'tools', '2water.trr')
        for i, _ in enumerate(read_gromacs_generic(any_trr)):
            max_i = i
        assert max_i == 5

    def test_single_step(self):
        """Test a single step using the MOCK GROMACS engine."""
        with tempfile.TemporaryDirectory() as tempdir:
            inputdir = copy_gmx_input(tempdir)
            eng = GromacsEngine(gmx=GMX,
                                mdrun=MDRUN,
                                input_path=inputdir,
                                timestep=0.002,
                                subcycles=3,
                                maxwarn=1,
                                gmx_format='gro',
                                write_vel=True,
                                write_force=False)
            rundir = os.path.join(tempdir, 'generategmx')
            # Create the directory for running:
            make_dirs(rundir)
            eng.exe_dir = rundir
            # Create the system:
            system = make_test_system((eng.input_files['conf'], 0))
            out = eng.step(system, 'gmx_mock_step')
            assert out == 'gmx_mock_step.gro'
            # Check that output files are present:
            for i in ('conf.gro', 'gmx_mock_step.gro'):
                assert os.path.isfile(os.path.join(rundir, i))
            # Check that output files contain the expected data:
            _, xyz1, vel1, _ = read_gromacs_gro_file(
                os.path.join(rundir, 'conf.gro')
            )
            _, xyz2, vel2, _ = read_gromacs_gro_file(
                os.path.join(rundir, 'gmx_mock_step.gro')
            )
            assert np.allclose(xyz2 - xyz1,
                               eng.subcycles * np.ones_like(xyz1))
            assert np.allclose(
                vel1,
                np.repeat([0.1111, 0.2222, 0.3333], 27).reshape(3, 27).T
            )
            assert np.allclose(vel2, np.ones_like(vel2))
            # Check the final state:
            assert eng.subcycles * 1.0 == pytest.approx(system.particles.ekin)
            assert eng.subcycles * -1.0 == \
                pytest.approx(system.particles.vpot)
            assert system.particles.get_pos()[0] == \
                os.path.join(rundir, 'gmx_mock_step.gro')
            eng.clean_up()

    def test_modify_velocities(self):
        """Test the modify velocities method."""
        with tempfile.TemporaryDirectory() as tempdir:
            inputdir = copy_gmx_input(tempdir)
            eng = GromacsEngine(gmx=GMX,
                                mdrun=MDRUN,
                                input_path=inputdir,
                                timestep=0.002,
                                subcycles=10,
                                maxwarn=1,
                                gmx_format='gro',
                                write_vel=False,
                                write_force=False)
            rundir = os.path.join(tempdir, 'generategmxvel')
            # Create the directory for running:
            make_dirs(rundir)
            eng.exe_dir = rundir
            # Create the system:
            system = make_test_system((eng.input_files['conf'], 0))
            ensemble = {'system': system}
            dek, kin_new = eng.modify_velocities(ensemble,
                                                 {'aimless': True,
                                                  'momentum': False})

            assert kin_new == pytest.approx(1234.5678)
            assert dek == float('inf')
            # Check that aiming fails:
            with pytest.raises(NotImplementedError):
                eng.modify_velocities(ensemble, {'aimless': False,
                                                 'momentum': False})
            # Check that rescaling fails:
            with pytest.raises(NotImplementedError):
                eng.modify_velocities(ensemble, {'aimless': False,
                                                 'momentum': True,
                                                 'rescale': 11})

            dek, kin_new = eng.modify_velocities(ensemble,
                                                 {'aimless': True,
                                                  'momentum': False})
            assert kin_new == pytest.approx(1234.5678)
            assert dek == pytest.approx(0.0)
            eng.clean_up()

    def test_propagate_forward(self):
        """Test the propagate method forward in time."""
        with tempfile.TemporaryDirectory() as tempdir:
            inputdir = copy_gmx_input(tempdir)
            eng = GromacsEngine(gmx=GMX,
                                mdrun=MDRUN,
                                input_path=inputdir,
                                timestep=0.002,
                                subcycles=7,
                                maxwarn=1,
                                gmx_format='gro',
                                write_vel=True,
                                write_force=False)
            rundir = os.path.join(tempdir, 'generategmxf')
            # Create the directory for running:
            make_dirs(rundir)
            eng.exe_dir = rundir
            # Create the system:
            system = make_test_system((eng.input_files['conf'], 0))
            # Propagate:
            orderp = Position(0, dim='x', periodic=False)
            path = Path(None, maxlen=8)
            ensemble = {'system': system,
                        'order_function': orderp,
                        'interfaces': [-0.45, 10.0, 14.0]}
            success, _ = eng.propagate(path, ensemble, reverse=False)
            initial_x = -0.422
            for i, point in enumerate(path.phasepoints):
                assert point.particles.ekin == \
                    pytest.approx(eng.subcycles * i)
                assert point.particles.vpot == \
                    pytest.approx(-1.0 * eng.subcycles * i)
                assert point.order[0] == pytest.approx(i * eng.subcycles +
                                                       initial_x, abs=1e-3)
            assert success
            eng.clean_up()

    def test_propagate_backward(self):
        """Test the propagate method backward in time."""
        with tempfile.TemporaryDirectory() as tempdir:
            inputdir = copy_gmx_input(tempdir)
            eng = GromacsEngine(gmx=GMX,
                                mdrun=MDRUN,
                                input_path=inputdir,
                                timestep=0.002,
                                subcycles=3,
                                maxwarn=1,
                                gmx_format='gro',
                                write_vel=True,
                                write_force=False)
            rundir = os.path.join(tempdir, 'generategmxb')
            # Create the directory for running:
            make_dirs(rundir)
            eng.exe_dir = rundir
            # Create the system:
            system = make_test_system((eng.input_files['conf'], 0))
            # Propagate:
            orderp = Position(0, dim='x', periodic=False)
            path = Path(None, maxlen=3)
            ensemble = {'system': system,
                        'order_function': orderp,
                        'interfaces': [-0.45, 10.0, 14.0]}
            counter.count = -1
            success, _ = eng.propagate(path, ensemble, reverse=True)
            assert not success
            # Check that velocities were reversed:
            _, _, vel1, _ = read_gromacs_gro_file(
                os.path.join(rundir, '0_conf.gro')
            )
            _, _, vel2, _ = read_gromacs_gro_file(
                os.path.join(rundir, 'r_0_conf.gro')
            )
            assert np.allclose(vel1, -1.0 * vel2)
            eng.clean_up()

    def test_integrate(self):
        """Test the integrate method."""
        with tempfile.TemporaryDirectory() as tempdir:
            inputdir = copy_gmx_input(tempdir)
            eng = GromacsEngine(gmx=GMX,
                                mdrun=MDRUN,
                                input_path=inputdir,
                                timestep=0.002,
                                subcycles=3,
                                maxwarn=1,
                                gmx_format='gro',
                                write_vel=True,
                                write_force=False)
            rundir = os.path.join(tempdir, 'gmxintegrate')
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
            for step in eng.integrate(ensemble, 3):
                assert step['order'][0] == pytest.approx(i * eng.subcycles +
                                                         initial_x, abs=1e-3)
                i += 1
            # Clean
            eng.clean_up()

    def test_select_energy_term(self):
        """Test the select_energy_term method."""
        with tempfile.TemporaryDirectory() as tempdir:
            inputdir = copy_gmx_input(tempdir)
            eng = GromacsEngine(gmx=GMX,
                                mdrun=MDRUN,
                                input_path=inputdir,
                                timestep=0.002,
                                subcycles=3,
                                gmx_format='gro')

            terms = 'full'
            expected = (
                r"b'Potential\nKinetic-En.\nTotal-Energy\nTemperature"
                r"\nPressure'"
            )
            assert str(eng.select_energy_terms(terms)) == expected

            terms = 'path'
            assert str(eng.select_energy_terms(terms)) == \
                r"b'Potential\nKinetic-En.'"

            terms = 'non_allowed_term'  # Will be equal to terms = 'path'.
            assert str(eng.select_energy_terms(terms)) == \
                r"b'Potential\nKinetic-En.'"

    def test_check_fails(self):
        """Test behavior if orderfunction not given."""
        with tempfile.TemporaryDirectory() as tempdir:
            inputdir = copy_gmx_input(tempdir)
            eng = GromacsEngine(gmx=GMX,
                                mdrun=MDRUN,
                                input_path=inputdir,
                                timestep=0.002,
                                subcycles=1,
                                gmx_format='gro')
            rundir = os.path.join(tempdir, 'gmxintegrate')
            # Create the directory for running:
            make_dirs(rundir)
            eng.exe_dir = rundir
            # Create the system:
            system = make_test_system((eng.input_files['conf'], 0))
            # Propagate:
            for step in eng.integrate({'system': system}, 1):
                assert "thermo" in step
                assert "order" not in step
            order_function = Position(0, dim='x', periodic=False)
            eng.ext = "notexistingextension"
            ensemble = {'system': system, 'order_function': order_function}
            # Propagate:
            with pytest.raises(Exception) as context:
                for _ in eng.integrate(ensemble, 1):
                    pass
            assert "GROMACS engine does not support reading" in \
                str(context.value)
            path = Path(None, maxlen=3)
            ensemble = {'system': system, 'order_function': order_function,
                        'interfaces':  [-1., 0.0, 1.0]}
            with pytest.raises(Exception) as context:
                eng.propagate(path, ensemble, reverse=True)
            assert "GROMACS engine does not support writing" in \
                str(context.value)

            # Add a non existing .ndx file and get an error
            eng.input_files['index'] = 'ufo.ndx'
            with pytest.raises(Exception) as context:
                eng._execute_grompp('should_be_a_grompp.mdp',
                                    'configuration.gro', 'a_label')
            assert "(GROMACS engine) failed" in str(context.value)

            eng.clean_up()

    def test_trjconv(self):
        """Test behavior when extracting frames."""
        # We only test the trajectory conversion here as the
        # other parts are tested elsewhere.
        with tempfile.TemporaryDirectory() as tempdir:
            inputdir = copy_gmx_input(tempdir)
            eng = GromacsEngine(gmx=GMX,
                                mdrun=MDRUN,
                                input_path=inputdir,
                                timestep=0.002,
                                subcycles=1,
                                gmx_format='gro')
            eng.exe_dir = tempdir
            # Use an already done TRR file:
            any_trr = os.path.join(HERE, '..', 'tools', '2water.trr')
            any_gro = os.path.join(HERE, '..', 'tools', '2water.gro')
            eng.input_files['conf'] = any_gro
            out_gro = os.path.join(tempdir, 'output.gro')
            eng._extract_frame(any_trr, 4, out_gro)
            # Check the created output file:
            assert os.path.isfile(out_gro)
            with open(out_gro, 'r', encoding="utf8") as infile:
                data = [i.strip().split() for i in infile]
            # Note: this is a fake output from a fake engine
            assert data[0][1] == 'froggy'

            # This tests only the creation of the conf.g96
            shutil.copytree(GMX_LOAD, tempdir+'2')
            anytrr = os.path.join(GMX_LOAD,
                                  '../loader/023/accepted/trajB.trr')
            eng = GromacsEngine(gmx=GMX,
                                mdrun=MDRUN,
                                input_path=tempdir+'2',
                                timestep=666,
                                subcycles=1,
                                gmx_format='g96')
            out_gro = os.path.join(tempdir+'2', 'my_uncle_jonny.g96')
            # Note, we use a fake engine and it will just copy the file
            eng.exe_dir = tempdir+'2'
            eng._extract_frame(anytrr, 0, out_gro)
            assert os.path.isfile(out_gro)

            shutil.rmtree(os.path.join(tempdir+'2'), 'conf.gro')
            # A partial test to avoid the need to install GROMACS to run
            with pytest.raises(RuntimeError) as err:
                GromacsEngine(gmx=GMX,
                              mdrun=MDRUN,
                              input_path=tempdir+'2',
                              timestep=0.002,
                              subcycles=1,
                              gmx_format='g96')
            assert 'GROMACS' in str(err.value)
