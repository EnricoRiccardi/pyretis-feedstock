# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the CP2KEngine class."""
from io import StringIO
import logging
import os
import pytest
import tempfile
from unittest.mock import patch
import numpy as np
from pyretis.core.common import counter
from pyretis.core.system import System
from pyretis.core.path import Path
from pyretis.core.random_gen import MockRandomGenerator
from pyretis.core.particles import ParticlesExt
from pyretis.core.pathensemble import PathEnsembleExt
from pyretis.core.particlefunctions import kinetic_energy
from pyretis.engines import CP2KEngine
from pyretis.inout.common import make_dirs
from pyretis.inout.formats.xyz import (
    read_xyz_file,
    convert_snapshot,
    write_xyz_file,
)
from pyretis.orderparameter.orderparameter import PositionVelocity
from .test_helpers.test_helpers import remove_dir
logging.disable(logging.DEBUG)

HERE = os.path.abspath(os.path.dirname(__file__))


def make_test_system(conf):
    """Just make a test system with particles."""
    system = System()
    system.particles = ParticlesExt(dim=3)
    system.particles.config = conf
    return system


class TestCP2KEngine:
    """Run the tests for the CP2KEngine."""

    def test_init(self):
        """Test that we can initiate the engine."""
        dir_name = os.path.join(HERE, 'cp2k_input')
        extra_files = ['extra_file']
        engine = CP2KEngine(cp2k='cp2k',
                            input_path=dir_name,
                            timestep=0.002,
                            subcycles=10,
                            extra_files=extra_files)
        assert len(engine.extra_files) == 1
        # Check that we get an error if we are missing files:
        with pytest.raises(ValueError):
            CP2KEngine(cp2k='cp2k',
                       input_path='',
                       timestep=0.002,
                       subcycles=10,
                       extra_files=extra_files)
        # Check that non-existing extra files are not added:
        extra_files.append('a-non-existing-file-please')
        engine = CP2KEngine(cp2k='cp2k',
                            input_path=dir_name,
                            timestep=0.002,
                            subcycles=10,
                            extra_files=extra_files)
        assert len(engine.extra_files) == 1

    def test_single_step(self):
        """Test that the single step method work as we intend."""
        cmd = os.path.join(HERE, 'mockcp2k.py')
        dir_name = os.path.join(HERE, 'cp2k_input')
        extra_files = ['extra_file']
        engine = CP2KEngine(cp2k=cmd,
                            input_path=dir_name,
                            timestep=0.002,
                            subcycles=10,
                            extra_files=extra_files)
        rundir = os.path.join(HERE, 'generate')
        # Create the directory for running:
        make_dirs(rundir)
        engine.exe_dir = rundir
        # Create the system:
        system = make_test_system((engine.input_files['conf'], 0))
        # Run a single step:
        out = engine.step(system, 'cp2k_step')
        # Check that we have the expected files after the step:
        for i in ('extra_file', 'step.inp', 'conf.xyz', 'cp2k_step.xyz'):
            assert os.path.isfile(os.path.join(rundir, i))
        assert system.particles.ekin == pytest.approx(0.9)
        assert system.particles.vpot == pytest.approx(-0.9)
        # Get snapshot:
        box, xyz, vel, names = convert_snapshot(next(read_xyz_file(out)))
        for i, j in zip(names, ['H', 'H']):
            assert i == j
        assert np.allclose(box, np.array([1., 2., 3.]))
        assert np.allclose(xyz,
                           np.array([[0.9, 1.8, 2.7], [1.9, 2.8, 3.7]]))
        assert np.allclose(vel,
                           np.array([[9.9, 10.8, 11.7], [10.9, 11.8, 12.7]]))
        engine.clean_up()
        remove_dir(rundir)

    def test_modify_velocities(self):
        """Test the modify velocities method."""
        cmd = os.path.join(HERE, 'mockcp2k.py')
        dir_name = os.path.join(HERE, 'cp2k_input')
        extra_files = ['extra_file']
        engine = CP2KEngine(cp2k=cmd,
                            input_path=dir_name,
                            timestep=0.002,
                            subcycles=10,
                            extra_files=extra_files)
        with tempfile.TemporaryDirectory() as tempdir:
            # Create the directory for running:
            make_dirs(tempdir)
            engine.exe_dir = tempdir
            # Create the system:
            system = make_test_system((engine.input_files['conf'], 0))
            mass = np.array([[1.], [1.]])
            system.particles.mass = mass
            system.particles.imass = np.reciprocal(mass*1.)
            system.temperature['beta'] = 1.
            pos = engine.input_files['conf']
            _, _, vel, _ = engine._read_configuration(pos)
            kin_old = kinetic_energy(vel=vel, mass=mass)[0]

            # Modify velocities:
            ensemble = {'system': system, 'rgen': MockRandomGenerator(seed=0)}
            vel_settings = {'aimless': True, 'momentum': False}
            dek, kin_new = engine.modify_velocities(ensemble, vel_settings)
            kin_new2 = kinetic_energy(vel=system.particles.vel, mass=mass)[0]
            assert kin_new == pytest.approx(kin_new2)
            assert kin_new == pytest.approx(1.7070285604149718)
            assert dek == pytest.approx(kin_new - kin_old)
            assert dek == pytest.approx(-43.79297143958503)

            # Test aim:
            vel_settings['aimless'] = False
            vel_settings['sigma_v'] = 1
            system.particles.config = (engine.input_files['conf'], 0)
            ensemble['rgen'] = MockRandomGenerator(seed=0)
            dek, kin_new = engine.modify_velocities(ensemble, vel_settings)
            kin_new2 = kinetic_energy(vel=system.particles.vel, mass=mass)[0]
            assert kin_new == pytest.approx(kin_new2)
            assert kin_new == pytest.approx(62.89068569041497)
            assert dek == pytest.approx(kin_new - kin_old)
            assert dek == pytest.approx(17.39068569041497)

            # Test rescaling:
            system.particles.vpot = 2.
            vel_settings['rescale'] = 10
            system.particles.config = (engine.input_files['conf'], 0)
            dek, kin_new = engine.modify_velocities(ensemble, vel_settings)
            assert kin_new == pytest.approx(8.)
            assert dek == pytest.approx(0.)

            # Test negative rescaling:
            vel_settings['rescale'] = -2
            system.particles.config = (engine.input_files['conf'], 0)
            dek, kin_new = engine.modify_velocities(ensemble, vel_settings)
            # Velocity unchanged, so kin_new == kin_old.
            assert kin_new == pytest.approx(kin_old)
            assert dek == pytest.approx(0.)

            # Test reset_momentum():
            vel_settings['momentum'] = True
            vel_settings['aimless'] = True
            vel_settings['rescale'] = False
            ensemble['rgen'] = MockRandomGenerator(seed=0)
            dek, kin_new = engine.modify_velocities(ensemble, vel_settings)
            kin_new2 = kinetic_energy(vel=system.particles.vel, mass=mass)[0]
            assert kin_new == kin_new2
            assert dek == kin_new - kin_old
            assert kin_new == pytest.approx(0.08269801069751077)
            assert dek == pytest.approx(-45.41730198930249)

            # Test when kin_old = 0.0
            vel = np.array([[0.]*3, [0.]*3])
            zero_vel_path = os.path.join(tempdir, 'zero_vel.xyz')
            write_xyz_file(zero_vel_path, pos=vel, vel=vel)
            system.particles.config = (zero_vel_path, None)
            dek, kin_new = engine.modify_velocities(ensemble, vel_settings)
            kin_new2 = kinetic_energy(vel=system.particles.vel, mass=mass)[0]
            assert kin_new == pytest.approx(kin_new2)
            assert kin_new == pytest.approx(0.23911029289543265)
            assert dek == pytest.approx(float('inf'))

    def test_path_ensemble_names(self):
        """Test the path naming in path ensemble simulations."""
        cmd = os.path.join(HERE, 'mockcp2k.py')
        dir_name = os.path.join(HERE, 'cp2k_input')
        extra_files = ['extra_file']
        engine = CP2KEngine(cp2k=cmd,
                            input_path=dir_name,
                            timestep=0.002,
                            subcycles=10,
                            extra_files=extra_files)
        with tempfile.TemporaryDirectory() as tempdir:
            # Propagate in temp folder:
            engine.exe_dir = tempdir
            # Create the system:
            system = make_test_system((engine.input_files['conf'], 0))
            # Propagate:
            orderp = PositionVelocity(0, dim='x', periodic=False)
            path = Path(None, maxlen=4)
            interfaces = [0.2, 8.0, 9.0]
            p_ens = PathEnsembleExt(ensemble_number=666, interfaces=interfaces)
            ensemble = {'system': system, 'order_function': orderp,
                        'interfaces': interfaces, 'path_ensemble': p_ens}
            with patch('sys.stdout', new=StringIO()):
                success, _ = engine.propagate(path, ensemble, reverse=False)
                assert not success
                assert '666_' in path.phasepoints[0].particles.config[0]

    def test_propagate_forward(self):
        """Test the propagate method forward in time."""
        cmd = os.path.join(HERE, 'mockcp2k.py')
        dir_name = os.path.join(HERE, 'cp2k_input')
        extra_files = ['extra_file']
        engine = CP2KEngine(cp2k=cmd,
                            input_path=dir_name,
                            timestep=0.002,
                            subcycles=10,
                            extra_files=extra_files)
        rundir = os.path.join(HERE, 'generatef')
        # Create the directory for running:
        make_dirs(rundir)
        engine.exe_dir = rundir
        # Create the system:
        system = make_test_system((engine.input_files['conf'], 0))
        # Propagate:
        orderp = PositionVelocity(0, dim='x', periodic=False)
        path = Path(None, maxlen=4)
        ensemble = {'system': system, 'order_function': orderp,
                    'interfaces': [0.2, 8.0, 9.0]}
        with patch('sys.stdout', new=StringIO()):
            success, _ = engine.propagate(path, ensemble, reverse=False)
            assert not success
        engine.clean_up()
        remove_dir(rundir)

    def test_propagate_backward(self):
        """Test the propagate method forward in time."""
        cmd = os.path.join(HERE, 'mockcp2k.py')
        dir_name = os.path.join(HERE, 'cp2k_input')
        extra_files = ['extra_file']
        engine = CP2KEngine(cmd,
                            input_path=dir_name,
                            timestep=0.002,
                            subcycles=10,
                            extra_files=extra_files)
        rundir = os.path.join(HERE, 'generateb')
        # Create the directory for running:
        make_dirs(rundir)
        engine.exe_dir = rundir
        # Create the system:
        system = make_test_system((engine.input_files['conf'], 0))
        # Propagate:
        orderp = PositionVelocity(0, dim='x', periodic=False)
        path = Path(None, maxlen=4)
        counter.count = -1
        ensemble = {'system': system, 'order_function': orderp,
                    'interfaces': [0.2, 0.5, 0.8]}
        with patch('sys.stdout', new=StringIO()):
            success, _ = engine.propagate(path, ensemble, reverse=True)
            assert success
        # Check that initial velocities were reversed:
        infile = os.path.join(rundir, '0_conf.xyz')
        _, _, vel, _ = convert_snapshot(next(read_xyz_file(infile)))
        outfile = os.path.join(rundir, 'r_0_conf.xyz')
        _, _, rvel, _ = convert_snapshot(next(read_xyz_file(outfile)))
        assert np.allclose(vel, -1.0 * rvel)
        engine.clean_up()
        remove_dir(rundir)

    def test_integrate(self):
        """Test the integrate method."""
        cmd = os.path.join(HERE, 'mockcp2k.py')
        dir_name = os.path.join(HERE, 'cp2k_input')
        extra_files = ['extra_file']
        engine = CP2KEngine(cmd,
                            input_path=dir_name,
                            timestep=0.002,
                            subcycles=10,
                            extra_files=extra_files)
        rundir = os.path.join(HERE, 'cp2kintegrate')
        # Create the directory for running:
        make_dirs(rundir)
        engine.exe_dir = rundir
        # Create the system:
        system = make_test_system((engine.input_files['conf'], 0))
        # Propagate:
        orderp = PositionVelocity(0, dim='x', periodic=False)
        correct_order = [
            [0.77317158100000005, 1.0],
            [0.9, 9.9],
            [0.9, 9.9],
        ]
        ensemble = {'system': system, 'order_function': orderp}
        for i, step in enumerate(engine.integrate(ensemble, 2)):
            assert np.allclose(step['order'], correct_order[i])
        # Call once more without specifying order_function
        ensemble = {'system': system}
        for step in engine.integrate(ensemble, 3):
            assert "order" not in step
        # Call once more, after first creating a fake backup-file:
        bfile = os.path.join(rundir, 'backup.bak-1')
        with open(bfile, 'w', encoding="utf8") as output:
            output.write('Just a line')
        listoffiles = [entry.name for entry in os.scandir(rundir)]
        assert "backup.bak-1" in listoffiles
        for step in engine.integrate(ensemble, 1):
            assert "order" not in step
        # Check that PyRETIS has removed the backup-file.
        listoffiles = [entry.name for entry in os.scandir(rundir)]
        assert "backup.bak-1" not in listoffiles
        engine.clean_up()
        remove_dir(rundir)

    def test_dump_config(self):
        """Test the dump_config method."""
        cmd = os.path.join(HERE, 'mockcp2k.py')
        dir_name = os.path.join(HERE, 'cp2k_input')
        extra_files = ['extra_file']
        engine = CP2KEngine(cmd,
                            input_path=dir_name,
                            timestep=0.002,
                            subcycles=10,
                            extra_files=extra_files)
        trajdir = os.path.join(HERE, 'cp2k_fake_trajectory')
        # Create the directory for storing a fake trajectory file:
        make_dirs(trajdir)
        engine.exe_dir = trajdir
        # Set up for a trajectory file:
        trajfile = os.path.join(trajdir, 'traj.xyz')
        # Set up for a config file:
        conffile = os.path.join(trajdir, 'conf.xyz')
        lines = [
            "2",
            "# Step: 0 Box:    1.0000    1.1    1.2",
            "O         0.71    -0.361     0.332  1.1 2.2 3.6",
            "H         0.72    -0.362     0.333  1.2 2.3 3.7",
            "2",
            "# Step: 1 Box:    1.0000    1.1    1.2",
            "O         0.61    -0.331     0.232  1.5 2.2 3.7",
            "H         0.62    -0.366     0.133  2.2 3.3 4.7",
        ]
        with open(trajfile, 'w', encoding="utf8") as output:
            output.write('\n'.join(lines))
        for i in range(3):
            config = (trajfile, i)
            if i == 0:  # Config file should not exist.
                assert not os.path.isfile(conffile)
            if i == 1:  # Config file should exist, and will be overwritten.
                assert os.path.isfile(conffile)
            engine.dump_config(config)
            if i == 1:
                #  Note: is it possible to check the log-file here?
                #  logger.warning('CP2K will overwrite %s', out_file)
                continue
            if i == 2:
                # Note: is it possible to check the log-file here?
                # logger.error('CP2K could not extract index %i ....
                break
            lines2 = None
            with open(conffile, 'r', encoding="utf8") as infile:
                lines2 = infile.readlines()
            len2 = len(lines2)
            for j in range(len2):
                words1 = lines[i*len2+j].split()
                words2 = lines2[j].split()
                if j == 0:
                    assert 1 == len(words2)
                    assert words1 == words2
                if j == 1:
                    assert words1 != words2
                    assert words2[1] == "Box:"
                if j > 1:
                    # First item is the particle type/elements:
                    assert words1[0] == words2[0]
                    # The rest is just numbers:
                    for k in range(1, len(words2)):
                        assert (float(words1[k]) ==
                                pytest.approx(float(words2[k])))
        engine.clean_up()
        remove_dir(trajdir)
