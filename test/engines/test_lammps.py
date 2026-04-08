# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the LAMMPSEngine class."""
import logging
import os
import shutil
import tempfile
import pytest
import numpy as np
from pyretis.core.random_gen import create_random_generator
from pyretis.core.system import System
from pyretis.core.particles import ParticlesExt
from pyretis.core.path import Path
from pyretis.engines.lammps import (
    create_lammps_md_input,
    read_lammps_input,
    read_lammps_log,
    system_to_lammps,
    LAMMPSEngine,
)


logging.disable(logging.CRITICAL)
HERE = os.path.abspath(os.path.dirname(__file__))
INPUT_PATH = os.path.join(HERE, 'lammps_input')
# Define some error messages we expect:
MISSING_FILE_MSG = 'Missing input file'
MISSING_REVVEL_MSG = "'reverse_velocities'"
AIMLESS_ERROR_MSG = 'LAMMPS only support the aimless velocity modification.'


LAMMPS_IN = [
    ('units', 'real'),
    ('timestep', '1.0'),
    ('atom_style', 'full'),
    ('bond_style', 'hybrid morse harmonic'),
    ('angle_style', 'hybrid class2 harmonic'),
    ('dihedral_style', 'class2'),
    ('improper_style', 'class2'),
    ('read_data', 'system.data'),
    ('include', '"groups.in"'),
    ('include', '"md.in.settings"'),
    ('thermo', '100'),
    ('variable', 'SET_TEMP equal 343'),
    ('include', '"order.in"'),
    ('fix', '1 mg nvt temp 343 343 100'),
    ('fix', '2 dummy nvt temp 343 343 10'),
    ('fix', '3 water rigid/nvt molecule temp 343 343 100'),
]


LAMMPS_ENERGY_DATA = {
    'Step': [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
    'Temp': [354.83622, 354.87936, 353.55139, 346.73993, 355.86763,
             356.02159, 351.24368, 338.50584, 350.65277, 341.57213,
             337.38405],
    'Press': [1503.3635, -1496.9831, -2155.9235, 2692.286,
              3461.4022, 522.80391, -190.30433, 55.010487,
              138.91587, 1222.3388, 343.12519],
    'PotEng': [-5047.0936, -4349.3609, -4642.5397, -1634.8362,
               -1475.5791, -3255.0873, -3918.0104, -3858.1462,
               -3498.4035, -3128.2873, -3484.1554],
    'KinEng': [1080.9687, 1081.1002, 1077.0547, 1056.3043,
               1084.1108, 1084.5799, 1070.0245, 1031.2201,
               1068.2243, 1040.5612, 1027.8027],
    'TotEng': [-3966.1248, -3268.2607, -3565.4851, -578.53193,
               -391.46827, -2170.5074, -2847.9859, -2826.9261,
               -2430.1791, -2087.7261, -2456.3527],
}


def get_setting(keyword, settings):
    """Return settings for a given keyword.

    Parameters
    ----------
    keyword : string
        The keyword to find settings for.
    settings : list of tuples
        The settings to search for the keyword in. It is assumed
        that the first item in the tuple is the keyword and the
        second the setting.

    Returns
    -------
    out : list of strings
        The settings found for a given keyword.

    """
    found = []
    for (key, val) in settings:
        if key.lower() == keyword.lower():
            found.append(val)
    return found


def create_system_for_lammps():
    """Create a test system for LAMMPS."""
    system = System()
    system.particles = ParticlesExt(dim=3)
    system.particles.set_pos(('system.data', 0))
    system.particles.set_vel(False)
    return system


def read_seed(lammps_input):
    """Read the seed from the given LAMMPS input file."""
    settings = read_lammps_input(lammps_input)
    seed = None
    for var in settings:
        if var[1].startswith('VEL_SEED'):
            seed = int(var[1].split('equal')[-1])
    return seed


def read_raw_file(filename):
    """Return the lines from a file."""
    lines = []
    with open(filename, 'r', encoding="utf8") as infile:
        lines = [i.strip() for i in infile]
    return lines


class TestLAMMPSEngineMethods:
    """Run the tests for the methods defined in LAMMPS module."""

    def test_read_lammps_input(self):
        """Test that we can read LAMMPS input."""
        infile = os.path.join(HERE, 'lammps_input', 'lammps.in')
        settings = read_lammps_input(infile)
        assert len(LAMMPS_IN) == len(settings)
        for i, j in zip(LAMMPS_IN, settings):
            assert i == j
        timestep = get_setting('timestep', settings)
        assert timestep == ['1.0']
        fix = get_setting('fix', settings)
        correct_fix = [item[1] for item in settings[-3:]]
        assert len(correct_fix) == len(fix)
        for i, j in zip(correct_fix, fix):
            assert i == j

    def test_create_md_input(self):
        """Test that we can create input files for LAMMPS."""
        infile = os.path.join(HERE, 'lammps_input', 'lammps.in')
        system = create_system_for_lammps()
        settings = {
            'subcycles': 100,
            'traj': 'test.lammpstrj',
            'steps_subcycles': 12345,
            'timestep': 0.006789,
            'generate_vel': {'seed': 9876, 'rescale': 1},
            'interfaces': [-1, 0, 1],
            'dimension': 3,
            'reverse_velocities': True,
        }
        with tempfile.NamedTemporaryFile() as temp:
            written = create_lammps_md_input(
                system, infile, temp.name, settings
            )
            temp.flush()
            data = []
            with open(temp.name, 'r', encoding="utf8") as infile:
                data = [i.strip() for i in infile]
            # Here written and data will differ, since data contains
            # the full output. But written should be contained in data,
            # check that the last lines are equal:
            offset = len(written)
            for i, line in enumerate(written):
                assert line == data[-(offset - i)]

    def test_read_lammps_log(self):
        """Test that we can read some data from a LAMMPS log file."""
        infile = os.path.join(HERE, 'lammps_input', 'pyretis_md.log')
        logdata = read_lammps_log(infile)
        assert 'energy' in logdata
        assert len(logdata['energy']) == len(LAMMPS_ENERGY_DATA)
        for i in LAMMPS_ENERGY_DATA:
            assert i in logdata['energy']
        for key, val in LAMMPS_ENERGY_DATA.items():
            assert np.allclose(val, logdata['energy'][key])

    def test_system_to_lammps(self):
        """Test the conversion of systems to LAMMPS commands."""
        system = System()
        system.particles = ParticlesExt(dim=3)

        # 1) Test for a restart file:
        system.particles.set_pos(('some.restart', 0))
        system.particles.set_vel(False)
        correct = [
            '# PyRETIS requested a LAMMPS restart file:',
            'read_restart some.restart',
        ]
        cmd = system_to_lammps(system, False, 3)
        assert correct == cmd

        # 2) Test for a restart file & reversing velocities:
        correct = [
            '# PyRETIS requested a LAMMPS restart file:',
            'read_restart some.restart',
            '# PyRETIS requested reversing of velocities:',
            'variable        vx atom -vx',
            'variable        vy atom -vy',
            'variable        vz atom -vz',
            'velocity        all set v_vx v_vy v_vz',
            'run 0',
        ]
        system.particles.set_pos(('some.restart', 0))
        system.particles.set_vel(True)
        cmd = system_to_lammps(system, True, 3)
        assert correct == cmd

        # 3) Test if we give a file with an index:
        correct = [
            '# PyRETIS requested the following snapshot:',
            'read_dump notsome.restart 10 x y z vx vy vz ix iy iz box yes',
            '# Reset time step:',
            'reset_timestep 0',
        ]
        system.particles.set_pos(('notsome.restart', 10))
        system.particles.set_vel(False)
        cmd = system_to_lammps(system, False, 3)
        assert correct == cmd

        # 4) Test index + reversing velocities:
        system.particles.set_pos(('dump.lammpstrj', 10))
        correct = [
            '# PyRETIS requested the following snapshot:',
            'read_dump dump.lammpstrj 10 x y z vx vy vz ix iy iz box yes',
            '# Reset time step:',
            'reset_timestep 0',
            '# PyRETIS requested reversing of velocities:',
            'variable        vx atom -vx',
            'variable        vy atom -vy',
            'variable        vz atom -vz',
            'velocity        all set v_vx v_vy v_vz',
            'run 0',
        ]
        system.particles.set_vel(True)
        cmd = system_to_lammps(system, True, 3)
        assert correct == cmd


class TestLAMMPSEngine:
    """Run the tests for the LAMMPSEngine class."""

    def test_lammps_initiate(self):
        """Test that we can initiate the LAMMPS engine."""
        # This should not fail:
        engine = LAMMPSEngine('lmp_serial', INPUT_PATH, 10, extra_files=None)
        assert engine.lmp == 'lmp_serial'
        assert engine.engine_type == 'external'
        assert engine.description == 'LAMMPS Engine'
        # This should also fail - input files not found:
        extra_files = ['Oh Lawd He Comin']
        with pytest.raises(ValueError) as context:
            LAMMPSEngine('lmp_serial', INPUT_PATH, 10,
                         extra_files=extra_files)
        except_str = str(context.value)
        assert except_str.startswith(MISSING_FILE_MSG)
        assert extra_files[0] in except_str

    def test_run_lammps(self):
        """Test the method to run LAMMPS."""
        # Note: More of the functionality of this method
        # is implicitly tested in the methods that test
        # integration and propagation.
        system = create_system_for_lammps()
        with tempfile.TemporaryDirectory() as tempdir:
            shutil.copy(os.path.join(HERE, 'aprogram.py'), tempdir)
            engine = LAMMPSEngine(os.path.join(HERE, 'aprogram.py'),
                                  INPUT_PATH, 10, extra_files=None)
            engine.exe_dir = tempdir
            # First we test with missing settings:
            settings = {}
            # None of the required settings were given, engine should fail:
            with pytest.raises(KeyError) as context:
                engine.run_lammps(system, settings, 'lammps_test')
            except_str = str(context.value)
            assert except_str == MISSING_REVVEL_MSG
            # Add the missing required setting(s):
            settings['steps_subcycles'] = 1
            settings['reverse_velocities'] = False
            # By construction of "aprogram.py" this just gives a
            # RuntimeError:
            with pytest.raises(RuntimeError):
                engine.run_lammps(system, settings, 'lammps_test')
            # Check the standard out and error:
            correct = [
                ['This is a program for testing external commands.'],
                ['ERROR: Program got arguments:',
                 (os.path.join(HERE, 'aprogram.py') +
                  ' -in lammps_test.in '
                  '-l lammps_test.log -screen '
                  'lammps_test.screen')]
            ]
            for filei, correcti in zip(('stdout.txt', 'stderr.txt'), correct):
                with open(os.path.join(tempdir, filei), 'r',
                          encoding="utf8") as output:
                    lines = output.readlines()
                    assert len(lines) == len(correcti)
                    for i, j in zip(lines, correcti):
                        assert i.strip() == j

    def test_add_files(self):
        """Test that the LAMMPS engine adds input files."""
        extra_files = ['order.in']
        engine = LAMMPSEngine('lammps', INPUT_PATH, 10,
                              extra_files=extra_files)
        with tempfile.TemporaryDirectory() as tempdir:
            engine.add_input_files(tempdir)
            files = [i.name for i in os.scandir(tempdir) if i.is_file()]
            expected = ['system.data'] + extra_files
            assert len(files) == len(expected)
            for i in files:
                assert i in expected

    def test_integrate(self):
        """Test the integrate method.

        Note that we do not run LAMMPS here, we just test that we
        can create files etc. for LAMMPS.

        """
        # We expect that the following files will be created:
        expected_files = [
            'pyretis_md.screen',
            'pyretis_md.log',
            'pyretis_md.lammpstrj',
            'order_pyretis_md.txt',
        ]
        steps = 10
        with tempfile.TemporaryDirectory() as tempdir:
            shutil.copy(os.path.join(HERE, 'mocklammps.py'), tempdir)
            engine = LAMMPSEngine(os.path.join(HERE, 'mocklammps.py'),
                                  INPUT_PATH, 10, extra_files=None)
            engine.exe_dir = tempdir
            engine.add_input_files(tempdir)
            system = create_system_for_lammps()
            # The expected files should not be present yet:
            for i in expected_files:
                assert not os.path.isfile(os.path.join(tempdir, i))
            engine.integrate(system, steps)
            # We now expect to find the files:
            for i in expected_files:
                assert os.path.isfile(os.path.join(tempdir, i))
            # And we expect the system to have an order parameter:
            assert (os.path.join(tempdir, 'pyretis_md.lammpstrj') ==
                    system.particles.get_pos()[0])
            assert engine.subcycles * steps == system.particles.get_pos()[1]
            assert len(system.order) == 2
            assert system.order[0] == engine.subcycles * steps + 1
            assert system.order[1] == -(engine.subcycles * steps + 1)

    def test_modify_velocities(self):
        """Test the modify_velocities method."""
        with tempfile.TemporaryDirectory() as tempdir:
            shutil.copy(os.path.join(HERE, 'mocklammps.py'), tempdir)
            engine = LAMMPSEngine(os.path.join(HERE, 'mocklammps.py'),
                                  INPUT_PATH, 10, extra_files=None)
            engine.exe_dir = tempdir
            engine.add_input_files(tempdir)
            system = create_system_for_lammps()
            # Add some fake data:
            system.order = [-1, -1]
            system.particles.ekin = None
            ensemble = {'system': system, 'rgen': None}
            vel_settings = {'aimless': True}
            dek, kin_new = engine.modify_velocities(ensemble, vel_settings)
            assert 4.0 == pytest.approx(kin_new)
            assert float('inf') == dek
            assert np.allclose([1., -1.], system.order)
            # Check the calculate order:
            assert system.order is engine.calculate_order(ensemble)
            # Check the created input file for LAMMPS:
            seed = read_seed(os.path.join(tempdir, 'generate_vel.in'))
            assert seed == 1
            # Test that we can set the seed:
            rgen = create_random_generator({'rgen': 'mock'})
            ensemble = {'system': system, 'rgen': rgen}
            vel_settings = {'aimless': True}
            dek, kin_new = engine.modify_velocities(ensemble, vel_settings)
            assert 4.0 == pytest.approx(kin_new)
            assert 0.0 == pytest.approx(dek)
            assert np.allclose([1., -1.], system.order)
            seed = read_seed(os.path.join(tempdir, 'generate_vel.in'))
            assert 1675209429 == pytest.approx(seed)
            with pytest.raises(ValueError) as context:
                vel_settings = {'aimless': False}
                engine.modify_velocities(ensemble, vel_settings)
            except_str = str(context.value)
            assert except_str == AIMLESS_ERROR_MSG

    def test_propagate_forward(self):
        """Test the propagate method."""
        with tempfile.TemporaryDirectory() as tempdir:
            shutil.copy(os.path.join(HERE, 'mocklammps.py'), tempdir)
            engine = LAMMPSEngine(os.path.join(HERE, 'mocklammps.py'),
                                  INPUT_PATH, 10, extra_files=None)
            engine.exe_dir = tempdir
            engine.add_input_files(tempdir)
            system = create_system_for_lammps()
            # Run "forward" until we cross the last interface:
            path = Path(rgen=None, maxlen=10)
            ensemble = {
                'system': system,
                'order_function': None,
                'interfaces': [-1.0, 0.0, 80.0]
            }
            success, _ = engine.propagate(
                path=path, ensemble=ensemble, reverse=False
            )
            # Check that the stopping condition is correctly written
            # to the input file.
            correct = (
                'fix op_stop_right all halt 10 v_op_1 > 80.0',
                'fix op_stop_left all halt 10 v_op_1 < -1.0'
            )
            pref = [i for i in os.listdir(tempdir) if i[0].isdigit()][0][0]
            lmp_in = read_raw_file(os.path.join(tempdir, pref + '_trajF.in'))
            for i in correct:
                assert i in lmp_in
            assert success
            assert path.length == 9
            traj = [i.particles.get_pos() for i in path.phasepoints]
            for i, point in enumerate(traj):
                assert point[1] == i*engine.subcycles
                assert tempdir in point[0]
                assert 'trajF.lammpstrj' in point[0]
            assert not any(i.particles.get_vel() for i in path.phasepoints)
            # Run "forward" until we exceed the max length:
            path = Path(rgen=None, maxlen=4)
            ensemble = {
                'system': system,
                'order_function': None,
                'interfaces': [-1., 0., float('inf')],
            }
            success, _ = engine.propagate(
                path=path, ensemble=ensemble, reverse=False,
            )
            assert not success
            assert path.length == path.maxlen
            traj = [i.particles.get_pos() for i in path.phasepoints]
            for i, point in enumerate(traj):
                assert point[1] == i*engine.subcycles
                assert tempdir in point[0]
                assert 'trajF.lammpstrj' in point[0]
            assert not any(i.particles.get_vel() for i in path.phasepoints)

    def test_propagate_backward(self):
        """Test the propagate method."""
        with tempfile.TemporaryDirectory() as tempdir:
            shutil.copy(os.path.join(HERE, 'mocklammps.py'), tempdir)
            engine = LAMMPSEngine(os.path.join(HERE, 'mocklammps.py'),
                                  INPUT_PATH, 10, extra_files=None)
            engine.exe_dir = tempdir
            engine.add_input_files(tempdir)
            system = create_system_for_lammps()
            # Run "forward" until we cross the last interface:
            path = Path(rgen=None, maxlen=4)
            ensemble = {
                'system': system,
                'order_function': None,
                'interfaces': [-1, 0.0, 42.0]
            }
            success, _ = engine.propagate(
                path=path, ensemble=ensemble, reverse=True
            )
            assert not success
            assert path.length == 4
            traj = [i.particles.get_pos() for i in path.phasepoints]
            for i, point in enumerate(traj):
                assert point[1] == i*engine.subcycles
                assert 'trajB.lammpstrj' in point[0]
                assert tempdir in point[0]
            assert all(i.particles.get_vel() for i in path.phasepoints)

    def test_dump_phasepoint(self):
        """Test the method to dump a phase point."""
        system = create_system_for_lammps()
        with tempfile.TemporaryDirectory() as tempdir:
            shutil.copy(os.path.join(HERE, 'mocklammps.py'), tempdir)
            engine = LAMMPSEngine(os.path.join(HERE, 'mocklammps.py'),
                                  INPUT_PATH, 10, extra_files=None)
            engine.exe_dir = tempdir
            engine.dump_phasepoint(system, deffnm='dumped')
            # Check that the dumped file exists:
            assert os.path.isfile(os.path.join(tempdir, 'dumped.lammpstrj'))
