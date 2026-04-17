# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test initiate load methods."""
import logging
import os
import pathlib
import shutil
import tempfile
import pytest
from io import StringIO
from unittest.mock import patch
import numpy as np
from pyretis.core.box import create_box
from pyretis.core.path import Path
from pyretis.core.particles import Particles, ParticlesExt
from pyretis.core.pathensemble import PathEnsemble, PathEnsembleExt
from pyretis.core.random_gen import RandomGenerator
from pyretis.core.system import System
from pyretis.core.units import units_from_settings
from pyretis.engines.gromacs import GromacsEngine
from pyretis.engines.internal import Langevin
from pyretis.initiation.initiate_load import (
    clean_path,
    read_path_files,
    read_path_files_ext,
    reorderframes,
    _check_path,
    _do_the_dirty_load_job,
    _load_external_trajectory,
)
from pyretis.inout.common import make_dirs
from pyretis.inout.formats.path import PathExtFile, PathIntFile
from pyretis.inout.formats.order import OrderPathFile
from pyretis.inout.formats.energy import EnergyPathFile
from pyretis.inout.settings import (fill_up_tis_and_retis_settings,
                                    parse_settings_file)
from pyretis.orderparameter import Position, PositionVelocity
from pyretis.setup import create_simulation

logging.disable(logging.CRITICAL)
HERE = os.path.abspath(os.path.dirname(__file__))
CP2K = os.path.join(HERE, 'cp2k', 'mockcp2k.py')
GMX = os.path.join(HERE, 'gromacs', 'mockgmx.py')
MDRUN = os.path.join(HERE, 'gromacs', 'mockmdrun.py')
CP2K_DIR = os.path.join(HERE, 'cp2k', 'cp2k_input')
CP2K_LOAD = os.path.join(HERE, 'cp2k', 'pippo')
GMX_DIR = os.path.join(HERE, 'gromacs', 'gmx_input')
GMX_LOAD = os.path.join(HERE, 'gromacs', 'loader')
GMX_LOAD_TRJ = os.path.join(HERE, 'gromacs', 'load-trj')
INTERNAL_LOAD = os.path.join(HERE, 'internal', 'tis', 'loader')
INTERNAL_LOAD_2 = os.path.join(HERE, 'internal', 'retis-load', 'load-sparse')


def remove_file(filename):
    """Remove a file."""
    if os.path.isfile(filename):
        try:
            os.remove(filename)
        except OSError:
            pass


def set_up_system(order, pos_vel, vpot=None, ekin=None, internal=True):
    """Create a system for testing."""
    pos, vel = pos_vel
    system = System()
    if internal:
        system.particles = Particles(dim=3)
        for posi, veli in zip(pos, vel):
            system.add_particle(posi, vel=veli)
    else:
        system.particles = ParticlesExt(dim=3)
        system.add_particle(pos, vel=vel)
    system.order = order
    system.particles.vpot = vpot
    system.particles.ekin = ekin
    return system


def fill_path(path=None, npoints=20, scale=1.):
    """Fill in data for a arbitrary path."""
    if path is None:
        path = Path()
    for i in range(npoints):
        path.append(
            set_up_system([(-20.0 + 2.5 * i) * scale],
                          [[np.ones(3) * i],
                           [np.ones(3) * i]],
                          vpot=i, ekin=i,
                          internal=True)
        )
    return path


def find_files_in_dir(dirname):
    """Return a list of files in the given directory."""
    return [i.name for i in os.scandir(dirname) if i.is_file()]


def read_traj_order_energy(dirname, path_reader):
    """Read trajectory, order parameter and energy files."""
    traj, order, energy = [], [], []
    traj_file_name = os.path.join(dirname, 'traj.txt')
    with path_reader(traj_file_name, 'r') as trajfile:
        traj = list(next(trajfile.load())['data'])
    order_file_name = os.path.join(dirname, 'order.txt')
    with OrderPathFile(order_file_name, 'r') as orderfile:
        order = list(next(orderfile.load())['data'])
    energy_file_name = os.path.join(dirname, 'energy.txt')
    with EnergyPathFile(energy_file_name, 'r') as energyfile:
        for i in energyfile.load():
            energy = i['data']
            break
    return traj, order, energy


def path_traj_order_energy(path):
    """Get trajectory, order parameters and energies from a path."""
    traj, order, energy = [], [], {'time': [], 'vpot': [], 'ekin': []}
    for i, phasepoint in enumerate(path.phasepoints):
        vel = -1 if phasepoint.particles.get_vel() else 1
        traj.append([i, os.path.basename(phasepoint.particles.get_pos()[0]),
                     phasepoint.particles.get_pos()[1], vel])
        order.append(np.array([i, *phasepoint.order], dtype=float))
        energy['time'].append(i)
        energy['vpot'].append(phasepoint.particles.vpot)
        energy['ekin'].append(phasepoint.particles.ekin)
    for key, val in energy.items():
        energy[key] = np.array(val)
    return traj, order, energy


def path_traj_order_energy_int(path):
    """Get trajectory, order parameters and energies from a path."""
    traj, order, energy = [], [], {'time': [], 'vpot': [], 'ekin': []}
    for i, phasepoint in enumerate(path.phasepoints):
        traj.append({'pos': phasepoint.particles.get_pos(),
                     'vel': phasepoint.particles.get_vel()})
        order.append(np.array([i, *phasepoint.order]))
        energy['time'].append(i)
        energy['vpot'].append(phasepoint.particles.vpot)
        energy['ekin'].append(phasepoint.particles.ekin)
    for key, val in energy.items():
        energy[key] = np.array(val)
    return traj, order, energy


def compare_trajs(path_load, path):
    """Compare loaded trajectory with source."""
    traj, order, energy = read_traj_order_energy(path_load, PathExtFile)
    trajp, orderp, energyp = path_traj_order_energy(path)

    if not (len(order) == len(orderp)
            or len(traj) == len(trajp)):
        return False
    for i, j in zip(order, orderp):
        if not np.allclose(i, j):
            return False
    for traji, trajj in zip(traj, trajp):
        trajjj = [str(i) for i in trajj]
        for i, j in zip(traji, trajjj):
            if not i == j:
                return False
    for term in ['time', 'ekin', 'vpot']:
        if not len(energy[term]) == len(energyp[term]):
            return False
        if not np.allclose(energy[term], energyp[term]):
            return False
    return True


def compare_trajs_int(path_load, path):
    """Compare internal trajectory with source."""
    traj, order, energy = read_traj_order_energy(path_load, PathIntFile)
    trajp, orderp, energyp = path_traj_order_energy_int(path)
    status = True
    if not (len(order) == len(orderp)
            or len(traj) == len(trajp)):
        status = False
    for i, j in zip(order, orderp):
        if not np.allclose(i, j, atol=1e-6):
            status = False
    for traji, trajj in zip(traj, trajp):
        for i in ('pos', 'vel'):
            if not np.allclose(traji[i], trajj[i]):
                status = False
    term = 'time'
    if not len(energy[term]) == len(energyp[term]):
        status = False
    if not np.allclose(energy[term], energyp[term]):
        status = False
    for term in ('ekin', 'vpot'):
        if not len(energy[term]) == len(energyp[term]):
            status = False
        try:
            if not np.allclose(energy[term], energyp[term]):
                status = False
        except TypeError:  # Infinite errors are ok
            pass
    return status


class TestReadPathFilesExt:
    """Run the tests for the reading external trajectories."""

    def test_read_path_ext(self):
        """Test reading of external trajectories."""
        with tempfile.TemporaryDirectory() as tempdir:
            gmx_dir = os.path.join(tempdir, 'gmx_input')
            shutil.copytree(GMX_DIR, gmx_dir)
            load_dir = os.path.join(tempdir, 'load')
            shutil.copytree(GMX_LOAD, load_dir)
            path_ensemble = PathEnsembleExt(23, (-0.26, 0.02, 0.2),
                                            exe_dir=tempdir)
            path_load = os.path.join(
                load_dir, path_ensemble.ensemble_name_simple)
            for i in path_ensemble.directories():
                make_dirs(i)
            files = find_files_in_dir(path_ensemble.directory['accepted'])
            assert 0 == len(files)
            gmx = GromacsEngine(GMX, MDRUN, gmx_dir,
                                timestep=0.002,
                                subcycles=3,
                                exe_path=os.path.join(HERE, 'gromacs'),
                                maxwarn=1,
                                gmx_format='gro', write_vel=True,
                                write_force=False)
            gmx.exe_dir = path_ensemble.directory['accepted']
            path = Path(None)
            system = System(box=None)
            system.particles = ParticlesExt()
            ensemble = {'system': system,
                        'order_function': None,
                        'path_ensemble': path_ensemble,
                        'engine': gmx}
            with patch('sys.stdout', new=StringIO()):
                path_ok = read_path_files_ext(path, path_load, ensemble)
                assert path_ok == (True, 'ACC')
            files = find_files_in_dir(path_ensemble.directory['accepted'])
            assert 2 == len(files)
            assert 'trajF.trr' in files
            assert 'trajB.trr' in files
            cmpok = compare_trajs(path_load, path)
            assert cmpok

    def test_read_path_ext_no_order(self):
        """Test reading of external trajectories with missing order.txt."""
        orderparameter = PositionVelocity(1472, dim='x', periodic=True)
        with tempfile.TemporaryDirectory() as tempdir:
            gmx_dir = os.path.join(tempdir, 'gmx_input')
            shutil.copytree(GMX_DIR, gmx_dir)
            load_dir = os.path.join(tempdir, 'load')
            shutil.copytree(GMX_LOAD, load_dir)
            path_ensemble = PathEnsembleExt(23, [1.78, 1.795, 1.80],
                                            exe_dir=tempdir)
            path_load = os.path.join(
                load_dir, path_ensemble.ensemble_name_simple)
            remove_file(os.path.join(path_load, 'order.txt'))
            for i in path_ensemble.directories():
                make_dirs(i)
            files = find_files_in_dir(path_ensemble.directory['accepted'])
            assert 0 == len(files)
            files = find_files_in_dir(path_load)
            assert 2 == len(files)
            assert 'energy.txt' in files
            assert 'traj.txt' in files
            assert 'order.txt' not in files
            gmx = GromacsEngine(GMX, MDRUN, gmx_dir,
                                timestep=0.002,
                                subcycles=3,
                                exe_path=os.path.join(HERE, 'gromacs'),
                                maxwarn=1,
                                gmx_format='gro', write_vel=True,
                                write_force=False)
            gmx.exe_dir = path_ensemble.directory['accepted']
            path = Path(None)
            system = System(box=None)
            system.particles = ParticlesExt()
            ensemble = {'system': system,
                        'path_ensemble': path_ensemble,
                        'order_function': orderparameter,
                        'engine': gmx}
            with patch('sys.stdout', new=StringIO()):
                path_ok = read_path_files_ext(path, path_load, ensemble)
                assert path_ok == (True, 'ACC')
            files = find_files_in_dir(path_ensemble.directory['accepted'])
            assert 2 == len(files)
            assert 'trajF.trr' in files
            assert 'trajB.trr' in files

    def test_read_path_ext_no_energy(self):
        """Test reading of external trajectories with missing order.txt."""
        with tempfile.TemporaryDirectory() as tempdir:
            gmx_dir = os.path.join(tempdir, 'gmx_input')
            shutil.copytree(GMX_DIR, gmx_dir)
            load_dir = os.path.join(tempdir, 'load')
            shutil.copytree(GMX_LOAD, load_dir)
            path_ensemble = PathEnsembleExt(23, (-0.26, 0.02, 0.2),
                                            exe_dir=tempdir)
            path_load = os.path.join(
                load_dir, path_ensemble.ensemble_name_simple)
            remove_file(os.path.join(path_load, 'energy.txt'))
            for i in path_ensemble.directories():
                make_dirs(i)
            files = find_files_in_dir(path_ensemble.directory['accepted'])
            assert 0 == len(files)
            files = find_files_in_dir(path_load)
            assert 2 == len(files)
            assert 'energy.txt' not in files
            assert 'traj.txt' in files
            assert 'order.txt' in files
            gmx = GromacsEngine(GMX, MDRUN, gmx_dir,
                                timestep=0.002,
                                subcycles=3,
                                exe_path=os.path.join(HERE, 'gromacs'),
                                maxwarn=1,
                                gmx_format='gro', write_vel=True,
                                write_force=False)
            gmx.exe_dir = path_ensemble.directory['accepted']
            path = Path(None)
            system = System(box=None)
            system.particles = ParticlesExt()
            ensemble = {'system': system,
                        'path_ensemble': path_ensemble,
                        'order_function': None,
                        'engine': gmx}
            with patch('sys.stdout', new=StringIO()):
                path_ok = read_path_files_ext(path, path_load, ensemble)
                assert path_ok == (True, 'ACC')
            files = find_files_in_dir(path_ensemble.directory['accepted'])
            assert 2 == len(files)
            assert 'trajF.trr' in files
            assert 'trajB.trr' in files
            for i in path.phasepoints:
                assert i.particles.ekin is None
                assert i.particles.vpot is None

    def test_check_path(self):
        """Test checking of the path."""
        with tempfile.TemporaryDirectory() as tempdir:
            gmx_dir = os.path.join(tempdir, 'gmx_input')
            shutil.copytree(GMX_DIR, gmx_dir)
            load_dir = os.path.join(tempdir, 'load')
            shutil.copytree(GMX_LOAD, load_dir)
            path_ensemble = PathEnsembleExt(23, [-0.26, 0.02, 0.2],
                                            exe_dir=tempdir)
            path_load = os.path.join(
                load_dir, path_ensemble.ensemble_name_simple)
            for i in path_ensemble.directories():
                make_dirs(i)
            gmx = GromacsEngine(GMX, MDRUN, gmx_dir,
                                timestep=0.002,
                                subcycles=3,
                                exe_path=os.path.join(HERE, 'gromacs'),
                                maxwarn=1,
                                gmx_format='gro', write_vel=True,
                                write_force=False)
            gmx.exe_dir = path_ensemble.directory['accepted']
            path = Path(None)
            system = System(box=None)
            system.particles = ParticlesExt()
            ensemble = {'system': system,
                        'path_ensemble': path_ensemble,
                        'order_function': None,
                        'engine': gmx}
            with patch('sys.stdout', new=StringIO()):
                read_path_files_ext(path,
                                    path_load, ensemble)
            accept, status = _check_path(path, path_ensemble)
            assert accept
            assert status == 'ACC'
            path_ensemble = PathEnsembleExt(23, [-0.26, 0.2, 0.4],
                                            exe_dir=tempdir)
            accept, status = _check_path(path, path_ensemble)
            assert not accept
            assert status == 'EWI'
            path_r = path.reverse(
                Position(0, dim='x', periodic=False)
            )
            path_ensemble = PathEnsembleExt(23, [-0.26, 0.02, 0.2],
                                            exe_dir=tempdir)
            accept, status = _check_path(path_r, path_ensemble)
            assert not accept
            assert status == 'SWI'
            path.append(
                set_up_system([-0.26, 0.0], [('somefile.trr', 0), True],
                              internal=False)
            )
            path_ensemble = PathEnsembleExt(23, [-0.26, 0.3, 0.3],
                                            exe_dir=tempdir)
            accept, status = _check_path(path, path_ensemble)
            assert not accept
            assert status == 'NCR'


class TestReadPathFiles:
    """Run the tests for reading internal trajectories."""

    def test_read_path(self):
        """Test reading of internal trajectories."""
        with tempfile.TemporaryDirectory() as tempdir:
            load_dir = os.path.join(tempdir, 'load')
            shutil.copytree(INTERNAL_LOAD, load_dir)
            path_ensemble = PathEnsemble(1, [-0.9, -0.8, -0.7],
                                         exe_dir=tempdir)
            path_load = os.path.join(
                load_dir, path_ensemble.ensemble_name_simple)
            path = Path(None)
            engine = Langevin(0.002, 0.3)
            box = create_box(periodic=[False])
            system = System(box=box)
            system.particles = Particles()
            ensemble = {'system': system,
                        'path_ensemble': path_ensemble,
                        'order_function': None,
                        'engine': engine}
            with patch('sys.stdout', new=StringIO()):
                accept, status = read_path_files(path, path_load, ensemble)
                assert accept
                assert status == 'ACC'
            cmpt = compare_trajs_int(path_load, path)
            assert cmpt

    def test_read_path_noorder(self):
        """Test reading of internal trajectories with missing order.txt."""
        orderparameter = PositionVelocity(0, dim='x', periodic=False)
        with tempfile.TemporaryDirectory() as tempdir:
            load_dir = os.path.join(tempdir, 'load')
            shutil.copytree(INTERNAL_LOAD, load_dir)
            path_ensemble = PathEnsemble(1, [-0.9, -0.8, -0.7],
                                         exe_dir=tempdir)
            path_load = os.path.join(
                load_dir, path_ensemble.ensemble_name_simple)
            remove_file(os.path.join(path_load, 'order.txt'))
            path = Path(None)
            engine = Langevin(0.002, 0.3)
            box = create_box(periodic=[False])
            system = System(box=box)
            system.particles = Particles(dim=1)
            ensemble = {'system': system,
                        'path_ensemble': path_ensemble,
                        'order_function': orderparameter,
                        'engine': engine}
            with patch('sys.stdout', new=StringIO()):
                accept, status = read_path_files(path, path_load, ensemble)
                assert accept
                assert status == 'ACC'
            shutil.copyfile(
                os.path.join(INTERNAL_LOAD, '001', 'order.txt'),
                os.path.join(path_load, 'order.txt')
            )
            cmpt = compare_trajs_int(path_load, path)
            assert cmpt

    def test_read_path_noenergy(self):
        """Test reading of internal trajectories with missing energy.txt."""
        with tempfile.TemporaryDirectory() as tempdir:
            load_dir = os.path.join(tempdir, 'load')
            shutil.copytree(INTERNAL_LOAD, load_dir)
            path_ensemble = PathEnsemble(1, [-0.9, -0.8, -0.7],
                                         exe_dir=tempdir)
            path_load = os.path.join(
                load_dir, path_ensemble.ensemble_name_simple)
            remove_file(os.path.join(path_load, 'energy.txt'))
            path = Path(None)
            engine = Langevin(0.002, 0.3)
            box = create_box(periodic=[False])
            system = System(box=box)
            system.particles = Particles(dim=1)
            ensemble = {'system': system,
                        'order_function': None,
                        'path_ensemble': path_ensemble,
                        'engine': engine}
            with patch('sys.stdout', new=StringIO()):
                accept, status = read_path_files(path, path_load, ensemble)
                assert accept
                assert status == 'ACC'
            for i in path.phasepoints:
                assert i.particles.ekin is None
                assert i.particles.vpot is None


class TestInitiateLoad:
    """Test the full initiate load for a TIS simulation."""

    @staticmethod
    def _set_up_for_load_int(settings, tempdir):
        """Set up a system for testing the initiate load method."""
        settings['simulation']['exe_path'] = os.path.join(tempdir)
        load_dir = os.path.join(tempdir, 'loader')

        shutil.copytree(os.path.join(HERE, 'internal', 'tis'), tempdir)

        settings['initial-path']['load_folder'] = load_dir
        settings['engine']['exe_path'] = tempdir
        settings['engine']['input_path'] = tempdir

        del settings['ensemble']
        fill_up_tis_and_retis_settings(settings)

        units_from_settings(settings)
        with patch('sys.stdout', new=StringIO()):
            simulation = create_simulation(settings)
        simulation.set_up_output(settings)
        return simulation, load_dir

    @staticmethod
    def _test_dirs(tempdir):
        """Test the creation of the data structure."""
        for test1 in ['000', '001', '002', '003']:
            dirt = os.path.join(tempdir, test1)
            if not os.path.exists(dirt):
                return False
            for test2 in ['traj', 'generate', 'accepted']:
                if not os.path.exists(os.path.join(dirt, test2)):
                    return False
            for test3 in ['order.txt', 'energy.txt', 'pathensemble.txt']:
                if not os.path.exists(os.path.join(tempdir, test1, test3)):
                    return False
        return True

    @staticmethod
    def _set_up_for_load(settings, tempdir, copy_load=True):
        """Set up a system for testing the initiate load method."""
        settings['simulation']['exe_path'] = tempdir
        load_dir = os.path.join(
            tempdir, settings['initial-path']['load_folder']
        )
        settings['initial-path']['load_folder'] = load_dir
        ens_name = ['000', '001', '002', '003', '004']
        for i in range(len(settings['ensemble'])):
            os.makedirs(os.path.join(tempdir, ens_name[i]))
            settings['ensemble'][i]['simulation']['exe_path'] = tempdir
            settings['ensemble'][i]['initial-path']['load_folder'] =\
                load_dir
            settings['ensemble'][i]['engine']['exe_path'] =\
                os.path.join(HERE, 'internal', 'tis')
            settings['ensemble'][i]['engine']['input_path'] =\
                os.path.join(HERE, 'internal', 'tis')
        if copy_load:
            shutil.copytree(INTERNAL_LOAD, load_dir)
        shutil.copyfile(
            os.path.join(HERE, 'internal', 'tis', 'initial.xyz'),
            os.path.join(tempdir, 'initial.xyz')
        )

        units_from_settings(settings)
        with patch('sys.stdout', new=StringIO()):
            simulation = create_simulation(settings)
        simulation.set_up_output(settings)
        return simulation

    @staticmethod
    def _set_up_for_load_2(settings, tempdir, copy_load=True):
        """Set up a second system for testing the initiate load method."""
        settings['simulation']['exe_path'] = tempdir
        load_dir = os.path.join(
            tempdir, settings['initial-path']['load_folder'])

        settings['engine']['gmx'] = GMX
        settings['engine']['mdrun'] = MDRUN
        settings['engine']['exe_path'] = tempdir
        settings['engine']['input_path'] = tempdir
        settings['simulation']['exe_path'] = tempdir
        settings['initial-path']['load_folder'] = load_dir

        del settings['ensemble']
        fill_up_tis_and_retis_settings(settings)

        if copy_load:
            shutil.copytree(INTERNAL_LOAD_2, load_dir)
        units_from_settings(settings)
        simulation = create_simulation(settings)
        simulation.set_up_output(settings)
        return simulation

    def test_internal_load(self):
        """Test the load initialise for internal trajectories."""
        settings = parse_settings_file(os.path.join(HERE, 'internal', 'tis',
                                                    'tis.rst'))

        with tempfile.TemporaryDirectory() as tempdirbase:
            tempdir = os.path.join(tempdirbase, 'internal', 'tis')
            with patch('sys.stdout', new=StringIO()):
                simulation, _ = \
                    self._set_up_for_load_int(settings, tempdir)
                init = simulation.initiate(settings)
                assert init
            # We expect now that we have loaded for 001:
            assert len(simulation.ensembles) == 1
            assert simulation.ensembles[0]['path_ensemble'] \
                             .ensemble_number == 1
            # Check that the path contains what it should:
            path = simulation.ensembles[0]['path_ensemble'].last_path
            # Check that we created an ensemble dir
            ensemble_dir = os.path.join(tempdir, '001')
            assert os.path.isdir(ensemble_dir)
            # All memory is dumped to file.
            assert compare_trajs_int(ensemble_dir, path)
            # Just delete simulation to force closing of files:
            del simulation
            assert compare_trajs_int(ensemble_dir, path)

    def test_fail_wrong_engine(self):
        """Test that the test fail for an unsupported engine type."""
        settings = parse_settings_file(os.path.join(HERE, 'internal', 'tis',
                                                    'tis.rst'))
        with tempfile.TemporaryDirectory() as tempdir:
            simulation = self._set_up_for_load(settings, tempdir)
            ensemble = simulation.ensembles[0]
            ensemble['engine'].engine_type = 'this-should-fail!'
            with patch('sys.stdout', new=StringIO()):
                with pytest.raises(ValueError):
                    simulation.initiate(settings)

    def test_internal_load_sparse(self):
        """Test the load initialise for internal trajectories."""
        settings = parse_settings_file(os.path.join(HERE, 'internal',
                                                    'retis-load',
                                                    'reload.rst'))
        with tempfile.TemporaryDirectory() as tempdirbase:
            tempdir = os.path.join(tempdirbase, 'internal', 'retis-load')

            with patch('sys.stdout', new=StringIO()):
                simulation = self._set_up_for_load_2(settings, tempdir)
                init = simulation.initiate(settings)
                assert init

            paths = []
            for i in range(3):
                paths.append(
                    simulation.ensembles[i]['path_ensemble'].last_path)
            del simulation
            assert compare_trajs_int(os.path.join(tempdir, '000'),
                                     paths[0])

            assert os.path.exists(os.path.join(tempdir, '000'))
            assert os.path.exists(os.path.join(tempdir, '001',
                                                        'traj'))
            assert os.path.exists(os.path.join(tempdir, '002',
                                                        'generate'))
            assert os.path.exists(os.path.join(tempdir, '001',
                                                        'pathensemble.txt'))
            assert os.path.exists(os.path.join(tempdir, '002',
                                                        'order.txt'))
            assert paths[0].status == 'ACC'
            for i in range(1, 2):
                assert paths[i].status == 'ACC'

    def test_internal_explorer(self):
        """Test the explorer for internal trajectories."""
        settings = parse_settings_file(os.path.join(HERE, 'internal',
                                                    'retis-load',
                                                    'reload.rst'))
        settings['simulation']['task'] = 'explore'
        with tempfile.TemporaryDirectory() as tempdirbase:
            tempdir = os.path.join(tempdirbase, 'internal', 'retis-load')

            with patch('sys.stdout', new=StringIO()):
                simulation = self._set_up_for_load_2(settings, tempdir)
                init = simulation.initiate(settings)
                assert init

                # Does it work too?
                res = simulation.step()
                assert res['status-1'] == 'EXP'
                assert res['status-2'] == 'EXP'

    def test_internal_explorer_bad_init(self):
        """Test the explorer for internal trajectories."""
        settings = parse_settings_file(os.path.join(HERE, 'internal',
                                                    'retis-load',
                                                    'reload.rst'))
        settings['simulation']['task'] = 'explore'
        # A case where the initial configurations do not have frames after the
        # ensemble interface.
        with tempfile.TemporaryDirectory() as tempdirbase:
            tempdir = os.path.join(tempdirbase, 'internal', 'retis-load')
            with patch('sys.stdout', new=StringIO()):
                settings['tis']['maxlength'] = 100
                settings['simulation']['interfaces'][0] = -2
                settings['simulation']['interfaces'].append(1000)
                simulation = self._set_up_for_load_2(settings, tempdir)
                init = simulation.initiate(settings)
                assert init

                # Does it work too?
                res = simulation.step()
                assert res['status-0'] == 'EXP'
                # It shall try anyhow.
                assert res['status-2'] == 'EXP'

    def test_internal_load_sparse_and_kick(self):
        """Test the load initialise that fixes paths."""
        settings = parse_settings_file(os.path.join(HERE, 'internal',
                                                    'retis-load',
                                                    'reload.rst'))
        settings['simulation']['interfaces'] = [-1.2, -0.4, -0.2]
        # Load poor initial paths
        with tempfile.TemporaryDirectory() as tempdirbase:
            tempdir = os.path.join(tempdirbase, 'internal', 'retis-load')

            with patch('sys.stdout', new=StringIO()):
                simulation = self._set_up_for_load_2(settings, tempdir)
                # Check that we fail first
                with pytest.raises(ValueError) as cm:
                    simulation.initiate(settings)
                assert "load_and_kick" in str(cm.value)

        settings = parse_settings_file(os.path.join(HERE, 'internal',
                                                    'retis-load',
                                                    'reload.rst'))
        settings['initial-path']['load_and_kick'] = True
        del settings['ensemble']
        settings['simulation']['interfaces'] = [-1.2, -0.4, -0.2]
        fill_up_tis_and_retis_settings(settings)

        with tempfile.TemporaryDirectory() as tempdirbase:
            tempdir = os.path.join(tempdirbase, 'internal', 'retis-load')

            with patch('sys.stdout', new=StringIO()):
                simulation2 = self._set_up_for_load_2(settings, tempdir)
                init = simulation2.initiate(settings)
                assert init

    def test_initiate_load_fail(self):
        """Test that the initiate load fails for non-existing load folder."""
        settings = parse_settings_file(
            os.path.join(HERE, 'internal', 'tis', 'tis.rst')
        )
        with tempfile.TemporaryDirectory() as tempdir:
            with patch('sys.stdout', new=StringIO()):
                simulation = self._set_up_for_load(settings, tempdir,
                                                   copy_load=False)
                with pytest.raises(FileNotFoundError):
                    simulation.initiate(settings)

    def test_initiate_load_fail2(self):
        """Test that the initiate load fails for insufficient file inputs."""
        settings = parse_settings_file(os.path.join(HERE, 'internal',
                                                    'retis-load',
                                                    'reload.rst'))
        with tempfile.TemporaryDirectory() as tempdir:
            with patch('sys.stdout', new=StringIO()):
                simulation = self._set_up_for_load_2(settings, tempdir)
            # Instead to change the inputs, we change the requirements.
            simulation.ensembles[0]['path_ensemble'].interfaces = [1,
                                                                   100,
                                                                   100000]
            with patch('sys.stdout', new=StringIO()):
                with pytest.raises(ValueError) as err:
                    _ = simulation.initiate(settings)
                assert '000 does not satisfy the relative ensemble' \
                    in str(err.value)

    def test_reorderframes(self):
        """Test the reorderframes function."""
        interfaces = [float('-inf'), -10., 10]
        path_ensemble = PathEnsemble(0, interfaces)
        path = fill_path()

        left, _, right = path_ensemble.interfaces
        path_ensemble.start_condition = 'R'

        # Test that the re-order function generated path
        # satisfies the ensemble requirements (R, R case)
        path, _ = reorderframes(path, path_ensemble)
        assert path.ordermin[0] < interfaces[1]
        assert path.ordermax[0] > interfaces[1]
        assert path.get_start_point(left, right) == 'R'
        assert path.get_end_point(left, right) == 'R'

        # Test that the re-order function generated path
        # satisfies the ensemble requirements (L, L case)
        path_ensemble = PathEnsemble(1, [-18., -10., 60.])
        left, _, right = path_ensemble.interfaces
        path_ensemble.start_condition = 'L'
        path = fill_path()
        path, _ = reorderframes(path, path_ensemble)

        assert path.get_start_point(left, right) == 'L'
        assert path.get_end_point(left, right) == 'L'
        assert path.ordermin[0] < interfaces[1]
        assert path.ordermax[0] > interfaces[1]

        # Test that the re-order function generated path
        # satisfies the ensemble requirements (L, R case)
        path_ensemble = PathEnsemble(2, [-10., 1., 2.])
        left, _, right = path_ensemble.interfaces
        path_ensemble.start_condition = 'L'
        path = fill_path()
        path, _ = reorderframes(path, path_ensemble)

        assert path.get_start_point(left, right) == 'L'
        assert not path.get_end_point(left, right) == 'L'

    def test_clean_path(self):
        """Test the clean_path function."""
        interfaces = [-14., -10., 10]
        path_ensemble = PathEnsemble(0, interfaces)

        path = fill_path()
        clean_path(path, path_ensemble)
        assert path.length == 6
        assert _check_path(path, path_ensemble)

        # Same results?
        clean_path(path, path_ensemble)
        assert path.length == 6
        assert _check_path(path, path_ensemble)

        # Some intuitive numerical checks: only crossing points
        # shall remain (4 for the 0 interface and one for the B
        # interface).
        path = fill_path(scale=10000)
        interfaces = [-16., 0., 16]
        path_ensemble = PathEnsemble(0, interfaces)
        initial_length = path.length
        path, _ = reorderframes(path, path_ensemble)

        # Numerical checks
        counter_r, counter_l, counter_m = 0, 0, 0
        for _, phasepoint in enumerate(path.phasepoints):
            if phasepoint.order[0] > interfaces[2]:
                counter_r += 1
            if phasepoint.order[0] < interfaces[0]:
                counter_l += 1
            if interfaces[0] < phasepoint.order[0] < interfaces[2]:
                counter_m += 1

        assert path.length == (counter_r + counter_l + counter_m)
        assert (counter_r, counter_m, counter_l) == (11, 1, 8)
        assert path.length == initial_length
        assert path.length == 20
        assert _check_path(path, path_ensemble)

        # This creates a bad path, invalid start but some ok frames
        path = fill_path(scale=0.1)
        interfaces = [-16., 0., 16]
        path_ensemble = PathEnsemble(0, interfaces)
        initial_length = path.length
        path, _ = reorderframes(path, path_ensemble)
        clean_path(path, path_ensemble)

        # Numerical checks
        counter_r, counter_l, counter_m = 0, 0, 0
        for phasepoint in path.phasepoints:
            if phasepoint.order[0] >= interfaces[2]:
                counter_r += 1
            if phasepoint.order[0] < interfaces[0]:
                counter_l += 1
            if interfaces[0] <= phasepoint.order[0] < interfaces[2]:
                counter_m += 1

        assert path.length == (counter_r + counter_l + counter_m)
        assert (counter_r, counter_m, counter_l) == (0, 21, 0)
        assert path.length > initial_length
        assert path.length == 21
        assert _check_path(path, path_ensemble)

        # This creates a good path, bookkeeping case for 0 ens
        path = fill_path(scale=0.1)
        interfaces = [-1., 0., 16]
        path_ensemble = PathEnsemble(0, interfaces)
        initial_length = path.length
        path, _ = reorderframes(path, path_ensemble)
        clean_path(path, path_ensemble)

        # Numerical checks
        counter_r, counter_l, counter_m = 0, 0, 0
        for phasepoint in path.phasepoints:
            if phasepoint.order[0] >= interfaces[2]:
                counter_r += 1
            if phasepoint.order[0] < interfaces[0]:
                counter_l += 1
            if interfaces[0] <= phasepoint.order[0] < interfaces[2]:
                counter_m += 1

        assert path.length == (counter_r + counter_l + counter_m)
        assert (counter_r, counter_m, counter_l) == (0, 14, 0)
        assert path.length < initial_length
        assert path.length == 14
        assert _check_path(path, path_ensemble)

    def test_clean_path_2(self):
        """Test the clean_path function."""
        interfaces = [float('-inf'), -10., -5]
        path_ensemble = PathEnsemble(0, interfaces)
        rgen = RandomGenerator(seed=0)
        path = Path(rgen, maxlen=1000)

        fill_path(path, npoints=20000, scale=0.01)
        clean_path(path, path_ensemble)
        assert path.length == 5
        assert _check_path(path, path_ensemble)

        interfaces = [-14., -10., 10]
        path_ensemble = PathEnsemble(0, interfaces)
        rgen = RandomGenerator(seed=0)
        path = Path(rgen, maxlen=1000)

        fill_path(path)
        clean_path(path, path_ensemble)
        assert path.length == 6
        assert _check_path(path, path_ensemble)

        # Check that it doesn't fix the impossible.
        # Too large interfaces, all should be kept
        interfaces = [-14000., 0., 1000]
        path_ensemble = PathEnsemble(0, interfaces)

        path = fill_path(npoints=1000, scale=0.16)
        path, _ = reorderframes(path, path_ensemble)
        assert _check_path(path, path_ensemble) == \
            (False, 'EWI')
        clean_path(path, path_ensemble)
        assert path.length == 21

        # Too large orderp, almost nothing should remain
        interfaces = [-16., 0., 16]
        path_ensemble = PathEnsemble(0, interfaces)
        path = fill_path(npoints=1000, scale=16016)
        path, _ = reorderframes(path, path_ensemble)
        clean_path(path, path_ensemble)
        assert _check_path(path, path_ensemble) == (False, 'EWI')
        assert path.length == 5

    def test_clean_path3(self):
        """Test the clean_path function."""
        interfaces = [float('-inf'), -1., 5]
        path_ensemble = PathEnsemble(0, interfaces)
        rgen = RandomGenerator(seed=0)
        path = Path(rgen, maxlen=1000)
        fill_path(path, npoints=20, scale=0.01)
        path_ensemble.start_condition = ('L', 'R')
        clean_path(path, path_ensemble)
        assert path.length == 20
        assert _check_path(path, path_ensemble)

        interfaces = [-1., 0., 1]
        path_ensemble = PathEnsemble(0, interfaces)
        rgen = RandomGenerator(seed=0)
        path = Path(rgen, maxlen=1000)
        fill_path(path, npoints=16, scale=0.1)
        path_ensemble.start_condition = ('L', 'R')
        assert path.phasepoints[0].order[0] <= interfaces[0]
        assert path.phasepoints[path.ordermax[1]].order[0] > \
            interfaces[-1]
        clean_path(path, path_ensemble)
        assert path.phasepoints[-1].order[0] <= interfaces[0]
        assert path.length == 11

        interfaces = [-1., 0., 10]
        path_ensemble = PathEnsemble(0, interfaces)
        rgen = RandomGenerator(seed=0)
        path = Path(rgen, maxlen=1000)
        fill_path(path, npoints=160, scale=0.01)
        path_ensemble.start_condition = ('L', 'R')
        clean_path(path, path_ensemble)
        assert not path.phasepoints[0].order[0] <= interfaces[0]
        assert not path.phasepoints[-1].order[0] > interfaces[-1]

    def test_clean_reorder_path_repptis(self):
        """Test clean path and reorderframes for pptis ensembles"""
        # 1) body ensemble [i^+-] with i > 0
        # 1A) too long LMR path --> Trim down to LMR
        interfaces = [-1., 0, 1.]
        path_ensemble = PathEnsemble(3, interfaces)
        rgen = RandomGenerator(seed=0)
        path = Path(rgen, maxlen=1000)
        fill_path(path, npoints=20, scale=0.055)
        path_ensemble.start_condition = ('L', 'R')
        clean_path(path, path_ensemble, simtype='repptis')
        assert path.length == 17
        assert _check_path(path, path_ensemble)
        assert path.phasepoints[0].order[0] < interfaces[0]
        assert path.phasepoints[-1].order[0] > interfaces[-1]
        # 1B) LM* path --> Makes an LML path
        interfaces = [-1., 0, 2.]
        path_ensemble = PathEnsemble(3, interfaces)
        rgen = RandomGenerator(seed=0)
        path = Path(rgen, maxlen=1000)
        fill_path(path, npoints=20, scale=0.055)
        path_ensemble.start_condition = ('L', 'R')
        clean_path(path, path_ensemble, simtype='repptis')
        pathnew, (accept, status) = reorderframes(path, path_ensemble)
        assert accept
        assert status == 'ACC'
        assert pathnew.phasepoints[0].order[0] < interfaces[0]
        assert pathnew.phasepoints[-1].order[0] < interfaces[0]
        assert pathnew.ordermax[0] > interfaces[1]
        assert _check_path(pathnew, path_ensemble)[0]
        assert not _check_path(path, path_ensemble)[0]
        # 2) ensemble [0^+-'] should only give LMR, RML, or LML paths
        interfaces = [-1, -1, 0]
        path_ensemble = PathEnsemble(1, interfaces)
        rgen = RandomGenerator(seed=0)
        path = Path(rgen, maxlen=1000)
        fill_path(path, npoints=20, scale=0.055)
        path_ensemble.start_condition = ('L', 'R')
        clean_path(path, path_ensemble, simtype='repptis')
        assert path.phasepoints[0].order[0] < interfaces[0]
        assert path.phasepoints[-1].order[0] > interfaces[-1]
        assert _check_path(path, path_ensemble)
        # 2b) what if start_condition is only L
        interfaces = [-1, -1, 0]
        path_ensemble = PathEnsemble(1, interfaces)
        rgen = RandomGenerator(seed=0)
        path = Path(rgen, maxlen=1000)
        fill_path(path, npoints=20, scale=0.055)
        path_ensemble.start_condition = 'L'
        clean_path(path, path_ensemble, simtype='repptis')
        assert path.phasepoints[0].order[0] < interfaces[0]
        assert path.phasepoints[-1].order[0] < interfaces[0]
        assert _check_path(path, path_ensemble)

    def test_external_load(self):
        """Test the load initialise for external trajectories."""
        settings = parse_settings_file(
            os.path.join(HERE, 'gromacs', 'tis.rst')
        )
        settings['engine']['gmx'] = GMX
        settings['engine']['mdrun'] = MDRUN

        with tempfile.TemporaryDirectory() as tempdirbase:
            tempdir = os.path.join(tempdirbase, 'gromacs')
            settings['simulation']['exe_path'] = tempdir
            gmx_dir = os.path.join(tempdir, 'gmx_input')
            shutil.copytree(GMX_DIR, gmx_dir)
            load_dir = os.path.join(tempdir, 'loader')
            shutil.copytree(GMX_LOAD, load_dir)
            settings['engine']['input_path'] = gmx_dir
            units_from_settings(settings)
            settings['initial-path']['load_folder'] = load_dir
            settings['engine']['exe_path'] = tempdir

            with patch('sys.stdout', new=StringIO()):
                # Reset inputs for this example set up in the tempdir
                del settings['ensemble']
                fill_up_tis_and_retis_settings(settings)
                # Calculate the WF path weight:
                settings['ensemble'][0]['tis']['high_accept'] = True
                settings['ensemble'][0]['tis']['shooting_move'] = 'wf'
                simulation = create_simulation(settings)
                simulation.set_up_output(settings)
                init = simulation.initiate(settings)
                assert init
            for to_test in ['accepted', 'order.txt',
                            'energy.txt', 'pathensemble.txt']:
                assert os.path.exists(os.path.join(tempdir, '023',
                                                   to_test))
            assert simulation.ensembles[0]['path_ensemble'] \
                             .last_path.status == 'ACC'
            assert simulation.ensembles[0]['path_ensemble'] \
                             .last_path.weight == 42.0

    def test_external_load_cp2k(self):
        """Test the load initialise for external trajectories with cp2k."""
        settings = parse_settings_file(
            os.path.join(HERE, 'cp2k', 'retis.rst'))
        settings['engine']['cp2k'] = CP2K

        with tempfile.TemporaryDirectory() as tempdirbase:
            tempdir = os.path.join(tempdirbase, 'cp2k')
            cp2k_dir = os.path.join(tempdir, 'cp2k_input')
            shutil.copytree(CP2K_DIR, cp2k_dir)
            load_dir = os.path.join(tempdir, 'pippo')
            shutil.copytree(CP2K_LOAD, load_dir)
            settings['engine']['input_path'] = cp2k_dir
            units_from_settings(settings)
            settings['initial-path']['load_folder'] = load_dir
            settings['simulation']['exe_path'] = tempdir
            settings['engine']['exe_path'] = tempdir
            with patch('sys.stdout', new=StringIO()):
                # Reset inputs for this example set up in the tempdir
                del settings['ensemble']
                fill_up_tis_and_retis_settings(settings)
                for i in ['000', '001', '002', '003']:
                    edir = os.path.join(load_dir, i)
                    _do_the_dirty_load_job(load_dir, edir)
                simulation = create_simulation(settings)
            simulation.set_up_output(settings)
            with patch('sys.stdout', new=StringIO()):
                init = simulation.initiate(settings)
                assert init
            assert self._test_dirs(tempdir)

            for i in range(3):
                assert simulation.ensembles[i]['path_ensemble'] \
                    .last_path.status == 'ACC'
                path = simulation.ensembles[i]['path_ensemble'].last_path
                for j in range(int(path.length/3)):
                    point1 = path.phasepoints[j]
                    point2 = path.phasepoints[-1-j]
                    assert point1.order[0] == point2.order[0]

    def test_external_load_external(self):
        """Test the load initialise for external trajectories."""
        settings = parse_settings_file(
            os.path.join(HERE, 'gromacs', 'auto-load.rst'))
        orderp_file = os.path.join(HERE, 'gromacs', 'orderp.py')
        settings['engine']['gmx'] = GMX
        settings['engine']['mdrun'] = MDRUN
        top_file = settings['initial-path']['top_file']

        with tempfile.TemporaryDirectory() as tempdirbase:
            tempdir = os.path.join(tempdirbase, 'gromacs')
            settings['simulation']['exe_path'] = tempdir
            gmx_dir = os.path.join(tempdir, 'gmx_input')
            shutil.copytree(GMX_DIR, gmx_dir)
            shutil.copy(orderp_file, tempdir)
            tmp_top_file = os.path.join(tempdir, top_file)
            load_dir = os.path.join(tempdir, 'loader')
            shutil.copytree(GMX_LOAD, load_dir)
            settings['engine']['input_path'] = gmx_dir

            units_from_settings(settings)
            settings['initial-path']['load_folder'] = load_dir
            settings['initial-path']['top_file'] = tmp_top_file
            settings['engine']['exe_path'] = tempdir

            # Reset inputs for this example set up in the tempdir
            del settings['ensemble']
            fill_up_tis_and_retis_settings(settings)

            with patch('sys.stdout', new=StringIO()):
                simulation = create_simulation(settings)
                simulation.set_up_output(settings)
                init = simulation.initiate(settings)
                assert init
            assert self._test_dirs(tempdir)

            for i in range(4):
                assert simulation.ensembles[i]['path_ensemble'] \
                   .last_path.status == 'ACC'

            assert os.path.exists(os.path.join(load_dir, '004',
                                               'traj.txt'))

            for i in range(1, 4):
                path = simulation.ensembles[i]['path_ensemble'].last_path
                for j in range(int(path.length / 2) - 1):
                    point1 = path.phasepoints[j]
                    point2 = path.phasepoints[j + 1]
                    assert point1.order[0] <= point2.order[0]

    def test_external_load_external_2(self):
        """Test the load initialise for external trajectories."""
        settings = parse_settings_file(
            os.path.join(HERE, 'gromacs', 'auto-load.rst'))
        orderp_file = os.path.join(HERE, 'gromacs', 'orderp.py')
        settings['simulation']['interfaces'] = [2.78, 2.88, 3.01, 3.80, 14.00]
        settings['engine']['gmx'] = GMX
        settings['engine']['mdrun'] = MDRUN
        top_file = settings['initial-path']['top_file']

        with tempfile.TemporaryDirectory() as tempdirbase:
            tempdir = os.path.join(tempdirbase, 'gromacs')
            settings['simulation']['exe_path'] = tempdir
            gmx_dir = os.path.join(tempdir, 'gmx_input')
            shutil.copytree(GMX_DIR, gmx_dir)
            shutil.copy(orderp_file, tempdir)
            tmp_top_file = os.path.join(tempdir, top_file)
            shutil.copy(os.path.join(HERE, 'gromacs', top_file), tmp_top_file)
            load_dir = os.path.join(tempdir, 'loader')
            shutil.copytree(GMX_LOAD, load_dir)
            settings['engine']['input_path'] = gmx_dir
            settings['initial-path']['load_folder'] = load_dir
            units_from_settings(settings)

            # Reset inputs for this example set up in the tempdir
            del settings['ensemble']
            fill_up_tis_and_retis_settings(settings)

            with patch('sys.stdout', new=StringIO()):
                simulation = create_simulation(settings)
                simulation.set_up_output(settings)
                init = simulation.initiate(settings)
                assert init
            assert self._test_dirs(tempdir)

            for i in range(3):
                assert simulation.ensembles[i]['path_ensemble'] \
                   .last_path.status == 'ACC'

            assert os.path.exists(os.path.join(load_dir, '004',
                                               'traj.txt'))

            for i in range(1, 3):
                path = simulation.ensembles[i]['path_ensemble'].last_path
                for j in range(int(path.length/3)):
                    assert path.phasepoints[j].order[0] == \
                                     path.phasepoints[-1-j].order[0]

    def test_external_load_external_trr(self):
        """Test the load initialise for external trajectories."""
        settings = parse_settings_file(
            os.path.join(HERE, 'gromacs', 'auto-load.rst')
        )
        orderp_file = os.path.join(HERE, 'gromacs', 'orderp.py')
        settings['engine']['gmx'] = GMX
        settings['engine']['mdrun'] = MDRUN
        del settings['initial-path']['top_file']

        with tempfile.TemporaryDirectory() as tempdirbase:
            tempdir = os.path.join(tempdirbase, 'gromacs')
            settings['simulation']['exe_path'] = tempdir
            gmx_dir = os.path.join(tempdir, 'gmx_input')
            shutil.copytree(GMX_DIR, gmx_dir)
            shutil.copy(orderp_file, tempdir)
            load_dir = os.path.join(tempdir, 'load-trj')
            shutil.copytree(GMX_LOAD_TRJ, load_dir)
            settings['engine']['input_path'] = gmx_dir
            settings['engine']['exe_path'] = tempdir

            units_from_settings(settings)
            settings['initial-path']['load_folder'] = load_dir

            # Reset inputs for this example set up in the tempdir
            del settings['ensemble']
            fill_up_tis_and_retis_settings(settings)

            with patch('sys.stdout', new=StringIO()):
                simulation = create_simulation(settings)
                simulation.set_up_output(settings)
                init = simulation.initiate(settings)
                assert init
            assert self._test_dirs(tempdir)

            for i in range(4):
                assert simulation.ensembles[i]['path_ensemble'] \
                   .last_path.status == 'ACC'

            assert os.path.exists(os.path.join(load_dir, '004',
                                               'traj.txt'))

            exp_len = [None, 11, 11, 11]
            for i in range(1, 4):
                path = simulation.ensembles[i]['path_ensemble'].last_path
                assert path.length == exp_len[i]

    @staticmethod
    def _create_empty_gro_files(tempdir, nfiles):
        """Create a given number of empty .gro files."""
        grofiles = []
        for i in range(nfiles):
            gro = f'conf{i:04d}.gro'
            grofile = os.path.join(tempdir, gro)
            pathlib.Path(grofile).touch()
            grofiles.append(gro)
        return grofiles

    def test_load_job(self):
        """Test some aspects of the dirty load method."""
        settings = parse_settings_file(os.path.join(HERE, 'internal', 'tis',
                                                    'tis.rst'))

        with tempfile.TemporaryDirectory() as tempdirbase:
            tempdir = os.path.join(tempdirbase, 'gromacs')
            self._set_up_for_load(settings, tempdir)
            # Remove the ensemble directory to test what happens:
            edir = os.path.join(tempdir, '001')
            assert os.path.isdir(tempdirbase)
            assert os.path.isdir(tempdir)
            assert os.path.isdir(edir)
            shutil.rmtree(edir)
            assert not os.path.isdir(edir)
            # And remove the initial.xyz:
            remove_file(os.path.join(tempdir, 'initial.xyz'))
            with pytest.raises(FileNotFoundError):
                _do_the_dirty_load_job(tempdir, edir)
            # Create some fake .gro files to see what happens next:
            grofiles = self._create_empty_gro_files(tempdir, 3)
            accepted = os.path.join(edir, 'accepted')
            assert not os.path.isdir(edir)
            assert not os.path.isdir(accepted)
            _do_the_dirty_load_job(tempdir, edir)
            # Check that edir was created:
            assert os.path.isdir(edir)
            # Check that the accepted dir was crated:
            assert os.path.isdir(accepted)
            # And that the gro files were copied:
            for i in grofiles:
                assert os.path.isfile(os.path.join(accepted, i))
                remove_file(os.path.join(tempdir, i))
            # Remove the accepted dir and see what happens:
            shutil.rmtree(accepted)
            with pytest.raises(FileNotFoundError):
                _do_the_dirty_load_job(tempdir, edir)
            # Add some files in tempdir and see what happens:
            grofiles = self._create_empty_gro_files(tempdir, 3)
            assert not os.path.isdir(accepted)
            _do_the_dirty_load_job(tempdir, edir)
            assert os.path.isdir(accepted)
            # Check that files were copied:
            for i in grofiles:
                assert os.path.isfile(os.path.join(accepted, i))
                assert os.path.isfile(os.path.join(tempdir, i))
                remove_file(os.path.join(tempdir, i))
            shutil.rmtree(accepted)
            # Add some files in edir and see what happens:
            grofiles = self._create_empty_gro_files(edir, 3)
            assert not os.path.isdir(accepted)
            _do_the_dirty_load_job(tempdir, edir)
            assert os.path.isdir(accepted)
            # Check that files were moved:
            for i in grofiles:
                assert os.path.isfile(os.path.join(accepted, i))
                assert not os.path.isfile(os.path.join(tempdir, i))


class TestLoadExternalTrajectoryNoDuplicatePath:
    """Regression test: _load_external_trajectory must not double dirname.

    When traj.txt contains filenames with a path prefix (e.g. 'load/frame.gro'
    instead of just 'frame.gro'), _load_external_trajectory must strip the
    prefix before joining with dirname, so we get 'load/frame.gro' — not
    'load/load/frame.gro'.
    """

    GRO_CONTENT = (
        'MD of 2 waters, t= 0.0\n'
        '    6\n'
        '    1WATER  OW     1   0.126   1.624   1.679'
        '  0.1227 -0.0580  0.0434\n'
        '    1WATER  HW1    2   0.190   1.661   1.747'
        '  0.8085  0.3191 -0.7791\n'
        '    1WATER  HW2    3   0.177   1.568   1.613'
        ' -0.9045 -2.6469  1.3180\n'
        '    2WATER  OW     4   2.275   0.053   0.622'
        '  0.2519  0.3140 -0.1734\n'
        '    2WATER  HW1    5   2.337   0.002   0.680'
        ' -1.0641 -1.1349  0.0257\n'
        '    2WATER  HW2    6   2.326   0.120   0.568'
        '  1.9427 -0.8216 -0.0244\n'
        '   9.82060   9.82060   9.82060\n'
    )

    @staticmethod
    def _write_traj_txt(filepath, filenames):
        """Write a traj.txt with the given filenames."""
        with open(filepath, 'w', encoding='utf-8') as fh:
            fh.write('# Cycle: 0, status: ACC\n')
            fh.write('#     Step              Filename'
                     '       index    vel\n')
            for i, name in enumerate(filenames):
                fh.write(f'{i:>10}  {name:>20s}  {i:>10}  {1:>5}\n')

    def test_no_double_prefix_copy_false(self):
        """Paths must not be doubled when copy=False."""
        with tempfile.TemporaryDirectory() as tempdir:
            load_dir = os.path.join(tempdir, 'load')
            accepted_dir = os.path.join(load_dir, 'accepted')
            os.makedirs(accepted_dir)

            # Create frame files in accepted/
            frames = ['frame_001.gro', 'frame_002.gro']
            for frame in frames:
                with open(os.path.join(accepted_dir, frame),
                          'w', encoding='utf-8') as fh:
                    fh.write(self.GRO_CONTENT)

            # Write traj.txt with path-prefixed filenames (the bug trigger)
            prefixed = [os.path.join('load', f) for f in frames]
            self._write_traj_txt(os.path.join(load_dir, 'traj.txt'),
                                 prefixed)

            engine = type('E', (), {'exe_dir': None})()
            with patch('sys.stdout', new=StringIO()):
                traj = _load_external_trajectory(load_dir, engine,
                                                 copy=False)

            for snapshot in traj['data']:
                path = snapshot[1]
                # The path must be load/frame_XXX.gro, never
                # load/load/frame_XXX.gro
                assert path.count('load') == 1, \
                    f'Doubled dirname in path: {path}'
                assert os.path.basename(path) in frames

    def test_no_double_prefix_copy_true(self):
        """Paths must not be doubled when copy=True."""
        with tempfile.TemporaryDirectory() as tempdir:
            load_dir = os.path.join(tempdir, 'load')
            accepted_dir = os.path.join(load_dir, 'accepted')
            exe_dir = os.path.join(tempdir, 'exe')
            os.makedirs(accepted_dir)
            os.makedirs(exe_dir)

            frames = ['frame_001.gro', 'frame_002.gro']
            for frame in frames:
                with open(os.path.join(accepted_dir, frame),
                          'w', encoding='utf-8') as fh:
                    fh.write(self.GRO_CONTENT)

            # Write traj.txt with path-prefixed filenames (the bug trigger)
            prefixed = [os.path.join('load', f) for f in frames]
            self._write_traj_txt(os.path.join(load_dir, 'traj.txt'),
                                 prefixed)

            engine = type('E', (), {'exe_dir': exe_dir})()
            with patch('sys.stdout', new=StringIO()):
                traj = _load_external_trajectory(load_dir, engine,
                                                 copy=True)

            for snapshot in traj['data']:
                path = snapshot[1]
                # Filenames must point into exe_dir with no extra prefix
                assert path.startswith(exe_dir), \
                    f'Path does not start with exe_dir: {path}'
                assert os.path.basename(path) in frames
            # Files must have been copied to exe_dir
            for frame in frames:
                assert os.path.isfile(os.path.join(exe_dir, frame))

    def test_basenames_still_work(self):
        """Normal traj.txt with bare basenames must keep working."""
        with tempfile.TemporaryDirectory() as tempdir:
            load_dir = os.path.join(tempdir, 'load')
            accepted_dir = os.path.join(load_dir, 'accepted')
            os.makedirs(accepted_dir)

            frames = ['frame_001.gro', 'frame_002.gro']
            for frame in frames:
                with open(os.path.join(accepted_dir, frame),
                          'w', encoding='utf-8') as fh:
                    fh.write(self.GRO_CONTENT)

            # Write traj.txt with bare basenames (normal case)
            self._write_traj_txt(os.path.join(load_dir, 'traj.txt'),
                                 frames)

            engine = type('E', (), {'exe_dir': None})()
            with patch('sys.stdout', new=StringIO()):
                traj = _load_external_trajectory(load_dir, engine,
                                                 copy=False)

            for snapshot in traj['data']:
                path = snapshot[1]
                expected = os.path.join(load_dir, os.path.basename(path))
                assert path == expected, \
                    f'Unexpected path: {path} != {expected}'
