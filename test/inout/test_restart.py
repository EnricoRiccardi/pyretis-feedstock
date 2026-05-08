# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the restart methods."""
from io import StringIO
import logging
import os
import pytest
import tempfile
import numpy as np
from unittest.mock import patch
from pyretis.simulation.simulation import Simulation
from pyretis.core.box import create_box
from pyretis.core.system import System
from pyretis.core.common import big_fat_comparer
from pyretis.core.random_gen import create_random_generator
from pyretis.core.particles import Particles
from pyretis.forcefield.forcefield import ForceField
from pyretis.forcefield.potentials import PairLennardJonesCutnp
from pyretis.setup.createsystem import create_system
from pyretis.core.units import (units_from_settings,
                                create_conversion_factors)
from pyretis.engines.internal import VelocityVerlet
from pyretis.inout.settings import (add_default_settings,
                                    fill_up_tis_and_retis_settings,
                                    SECTIONS)
from pyretis.inout.restart import (write_restart_file,
                                   write_ensemble_restart,
                                   read_restart_file)
from pyretis.setup.createsimulation import create_retis_simulation
from pyretis.inout.common import make_dirs
from pyretis.tools.lattice import generate_lattice


logging.disable(logging.CRITICAL)

HERE = os.path.abspath(os.path.dirname(__file__))


def create_test_system():
    """Create a system we can use for testing."""
    create_conversion_factors('lj')
    xyz, size = generate_lattice('fcc', [2, 2, 2], density=0.9)
    size = np.array(size)
    box = create_box(low=size[:, 0], high=size[:, 1])
    system = System(units='lj', box=box, temperature=2.0)
    system.particles = Particles(dim=3)
    for pos in xyz:
        system.add_particle(pos, vel=np.zeros_like(pos),
                            force=np.zeros_like(pos),
                            mass=1.0, name='Ar', ptype=0)
    rgen = create_random_generator({'seed': 0})
    system.generate_velocities(rgen, distribution='maxwell', momentum=True)
    potentials = [
        PairLennardJonesCutnp(dim=3, shift=True, mixing='geometric'),
    ]
    parameters = [
        {0: {'sigma': 1, 'epsilon': 1, 'rcut': 2.5}},
    ]
    system.forcefield = ForceField(
        'Lennard Jones force field',
        potential=potentials,
        params=parameters,
    )

    return system


class TestRestartMethods:
    """Test methods defined in the module."""

    def test_write_and_read(self):
        """Test write/read for simulation restart files."""
        settings = {
            'simulation': {},
            'system': {
                'dimensions': 3,
                'units': 'reduced',
                'temperature': 1.0
            },
            'particles': {
                'position': {
                    'generate': 'fcc',
                    'repeat': [2, 2, 2],
                    'density': 0.9,
                }
            },
            'potential': [
                {
                    'class': 'PairLennardJonesCutnp',
                    'shift': True,
                    'parameter': {
                        0: {'sigma': 1.0, 'epsilon': 1.0, 'rcut': 2.5}
                    }
                }
            ]
        }
        simulation = Simulation(settings, {'steps': 100})
        add_default_settings(settings)
        units_from_settings(settings)
        system = create_system(settings)
        simulation.system = system

        with tempfile.NamedTemporaryFile() as tmp:
            write_restart_file(tmp.name, simulation)
            tmp.flush()
            read = read_restart_file(tmp.name)
            for key in ['simulation', 'system']:
                assert key in read
            assert big_fat_comparer(read, simulation.restart_info())

    def test_read_write_ensemble(self):
        """Test read/write for path ensemble restart files."""
        settings = {
            'simulation': {'task': 'retis',
                           'interfaces': [-1., 0., 1.],
                           'exe_path': '.',
                           'steps': 10,
                           'zero_ensemble': True,
                           'flux': True},
            'tis': SECTIONS['tis'],
            'retis': SECTIONS['retis'],
            'engine': {'obj': VelocityVerlet(0.002)},
            'particles': {'type': 'internal'},
            'system': {'type': 'internal', 'units': 'lj',
                       'obj': create_test_system()}}

        settings['orderparameter'] = {'class': 'Position',
                                      'dim': 'x', 'index': 0,
                                      'periodic': False}
        with tempfile.TemporaryDirectory() as tempdir:
            settings['simulation']['exe_path'] = tempdir
            fill_up_tis_and_retis_settings(settings)

            with patch('sys.stdout', new=StringIO()):
                simulation = create_retis_simulation(settings)

            filename = os.path.join(tempdir, 'devil.rst')

            write_restart_file(filename, simulation)
            assert os.path.exists(filename)

            info = simulation.restart_info()
            info_file = read_restart_file(filename)

            big_fat_comparer(info, info_file, hard=True)

            ensembles = simulation.ensembles
            for ens_set, ensemble in zip(settings['ensemble'], ensembles):
                for name in ensemble['path_ensemble'].directories():
                    make_dirs(os.path.join(tempdir, name))
                write_ensemble_restart(ensemble, ens_set)
                restart_file = os.path.join(
                    tempdir,
                    ensemble['path_ensemble'].directory['path_ensemble'],
                    'ensemble.restart')
                read = read_restart_file(restart_file)
                restart = {**ens_set,
                           **{'system': ensemble['system'].restart_info(),
                              'engine': ensemble['engine'].restart_info(),
                              'path_ensemble':
                                  ensemble['path_ensemble'].restart_info(),
                              'rgen': ensemble['rgen'].get_state(),
                              'order_function':
                              ensemble['order_function'].restart_info()
                              }
                           }

                assert big_fat_comparer(read, restart, hard=True)
