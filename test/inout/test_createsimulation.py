# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the methods in pyretis.setup.createsimulation"""
import os
import logging
import pytest
import numpy as np
from io import StringIO
from unittest.mock import patch
from pyretis.engines.internal import VelocityVerlet
from pyretis.core.units import create_conversion_factors
from pyretis.core.system import System
from pyretis.core.particles import Particles
from pyretis.core.pathensemble import (
    PathEnsemble,
    PathEnsembleExt,
)
from pyretis.setup.createsimulation import (
    create_ensemble,
    create_ensembles,
    create_nve_simulation,
    create_mdflux_simulation,
    create_md_simulation,
    create_umbrellaw_simulation,
    create_tis_simulation,
    create_retis_simulation,
    create_simulation,
)
from pyretis.tools.lattice import generate_lattice
from pyretis.forcefield.forcefield import ForceField
from pyretis.forcefield.potentials import PairLennardJonesCutnp
from pyretis.core.box import create_box
from pyretis.inout.settings import SECTIONS

from pyretis.simulation.simulation import Simulation
from pyretis.simulation.md_simulation import (
    SimulationNVE,
    SimulationMD,
    SimulationMDFlux,
)
from pyretis.simulation.mc_simulation import UmbrellaWindowSimulation
from pyretis.simulation.path_simulation import (
    SimulationTIS,
    SimulationRETIS
)
from pyretis.inout.settings import fill_up_tis_and_retis_settings
from pyretis.inout.checker import check_ensemble

logging.disable(logging.CRITICAL)


HERE = os.path.abspath(os.path.dirname(__file__))

ORDER_ERROR_MSG = (
    'No order parameter was defined,'
    ' but the engine *does* require it.'
)
MISSING_STEPS_ERROR_MSG = 'Simulation setting "steps" is missing!'
MISSING_INTERFACES_ERROR_MSG = 'Simulation setting "interfaces" is missing!'


def create_test_settings():
    """Create a settings we can use for testing."""
    settings = {'simulation': {'task': 'retis',
                               'interfaces': [-1., 0, 1],
                               'steps': 10,
                               'exe_path': HERE,
                               'zero_ensemble': True,
                               'flux': True},
                'orderparameter': {'class': 'Position',
                                   'dim': 'x',
                                   'index': 0,
                                   'periodic': False},
                'particles': {'type': 'internal'},
                'system': {'units': 'lj',
                           'obj': create_test_system()},
                'engine': {'obj': VelocityVerlet(0.002)},
                'tis': {'shooting_moves': ['sh', 'sh', 'sh']}}
    return settings


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
    gen_settings = {'distribution': 'maxwell', 'momentum': True, 'seed': 0}
    system.generate_velocities(**gen_settings)
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


class TestMethods:
    """Test some of the methods from .createsimulation."""
    def test_create_path_ensemble(self):
        """Test create_path_ensemble."""
        path_ensemble = PathEnsemble(ensemble_number=2, interfaces=[-1., 0, 1])
        assert path_ensemble.ensemble_number == 2
        assert isinstance(path_ensemble, PathEnsemble)
        path_ensemble = PathEnsembleExt(ensemble_number=2,
                                        interfaces=[-1., 0, 1])
        assert isinstance(path_ensemble, PathEnsembleExt)
        assert path_ensemble.ensemble_number == 2
        # Test with some "missing" settings:
        settings = {'simulation': {'interfaces': [-1., 0.]}}
        with patch('sys.stdout', new=StringIO()):
            with pytest.raises(KeyError):
                create_ensemble(settings)

    def test_create_ensemble(self):
        """Test create_ensemble."""
        settings = create_test_settings()

        fill_up_tis_and_retis_settings(settings)
        check_ensemble(settings)
        with patch('sys.stdout', new=StringIO()):
            ensembles = create_ensembles(settings)
        for i_ens, ensemble in enumerate(ensembles):
            path_ensemble = ensemble['path_ensemble']
            assert path_ensemble.ensemble_number == i_ens
            assert isinstance(path_ensemble, PathEnsemble)

        for i_ens in range(len(ensembles)):
            settings['ensemble'][i_ens]['particles']['type'] = 'external'

        with patch('sys.stdout', new=StringIO()):
            ensembles = create_ensembles(settings)
        for i_ens, ensemble in enumerate(ensembles):
            path_ensemble = ensemble['path_ensemble']
            assert isinstance(path_ensemble, PathEnsembleExt)
            assert path_ensemble.ensemble_number == i_ens

    def test_create_ensemble_permeability(self):
        """Test create_ensemble."""
        settings = create_test_settings()
        settings['simulation'] = {'task': 'retis',
                                  'interfaces': [-1., 0, 1],
                                  'exe_path': HERE,
                                  'zero_ensemble': True,
                                  'flux': True,
                                  'permeability': True,
                                  'zero_left': -2}

        fill_up_tis_and_retis_settings(settings)
        with patch('sys.stdout', new=StringIO()):
            ensembles = create_ensembles(settings)
        ens = ensembles[0]
        assert ens['path_ensemble'].interfaces == (-2, -1.5, -1)
        assert ["R", "L"] == ens['path_ensemble'].start_condition

    def test_create_ensemble_pptis(self):
        """Test create_ensemble for pptis."""
        # with flux=False, zero_ensemble=True
        # -------------------------------------
        settings = create_test_settings()
        settings['tis']['shooting_moves'] += ['sh', 'sh']
        settings['simulation'] = {'task': 'pptis',
                                  'interfaces': [-1., 0., 1., 2., 3.],
                                  'flux': False,
                                  'zero_ensemble': True,
                                  'exe_path': HERE}

        fill_up_tis_and_retis_settings(settings)
        with patch('sys.stdout', new=StringIO()):
            ensembles = create_ensembles(settings)
        n_ens = len(ensembles)
        for i in range(n_ens):
            assert ensembles[i]['path_ensemble'].start_condition == ['R', 'L']
        hoped = [(-1., -1., 0.), (-1., 0., 1.),
                 (0., 1., 2.), (1., 2., 3.)]
        assert len(hoped) == n_ens
        for ensemble, hope in zip(ensembles, hoped):
            assert ensemble['path_ensemble'].interfaces == hope

    def test_create_ensemble_pptis_1(self):
        """Test create_ensemble for pptis."""
        # with flux=True, zero_ensemble=True
        # -------------------------------------
        settings = create_test_settings()
        settings['tis']['shooting_moves'] += ['sh', 'sh']
        settings['simulation'] = {'task': 'pptis',
                                  'interfaces': [-1., 0., 1., 2., 3.],
                                  'flux': True,
                                  'zero_ensemble': True,
                                  'exe_path': HERE}

        fill_up_tis_and_retis_settings(settings)
        with patch('sys.stdout', new=StringIO()):
            ensembles = create_ensembles(settings)
        n_ens = len(ensembles)
        # start condition
        assert ensembles[0]['path_ensemble'].start_condition == 'R'
        for i in range(1, n_ens):
            assert ensembles[i]['path_ensemble'].start_condition == ['R', 'L']
        # interfaces
        hoped = [(-np.inf, -1., -1.), (-1., -1., 0.), (-1., 0., 1.),
                 (0., 1., 2.), (1., 2., 3.)]
        assert len(hoped) == n_ens
        for ensemble, hope in zip(ensembles, hoped):
            assert ensemble['path_ensemble'].interfaces == hope

    def test_create_ensemble_pptis_2(self):
        """Test create_ensemble for repptis."""
        # with flux=False, zero_ensemble=False
        # -------------------------------------
        settings = create_test_settings()
        settings['tis']['shooting_moves'] += ['sh', 'sh']
        settings['simulation'] = {'task': 'pptis',
                                  'interfaces': [-1., 0., 1., 2., 3.],
                                  'flux': False,
                                  'zero_ensemble': False,
                                  'exe_path': HERE}
        # just [1+], [2+], [3+]
        fill_up_tis_and_retis_settings(settings)
        with patch('sys.stdout', new=StringIO()):
            ensembles = create_ensembles(settings)
        n_ens = len(ensembles)
        # start condition
        for i in range(n_ens):
            assert ensembles[i]['path_ensemble'].start_condition == ['R', 'L']
        # interfaces
        hoped = [(-1., 0., 1.), (0., 1., 2.), (1., 2., 3.)]
        assert len(hoped) == n_ens
        for ensemble, hope in zip(ensembles, hoped):
            assert ensemble['path_ensemble'].interfaces == hope

    def test_create_ensemble_pptis_3(self):
        """Test create_ensemble for pptis."""
        # with flux=True, zero_ensemble=False
        # ------------------------------------
        settings = create_test_settings()
        settings['tis']['shooting_moves'] += ['sh', 'sh']
        settings['simulation'] = {'task': 'pptis',
                                  'interfaces': [-1., 0., 1., 2., 3.],
                                  'flux': True,
                                  'zero_ensemble': False,
                                  'exe_path': HERE}
        # [0-], [1+], [2+], [3+]
        fill_up_tis_and_retis_settings(settings)
        with patch('sys.stdout', new=StringIO()):
            ensembles = create_ensembles(settings)
        n_ens = len(ensembles)
        # start condition
        assert ensembles[0]['path_ensemble'].start_condition == 'R'
        for i in range(1, n_ens):
            assert ensembles[i]['path_ensemble'].start_condition == ['R', 'L']
        # interfaces
        hoped = [(-np.inf, -1., -1.), (-1., 0., 1.),
                 (0., 1., 2.), (1., 2., 3.)]
        assert len(hoped) == n_ens
        for ensemble, hope in zip(ensembles, hoped):
            assert ensemble['path_ensemble'].interfaces == hope

    def test_create_ensemble_repptis_0(self):
        """Test create_ensemble for repptis."""
        # with flux=False, zero_ensemble=False, no permeability
        settings = create_test_settings()
        settings['simulation'] = {'task': 'repptis',
                                  'interfaces': [-1., 0., 1., 2., 3.],
                                  'flux': False,
                                  'zero_ensemble': False,
                                  'exe_path': HERE}
        fill_up_tis_and_retis_settings(settings)
        with patch('sys.stdout', new=StringIO()):
            ensembles = create_ensembles(settings)
        n_ens = len(ensembles)
        for i in range(n_ens):
            assert ensembles[i]['path_ensemble'].start_condition == ['R', 'L']
        hoped = [(-1., 0., 1.), (0., 1., 2.), (1., 2., 3.)]
        assert len(hoped) == n_ens
        for ensemble, hope in zip(ensembles, hoped):
            assert ensemble['path_ensemble'].interfaces == hope

    def test_create_ensemble_repptis_1(self):
        """Test create_ensemble for repptis."""
        # with flux=False, zero_ensemble=False, with permeability
        settings = create_test_settings()
        settings['simulation'] = {'task': 'repptis',
                                  'interfaces': [-1., 0., 1., 2., 3.],
                                  'exe_path': HERE,
                                  'zero_ensemble': True,
                                  'flux': True,
                                  'zero_left': -2.,
                                  'permeability': True}
        # [0-], [0+], [1+], [2+], [3+]
        settings['tis']['shooting_moves'] += ['sh', 'sh']
        fill_up_tis_and_retis_settings(settings)
        with patch('sys.stdout', new=StringIO()):
            ensembles = create_ensembles(settings)
        # start_condition for paths
        n_ens = len(ensembles)
        assert ensembles[0]['path_ensemble'].start_condition == ['R', 'L']
        for i in range(1, n_ens):
            assert ensembles[i]['path_ensemble'].start_condition == ['R', 'L']
        hoped = [(-2, -1.5, -1), (-1., -1., 0.), (-1., 0., 1.),
                 (0., 1., 2.), (1., 2., 3.)]
        assert len(hoped) == n_ens
        for ensemble, hope in zip(ensembles, hoped):
            assert ensemble['path_ensemble'].interfaces == hope

    def test_create_ensemble_repptis_2(self):
        """Test create_ensemble for repptis."""
        settings = create_test_settings()
        settings['simulation'] = {'task': 'repptis',
                                  'interfaces': [-1., 0., 1., 2., 3.],
                                  'exe_path': HERE,
                                  'flux': False,
                                  'zero_ensemble': True,
                                  }
        settings['tis']['shooting_moves'] += ['sh']
        fill_up_tis_and_retis_settings(settings)
        with patch('sys.stdout', new=StringIO()):
            ensembles = create_ensembles(settings)
        assert ensembles[-1]['path_ensemble'].start_condition == ['R', 'L']
        hoped = [(-1., -1., 0.), (-1., 0., 1.),
                 (0., 1., 2.), (1., 2., 3.)]
        assert len(hoped) == len(ensembles)
        for ensemble, hope in zip(ensembles, hoped):
            assert ensemble['path_ensemble'].interfaces == hope

    def test_raise_on_wrong_mirror(self):
        """Test bad mirror."""
        settings = create_test_settings()
        settings['simulation'] = {'task': 'retis',
                                  'interfaces': [-1., 0, 1],
                                  'exe_path': HERE,
                                  'zero_ensemble': True,
                                  'flux': True,
                                  'zero_left': -2,
                                  'permeability': True}
        settings['orderparameter'] = {'class': 'Permeability',
                                      'dim': 'x',
                                      'index': 0,
                                      'mirror_pos': 0}

        fill_up_tis_and_retis_settings(settings)
        with pytest.raises(ValueError, match="have a mirror at -1.5"):
            with patch('sys.stdout', new=StringIO()):
                create_ensembles(settings)

    def test_create_ensembles(self):
        """Test create_ensembles."""
        settings = create_test_settings()
        settings['simulation'] = {'task': 'tis',
                                  'interfaces': [-1., 0., 1.0, 2.0],
                                  'exe_path': HERE,
                                  'zero_ensemble': False,
                                  'flux': False}
        fill_up_tis_and_retis_settings(settings)
        check_ensemble(settings)
        with patch('sys.stdout', new=StringIO()):
            ensembles = create_ensembles(settings)
        detects = [1.0, 2.0]
        assert len(ensembles) == 2
        for iens, ens in enumerate(ensembles):
            assert isinstance(ens['path_ensemble'], PathEnsemble)
        for iens, ens in enumerate(settings['ensemble']):
            assert detects[iens] == ens['tis']['detect']

    def test_fail_simulation(self):
        """Test fail simulation."""
        settings = {'simulation': {'task': 'fake'}}
        with pytest.raises(ValueError) as err:
            create_simulation(settings)
        assert 'Unknown simulation task fake' in str(err.value)

    def test_load_simulation(self):
        """Test cycle counter."""
        settings = {}
        sim = Simulation(settings, controls={'steps': 16})
        restart = {'simulation': {'cycle': {'steps': 30,
                                            'endcycle': 30}}}
        sim.load_restart_info(restart)
        assert sim.cycle['steps'] == 16
        assert sim.cycle['step'] == 0

    def test_create_nve_simulation(self):
        """Test create_nve_simulation."""
        settings = {'simulation': {'task': 'md-nve', 'steps': 10},
                    'system': {'obj': create_test_system()},
                    'engine': {'obj': VelocityVerlet(0.002)}}
        sim = create_nve_simulation(settings)
        assert isinstance(sim, SimulationNVE)
        assert sim.simulation_type == SimulationNVE.simulation_type
        del settings['simulation']['steps']
        with pytest.raises(ValueError) as err:
            create_nve_simulation(settings)
        assert MISSING_STEPS_ERROR_MSG == str(err.value)

    def test_create_md_simulation(self):
        """Test create_nve_simulation."""
        settings = {'simulation': {'task': 'md', 'steps': 10},
                    'system': {'obj': create_test_system()},
                    'engine': {'obj': VelocityVerlet(0.002)}}
        sim = create_md_simulation(settings)
        assert isinstance(sim, SimulationMD)
        assert sim.simulation_type == SimulationMD.simulation_type
        del settings['simulation']['steps']
        with pytest.raises(ValueError) as err:
            create_md_simulation(settings)
        assert MISSING_STEPS_ERROR_MSG == str(err.value)

    def test_create_mdflux_simulation(self):
        """Test create_mdflux_simulation."""
        settings = {'simulation': {'task': 'mdflux'},
                    'system': {'obj': create_test_system()},
                    'engine': {'obj': VelocityVerlet(0.002)}}

        with pytest.raises(ValueError) as err:
            create_mdflux_simulation(settings)
        assert ORDER_ERROR_MSG == str(err.value)
        settings['orderparameter'] = {'class': 'Position',
                                      'dim': 'x', 'index': 0,
                                      'periodic': False}
        # Test that we fail because required setting "steps" are missing:
        with pytest.raises(ValueError) as err:
            create_mdflux_simulation(settings)
        assert MISSING_STEPS_ERROR_MSG == str(err.value)
        settings['simulation']['steps'] = 10
        # Test that we fail because required setting "interfaces" are
        # missing:
        with pytest.raises(ValueError) as err:
            create_mdflux_simulation(settings)
        assert MISSING_INTERFACES_ERROR_MSG in str(err.value)
        settings['simulation']['interfaces'] = [0, 1]
        sim = create_mdflux_simulation(settings)
        assert isinstance(sim, SimulationMDFlux)

    def test_create_umbrellawsimulation(self):
        """Test create_umbrellaw_simulation."""
        settings = {'simulation': {'task': 'umbrellawindow', 'seed': 3},
                    'system': {'obj': create_test_system}}

        with pytest.raises(ValueError) as err:
            create_umbrellaw_simulation(settings)
        assert 'Simulation setting "umbrella" is missing!' == str(err.value)

        settings['simulation']['umbrella'] = [-1.0, 0.0]
        with pytest.raises(ValueError) as err:
            create_umbrellaw_simulation(settings)
        assert 'Simulation setting "overlap" is missing!' == str(err.value)

        settings['simulation']['overlap'] = -0.5
        with pytest.raises(ValueError) as err:
            create_umbrellaw_simulation(settings)
        assert 'Simulation setting "maxdx" is missing!' in str(err.value)

        settings['simulation']['maxdx'] = 1.0
        with pytest.raises(ValueError) as err:
            create_umbrellaw_simulation(settings)
        assert 'Simulation setting "mincycle" is missing!' == str(err.value)
        settings['simulation']['mincycle'] = 10
        sim = create_umbrellaw_simulation(settings)
        assert isinstance(sim, UmbrellaWindowSimulation)

        del settings['system']
        with patch('sys.stdout', new=StringIO()):
            with pytest.raises(KeyError) as err:
                create_umbrellaw_simulation(settings)
        assert "'system'" == str(err.value)

    def test_create_tis_simulation(self):
        """Test create_tis_simulation."""
        settings = {'simulation': {'task': 'tis',
                                   'interfaces': [-1., 0., 1.],
                                   'exe_path': HERE,
                                   'zero_ensemble': False,
                                   'flux': False},
                    'tis': SECTIONS['tis'],
                    'particles': {'type': 'internal'},
                    'system': {'obj': create_test_system()},
                    'engine': {'obj': VelocityVerlet(0.002)},
                    'output': SECTIONS['output']}
        fill_up_tis_and_retis_settings(settings)

        # Test that the set up fails, when we did not define an
        # order parameter and the engine wants it:
        with patch('sys.stdout', new=StringIO()):
            with pytest.raises(ValueError) as err:
                create_tis_simulation(settings)
        assert ORDER_ERROR_MSG == str(err.value)
        # Continue testing with an order parameter defined:
        settings['orderparameter'] = {'class': 'Position',
                                      'dim': 'x', 'index': 0,
                                      'periodic': False}
        # Test that we fail when we are missing some settings:
        fill_up_tis_and_retis_settings(settings)
        with patch('sys.stdout', new=StringIO()):
            with pytest.raises(ValueError) as err:
                create_tis_simulation(settings)
        # it should raise also missing steps
        assert MISSING_STEPS_ERROR_MSG == str(err.value)
        settings['simulation']['steps'] = 10
        fill_up_tis_and_retis_settings(settings)
        with patch('sys.stdout', new=StringIO()):
            sim = create_tis_simulation(settings)
        assert isinstance(sim, SimulationTIS)
        del settings['ensemble']
        del settings['system']
        fill_up_tis_and_retis_settings(settings)
        with patch('sys.stdout', new=StringIO()):
            with pytest.raises(KeyError) as err:
                create_tis_simulation(settings)
        assert "'system'" == str(err.value)

    def test_create_tis_simulations(self):
        """Test create_tis_simulations."""
        settings = create_test_settings()
        settings['simulation'] = {'task': 'tis',
                                  'interfaces': [-1., 0., 1., 2],
                                  'steps': 10,
                                  'exe_path': HERE,
                                  'zero_ensemble': False,
                                  'flux': False}
        settings['tis'] = SECTIONS['tis']
        fill_up_tis_and_retis_settings(settings)
        with patch('sys.stdout', new=StringIO()):
            simulations = create_tis_simulation(settings)
        assert isinstance(simulations, SimulationTIS)
        del settings['ensemble']
        del settings['engine']
        fill_up_tis_and_retis_settings(settings)
        with patch('sys.stdout', new=StringIO()):
            with pytest.raises(ValueError) as err:
                create_tis_simulation(settings)
        assert "Could not create engine from settings!" == str(err.value)

    def test_create_retis_simulation(self):
        """Test create_retis_simulation."""
        settings = {'simulation': {'task': 'retis',
                                   'interfaces': [-1., 0., 1.],
                                   'exe_path': '.',
                                   'zero_ensemble': True,
                                   'flux': True},
                    'tis': SECTIONS['tis'],
                    'retis': SECTIONS['retis'],
                    'engine': {'obj': VelocityVerlet(0.002)},
                    'particles': {'type': 'internal'},
                    'system': {'type': 'internal',
                               'units': 'lj',
                               'obj': create_test_system()}}

        del settings['tis']['ensemble_number']
        fill_up_tis_and_retis_settings(settings)

        # Test that we fail without an order parameter defined:
        with patch('sys.stdout', new=StringIO()):
            with pytest.raises(ValueError) as err:
                create_retis_simulation(settings)
        assert ORDER_ERROR_MSG == str(err.value)
        settings['orderparameter'] = {'class': 'Position',
                                      'dim': 'x', 'index': 0,
                                      'periodic': False}

        fill_up_tis_and_retis_settings(settings)
        # Test that we fail when we are missing some settings:
        with patch('sys.stdout', new=StringIO()):
            with pytest.raises(ValueError) as err:
                create_retis_simulation(settings)
        assert MISSING_STEPS_ERROR_MSG == str(err.value)
        settings['simulation']['steps'] = 10
        with patch('sys.stdout', new=StringIO()):
            sim = create_retis_simulation(settings)
        assert isinstance(sim, SimulationRETIS)
        del settings['ensemble']
        del settings['engine']
        fill_up_tis_and_retis_settings(settings)
        with patch('sys.stdout', new=StringIO()):
            with pytest.raises(ValueError) as err:
                create_tis_simulation(settings)
        assert "Could not create engine from settings!" == str(err.value)

    def test_create_retis_simulation2(self):
        """Test create_simulation for SimulationRETIS."""
        settings = create_test_settings()
        settings['tis'] = SECTIONS['tis']
        settings['retis'] = SECTIONS['retis']

        fill_up_tis_and_retis_settings(settings)
        with patch('sys.stdout', new=StringIO()):
            sim = create_retis_simulation(settings)
        assert isinstance(sim, SimulationRETIS)
        del settings['ensemble']
        del settings['system']
        fill_up_tis_and_retis_settings(settings)
        with patch('sys.stdout', new=StringIO()):
            with pytest.raises(KeyError) as err:
                create_tis_simulation(settings)
        assert "'system'" == str(err.value)

    def test_create_simulation(self):
        """Test create_simulation."""
        settings = {'simulation': {'task': 'does-not-exist'}}
        with pytest.raises(ValueError) as err:
            create_simulation(settings)
        assert 'Unknown simulation task does-not-exist' in str(err.value)

    def test_create_simulationnve(self):
        """Test create_simulation for NVE"""
        settings = {'simulation': {'steps': 10, 'task': 'md-nve'},
                    'system': {'obj': create_test_system()},
                    'engine': {'obj': VelocityVerlet(0.002)}}
        sim = create_simulation(settings)
        assert isinstance(sim, SimulationNVE)

    def test_create_simulationmdflux(self):
        """Test create_simulation for MD-Flux"""
        settings = create_test_settings()
        settings['simulation'] = {'steps': 10,
                                  'task': 'md-flux',
                                  'interfaces': [-1.0]}
        sim = create_mdflux_simulation(settings)
        assert isinstance(sim, SimulationMDFlux)

    def test_create_simulationumb(self):
        """Test create_simulation for UmbrellaWindow"""
        settings = create_test_settings()
        settings['simulation'] = {'task': 'umbrellawindow',
                                  'umbrella': [-1.0, 0.0],
                                  'overlap': -0.5,
                                  'maxdx': 1.0,
                                  'mincycle': 10}
        sim = create_umbrellaw_simulation(settings)
        assert isinstance(sim, UmbrellaWindowSimulation)

    def test_create_simulationtis2(self):
        """Test create_simulation for SimulationTIS multiple."""
        # TIS-multiple:
        settings = create_test_settings()
        settings['simulation'] = {'task': 'tis',
                                  'interfaces': [-1., 0., 1., 2.0],
                                  'steps': 10,
                                  'exe_path': HERE,
                                  'zero_ensemble': False,
                                  'flux': False}
        settings['tis'] = SECTIONS['tis']
        fill_up_tis_and_retis_settings(settings)
        with patch('sys.stdout', new=StringIO()):
            simulations = create_tis_simulation(settings)
        assert isinstance(simulations, SimulationTIS)
