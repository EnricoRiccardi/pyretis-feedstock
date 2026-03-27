# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Define common methods and variables for the tests."""
from contextlib import contextmanager
import logging
import warnings
import os
from unittest.mock import patch
from io import StringIO
import numpy as np
from pyretis.core.box import create_box
from pyretis.core.particles import Particles, ParticlesExt
from pyretis.core.path import Path
from pyretis.core.random_gen import (RandomGenerator,
                                     MockRandomGenerator)
from pyretis.core.system import System
from pyretis.engines.internal import MDEngine
from pyretis.forcefield import ForceField, PotentialFunction
from pyretis.inout.checker import check_ensemble, check_engine
from pyretis.inout.settings import fill_up_tis_and_retis_settings
from pyretis.orderparameter import OrderParameter
from pyretis.setup import create_force_field, create_engine
from pyretis.setup.createsimulation import (create_ensembles,
                                            create_tis_simulation)


HERE = os.path.abspath(os.path.dirname(__file__))


def set_up_system(order, pos, vel, vpot, ekin, internal=True):
    """Create a system for testing."""
    system = System()
    if internal:
        system.particles = Particles(dim=3)
    else:
        system.particles = ParticlesExt(dim=3)
    system.add_particle(pos, vel=vel)
    system.order = order
    system.particles.vpot = vpot
    system.particles.ekin = ekin
    return system


RGEN0 = RandomGenerator(seed=0)
PATHTEST0 = Path(RGEN0)
for _ in range(4):
    for k in range(16):
        PATHTEST0.append(
            set_up_system([k, k], np.zeros(3), np.zeros(3), 0.0, 0.0)
        )

PATHTEST1 = Path(RGEN0)
for k in [4, 5, 6]:
    PATHTEST1.append(
        set_up_system([k, k], np.zeros(3), np.zeros(3), 0.0, 0.0)
    )

PATHTEST2 = Path(RGEN0)
for _ in range(4):
    for k in [4, 5, 6]:
        PATHTEST2.append(
            set_up_system([k, k], np.zeros(3), np.zeros(3), 0.0, 0.0)
        )


@contextmanager
def turn_on_logging():
    """Turn on logging so that tests can detect it."""
    logging.disable(logging.NOTSET)
    try:
        yield
    finally:
        logging.disable(logging.CRITICAL)


class MockEngine(MDEngine):
    """An engine used for testing. It will not do actual dynamics."""

    def __init__(self, timestep, turn_around=20):
        """Just set the time step."""
        super().__init__(timestep, 'Engie McEngineface', dynamics='Fake')
        self.time = 0
        self.total_eclipse = turn_around  # every now and then
        self.delta_v = 0.0123456

    def integration_step(self, system):
        """Do a fake integration step."""
        self.time += 1
        particles = system.particles
        particles.pos += self.timestep * self.delta_v
        system.potential_and_force()
        if self.time > self.total_eclipse:
            self.time = 0
            self.delta_v *= -1.


class MockEngine2(MDEngine):
    """An engine used for testing. It will not do actual dynamics."""

    def __init__(self, timestep, interfaces):
        """Just set the time step."""
        super().__init__(timestep, 'Engie McEngineface', dynamics='Fake')
        self.time = 0
        self.delta_v = 0.0123456
        self.cross_left = False
        self.interfaces = interfaces

    def integration_step(self, system):
        """Do a fake integration step."""
        self.time += 1
        particles = system.particles
        if not self.cross_left:
            particles.pos += self.timestep * self.delta_v
            if particles.pos[0][0] < self.interfaces[0]:
                self.cross_left = True
                self.delta_v *= -1.
        system.potential_and_force()


class MockOrder(OrderParameter):
    """Just return the position as the order parameter."""

    def __init__(self):
        super().__init__(description='Ordey McOrderface')
        # Use this internally, but allow for the change of index
        self._index = 0
        self.index = 0
        self.dim = 0
        self._mirrored = False

    def calculate(self, system):
        if not self._mirrored:
            return [system.particles.pos[self._index][self.dim]]
        return [-system.particles.pos[self._index][self.dim]]

    def mirror(self):
        self._mirrored = not self._mirrored


class MockOrder2(OrderParameter):
    """Just return the position as the order parameter."""

    def __init__(self):
        super().__init__(description='Ordey McOrderface')

    def calculate(self, system):
        return [11.]


class NegativeOrder(OrderParameter):
    """An order parameter which is just the negative of the input order."""

    def __init__(self):
        """Set up the order parameter."""
        super().__init__(description='Mock order parameter', velocity=True)

    def calculate(self, system):
        """Return the negative order parameter."""
        order = system.order
        if order is None:
            return order
        return [i * -1 for i in order]


class SameOrder(OrderParameter):
    """An order parameter which is just the same as the input order."""

    def __init__(self):
        """Set up the order parameter."""
        super().__init__(description='Mock order parameter', velocity=False)

    def calculate(self, system):
        """Return the order parameter."""
        return system.order


class MockPotential(PotentialFunction):
    """Set the potential equal to `x**2`."""

    def __init__(self):
        super().__init__(dim=1, desc='Potey McPotentialface')

    def potential(self, system):
        # pylint: disable=no-self-use,missing-docstring
        pos = system.particles.pos
        vpot = pos**2
        return vpot.sum()

    def force(self, system):
        # pylint: disable=missing-docstring
        pos = system.particles.pos
        forces = pos * 1.0
        virial = np.zeros((self.dim, self.dim))
        return forces, virial

    def potential_and_force(self, system):
        # pylint: disable=missing-docstring
        pot = self.potential(system)
        force, virial = self.force(system)
        return pot, force, virial


def make_internal_system(order=None, pos=None, vel=None,
                         vpot=None, ekin=None):
    """Create a system for testing internal paths."""
    box = create_box(periodic=[False])
    system = System(units='reduced', temperature=1.0, box=box)
    system.particles = Particles(dim=1)
    system.forcefield = ForceField('empty force field',
                                   potential=[MockPotential()])
    system.order = order
    if pos is not None and vel is not None:
        for posi, veli in zip(pos, vel):
            system.add_particle(posi, vel=veli)
    # This is to ensure that the virial is set:
    system.potential_and_force()
    # Update with the given potential and kinetic energy:
    system.particles.vpot = vpot
    system.particles.ekin = ekin
    return system


def make_internal_path(start, end, maxorder, interface=None, points=100):
    """Return a dummy path.

    Parameters
    ----------
    start : tuple of floats
        The starting point for the path.
    end : tuple of floats
        The ending point for the path.
    maxorder : tuple of floats
        The maximum order parameter the path should attain.
    interface : integer or None
        The interface can be used to remove points from the path so
        that the path will be valid.
    points : integer, optional
        The number of points to add to the path.

    """
    xxx = [start[0], maxorder[0], end[0]]
    yyy = [start[1], maxorder[1], end[1]]
    warnings.simplefilter('ignore', getattr(np, 'RankWarning', Warning))
    par = np.polyfit(xxx, yyy, 2)
    xre = np.linspace(0., xxx[-1], points)
    yre = np.polyval(par, xre)
    # Delete some points from "yre" so that the path will be ok:
    if interface is not None:
        idx = [0]
        if yre[0] < interface:
            for i in np.where(yre > interface)[0]:
                idx.append(i)
        else:
            for i in np.where(yre < interface)[0]:
                idx.append(i)
        idx.append(-1)
        yre2 = [yre[i] for i in idx]
    else:
        yre2 = yre
    path = Path(rgen=MockRandomGenerator(seed=0))
    previous = None
    for order in yre2:
        if previous is None:
            vel = 0.0
        else:
            vel = order - previous
        ekin = 0.5 * vel**2
        phasepoint = make_internal_system(
            order=[order],
            pos=np.array([[order, ]]),
            vel=np.array([[vel, ]]),
            vpot=order**2,
            ekin=ekin
        )
        path.append(phasepoint)
        previous = order
    path.generated = ('fake', 0, 0, 0)
    path.maxlen = 10000

    return path


def create_ensembles_and_paths(task='retis', total_eclipse=20, number=0):
    """Return some test data we can use."""
    settings = {
        'simulation': {
            'task': task,
            'exe_path': HERE,
            'interfaces': [-1., 0., 1., 2., 10] if task == 'retis'
            else [-1., 0., 1., 2., 5., 10.],
            'zero_ensemble': True,
            'flux': True
        },
        'system': {
            'units': 'lj',
            'temperature': 0.1
        },
        'particles': {
            'type': 'internal',
            'dim': 1,
            'npart': 1,
            'position': {'input_file': os.path.join(HERE, 'config.xyz')}
        },
        'engine': {
            'class': 'Verlet',
            'timestep': 1
        },
        'orderparameter': {
            'class': 'position',
            'index': 1
        },
        'output': {'restart-file': 10000},
        'tis': {'detect': 0,
                'sigma_v': -1,
                'aimless': True,
                'zero_momentum': False,
                'rescale_energy': False,
                'allowmaxlength': False},
        'retis': {
            'nullmoves': True,
            'swapsimul': True,
            'swapfreq': 0.5
        },
    }

    with patch('sys.stdout', new=StringIO()):
        fill_up_tis_and_retis_settings(settings)
        check_ensemble(settings)
        ensembles = create_ensembles(settings)

    for i_ens, ens in enumerate(ensembles):
        system = System()
        particles = Particles(dim=1)
        particles.add_particle(np.zeros((1, 1)), np.zeros((1, 1)),
                               np.zeros((1, 1)))
        system.particles = particles
        system.forcefield = ForceField('empty force field',
                                       potential=[MockPotential()])

        ens['order_function'] = MockOrder()
        ens['engine'] = MockEngine(timestep=1.0, turn_around=total_eclipse)
        ens['system'] = system
        ens['rgen'] = MockRandomGenerator(seed=0)
        settings['ensemble'][i_ens]['tis']['maxlength'] = 1000

    if task == 'repptis':
        if number == 0:
            # We have 6 ensembles, where we create
            # the following acceptable paths:
            # [0^-]---RMR
            # [0^+-']-LML
            # [1^+-]--LMR
            # [2^+-]--LML
            # [3^+-]--RMR
            # [4^+-]--LML
            starts = [(0, -0.91), (0, -1.12), (0, -1.13), (0, -.123),
                      (0, 5.05), (0, 1.98)]
            ends = [(100, -0.95), (100, -1.05), (100, 1.2), (100, -0.123),
                    (100, 5.04), (100, 1.96)]
            maxs = [(50, -5.02), (50, -0.053), (50, 0.122), (50, 1.051),
                    (50, 1.973), (50, 5.48)]
        elif number == 1:
            # We have 6 ensembles, where we create
            # the following acceptable paths:
            # [0^-]---RMR
            # [0^+-']-RML
            # [1^+-]--LMR
            # [2^+-]--LML
            # [3^+-]--RMR
            # [4^+-]--LML
            starts = [(0, -0.91), (0, 0.031), (0, -1.13), (0, -.123),
                      (0, 5.05), (0, 1.98)]
            ends = [(100, -0.95), (100, -1.05), (100, 1.2), (100, -0.123),
                    (100, 5.04), (100, 1.96)]
            maxs = [(50, -5.02), (50, -0.53), (50, 0.122), (50, 1.051),
                    (50, 1.973), (50, 5.48)]
        else:
            raise NotImplementedError("Only number=0 or number=1 setups exist")
        for i, j, k, ens in zip(starts, ends, maxs, ensembles):
            path_ensemble = ens['path_ensemble']
            path = make_internal_path(i, j, k)
            path_ensemble.add_path_data(path, 'ACC')
        # plot the paths
        # import matplotlib.pyplot as plt
        # fig, ax = plt.subplots()
        # for ens in ensembles:
        #     path_ = ens['path_ensemble'].last_path
        #     phpath = [ph.order[0] for ph in path_.phasepoints]
        #     ax.plot(phpath, label=ens['path_ensemble'].ensemble_name)
        # # plot the interfaces
        # for i in settings['simulation']['interfaces']:
        #     ax.axhline(i, ls='--', color='k')
        # ax.legend()
        # fig.savefig('test_paths.png')
        return settings, ensembles

    # Make some paths for the ensembles.
    starts = [(0, -0.9), (0, -1.1), (0, -1.05), (0, -1.123), (0, -1.01)]
    ends = [(100, -0.95), (100, -1.2), (100, -1.3), (100, -1.01),
            (100, 10.01)]
    maxs = [(50, -5), (50, -0.2), (50, 0.5), (50, 2.5), (100, 10.01)]
    for i, j, k, ens in zip(starts, ends, maxs, ensembles):
        path_ensemble = ens['path_ensemble']
        path = make_internal_path(i, j, k, path_ensemble.interfaces[1])
        path_ensemble.add_path_data(path, 'ACC')

    return settings, ensembles


def prepare_test_simulation(settings=None):
    """Prepare a small system we can integrate."""
    if settings is None:
        settings = prepare_test_settings()

    box = create_box(periodic=[False])
    system = System(temperature=settings['system']['temperature'],
                    units=settings['system']['units'], box=box)
    system.particles = Particles(dim=system.get_dim())
    system.forcefield = create_force_field(settings)
    system.add_particle(np.array([-1.0]), mass=1, name='Ar', ptype=0)
    check_engine(settings)
    engine = create_engine(settings)
    settings['system']['obj'] = system
    settings['engine']['obj'] = engine

    fill_up_tis_and_retis_settings(settings)
    check_ensemble(settings)
    with patch('sys.stdout', new=StringIO()):
        simulation = create_tis_simulation(settings)
    # here we do a hack so that the simulation and Langevin integrator
    # both use the same random generator:
    simulation.ensembles[0]['rgen'] =\
        MockRandomGenerator(settings['tis']['seed'])
    simulation.ensembles[0]['system'] = system
    simulation.ensembles[0]['order_function'] = MockOrder()
    system.particles.vel = np.array([[0.78008019924163818]])
    return simulation, settings


def prepare_test_settings(simtype='tis'):
    """Prepare a good set of settings"""
    settings = {}
    # Basic settings for the simulation:
    settings['simulation'] = {'task': simtype,
                              'steps': 10,
                              'exe_path': '.',
                              'interfaces': [-0.9, -0.9, 1.0],
                              'zero_ensemble': False,
                              'flux': False}
    settings['system'] = {'units': 'lj',
                          'temperature': 0.07}
    # Basic settings for the Langevin integrator:
    settings['engine'] = {'class': 'Langevin',
                          'gamma': 0.3,
                          'high_friction': False,
                          'seed': 1,
                          'exe_path': '.',
                          'rgen': 'rgen',
                          'timestep': 0.002}
    # Fake initial position for particles:
    settings['particles'] = {'type': 'internal',
                             'position': {'generate': 'fcc',
                                          'repeat': [1, 2, 3],
                                          'density': 0.9
                                          }}
    # Potential parameters:
    # The potential is: `V_\text{pot} = a x^4 - b (x - c)^2`.
    settings['potential'] = [{'a': 1.0, 'b': 2.0, 'c': 0.0,
                              'class': 'DoubleWell'}]
    # Settings for the order parameter:
    settings['orderparameter'] = {'class': 'Position',
                                  'dim': 'x', 'index': 0, 'name': 'Position',
                                  'periodic': False}
    # TIS specific settings:
    settings['tis'] = {'freq': 0.5,
                       'maxlength': 20000,
                       'aimless': True,
                       'allowmaxlength': False,
                       'sigma_v': -1,
                       'seed': 1,
                       'rgen': 'mock',
                       'zero_momentum': False,
                       'rescale_energy': False}
    settings['initial-path'] = {'method': 'kick'}

    return settings
