# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Define some common methods for the tests."""
from contextlib import contextmanager
import logging
import math
import os
import numpy as np
from numpy.random import randint, random
from pyretis.core.random_gen import MockRandomGenerator
from pyretis.core.system import System
from pyretis.core.particles import Particles, ParticlesExt
from pyretis.core.path import Path


DATA_CROSS = [[(773, 5, '+')], [(1326, 9, '-')], [(1889, 2, '+')],
              [(2761, 3, '+')], [(2966, 9, '+')], [(3233, 2, '+')],
              [(3708, 7, '-')], [(4580, 4, '+')], [(5045, 8, '+')],
              [(5672, 5, '+')]]

DATA_CROSS2 = [(773, 6, 1), (1326, 10, -1), (1889, 3, 1), (2761, 4, 1),
               (2966, 10, 1), (3233, 3, 1), (3708, 8, -1), (4580, 5, 1),
               (5045, 9, 1), (5672, 6, 1)]

CORRECT_CROSS = [
    '       773    6   +',
    '      1326   10   -',
    '      1889    3   +',
    '      2761    4   +',
    '      2966   10   +',
    '      3233    3   +',
    '      3708    8   -',
    '      4580    5   +',
    '      5045    9   +',
    '      5672    6   +',
]

DATA_ENERGY = [
    (0.873, 0.477, 0.262, 0.615),
    (0.219, 0.740, 0.786, 0.147),
    (0.235, 0.764, 0.454, 0.173),
    (0.679, 0.098, 0.130, 0.110),
    (0.442, 0.374, 0.606, 0.954),
    (0.634, 0.774, 0.168, 0.314),
    (0.732, 0.250, 0.619, 0.691),
    (0.920, 0.607, 0.515, 0.487),
    (0.856, 0.105, 0.694, 0.233),
    (0.384, 0.432, 0.528, 0.155),
]

CORRECT_ENERGY = [
    '         0       0.873000       0.477000       0.262000       0.615000',
    '         1       0.219000       0.740000       0.786000       0.147000',
    '         2       0.235000       0.764000       0.454000       0.173000',
    '         3       0.679000       0.098000       0.130000       0.110000',
    '         4       0.442000       0.374000       0.606000       0.954000',
    '         5       0.634000       0.774000       0.168000       0.314000',
    '         6       0.732000       0.250000       0.619000       0.691000',
    '         7       0.920000       0.607000       0.515000       0.487000',
    '         8       0.856000       0.105000       0.694000       0.233000',
    '         9       0.384000       0.432000       0.528000       0.155000',
]

DATA_ORDER = [
    (0.90415742, 0.99887650, 0.15510997, 0.50120026),
    (0.90798451, 0.74665316, 0.67743021, 0.88284717),
    (0.72793238, 0.26006559, 0.57081755, 0.80117034),
    (0.37518802, 0.58331773, 0.80752662, 0.19093917),
    (0.93618894, 0.01188808, 0.07401551, 0.92200559),
    (0.35298599, 0.51682127, 0.00243432, 0.63432985),
    (0.01989170, 0.86621268, 0.68238818, 0.34800920),
    (0.09492324, 0.52609810, 0.88487148, 0.56365333),
    (0.22137013, 0.63323726, 0.44261574, 0.86234165),
    (0.52925374, 0.08251509, 0.59516637, 0.87529531),
]

CORRECT_ORDER = [
    '         0     0.904157     0.998876     0.155110     0.501200',
    '         1     0.907985     0.746653     0.677430     0.882847',
    '         2     0.727932     0.260066     0.570818     0.801170',
    '         3     0.375188     0.583318     0.807527     0.190939',
    '         4     0.936189     0.011888     0.074016     0.922006',
    '         5     0.352986     0.516821     0.002434     0.634330',
    '         6     0.019892     0.866213     0.682388     0.348009',
    '         7     0.094923     0.526098     0.884871     0.563653',
    '         8     0.221370     0.633237     0.442616     0.862342',
    '         9     0.529254     0.082515     0.595166     0.875295',
]

CORRECT_PATH_EXT = [
    '# Cycle: 0, status: ACC',
    '#     Step              Filename       index    vel',
    '         0           initial.g96           0      1',
    '         1             trajB.trr           5     -1',
    '         2             trajB.trr           4     -1',
    '         3             trajB.trr           3     -1',
    '         4             trajB.trr           2     -1',
    '         5             trajB.trr           1     -1',
    '         6             trajF.trr           0      1',
    '         7             trajF.trr           1      1',
    '         8             trajF.trr           2      1',
    '         9             trajF.trr           3      1',
    '        10             trajF.trr           4      1',
]

PATH_EXT_RAW = [
    ['0', 'initial.g96', '0', '1'],
    ['1', 'trajB.trr', '5', '-1'],
    ['2', 'trajB.trr', '4', '-1'],
    ['3', 'trajB.trr', '3', '-1'],
    ['4', 'trajB.trr', '2', '-1'],
    ['5', 'trajB.trr', '1', '-1'],
    ['6', 'trajF.trr', '0', '1'],
    ['7', 'trajF.trr', '1', '1'],
    ['8', 'trajF.trr', '2', '1'],
    ['9', 'trajF.trr', '3', '1'],
    ['10', 'trajF.trr', '4', '1'],
]


@contextmanager
def turn_on_logging():
    """Turn on logging so that tests can detect it."""
    logging.disable(logging.NOTSET)
    try:
        yield
    finally:
        logging.disable(logging.CRITICAL)


def compare_lists(list1, list2, cmps='equal'):
    """Just comare some items from a list."""
    if len(list1) != len(list2):
        raise AssertionError('Different length of lists')
    if cmps == 'float':
        for i, j in zip(list1, list2):
            if not math.isclose(i, j, abs_tol=1e-8):
                raise AssertionError(f'Different {i} != {j}')
    else:
        for i, j in zip(list1, list2):
            if not i == j:
                raise AssertionError(f'Different {i} != {j}')


def assert_equal_path_dict(dict1, dict2):
    """Compare two path dictionaries."""
    comp = {
        'cycle': 'int',
        'generated': 'tuple',
        'interface': 'tuple',
        'length': 'int',
        'ordermax': 'tuple-float',
        'ordermin': 'tuple-float',
        'status': 'string',
        'weight': 'float'
    }
    for key in dict1:
        if key not in dict2:
            raise AssertionError('Different keys in dicts!')
    for key in dict2:
        if key not in dict1:
            raise AssertionError('Different keys in dics!')
    for key, val in comp.items():
        if key not in dict1 or key not in dict2:
            raise AssertionError(f'Missing key: {key}')
        if val in ('int', 'string', 'float'):
            if not dict1[key] == dict2[key]:
                raise AssertionError(f'Different {key}')
        elif val in ('tuple', 'tuple-float'):
            if val != 'tuple-float':
                compare_lists(dict1[key], dict2[key])
            else:
                compare_lists(dict1[key], dict2[key], cmps='float')


def _add_path_data(path_data, phasepoint):
    """Add path data for storage."""
    if not path_data:
        step = 0
    else:
        step = path_data[-1][0] + 1
    new_data = [step, phasepoint.particles.get_pos()[0],
                phasepoint.particles.get_pos()[1],
                phasepoint.particles.get_vel()]
    path_data.append(new_data)


def set_up_system(order, pos, vel, vpot=None, ekin=None, internal=True):
    """Create a system for testing."""
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


def create_external_path(random_length=False):
    """Create an external path for testing."""
    path_data = []
    path = Path(None)
    phasepoint = set_up_system([random(), random()],
                               ('initial.g96', None),
                               False, internal=False)
    path.append(phasepoint)
    _add_path_data(path_data, phasepoint)
    if random_length:
        length = randint(10, 100)
    else:
        length = 5
    for i in range(length, 0, -1):
        phasepoint = set_up_system([random(), random()],
                                   ('trajB.trr', i), True,
                                   vpot=random(), ekin=random(),
                                   internal=False)
        path.append(phasepoint)
        _add_path_data(path_data, phasepoint)
    for i in range(0, length):
        phasepoint = set_up_system([random(), random()],
                                   ('trajF.trr', i), False,
                                   vpot=random(), ekin=random(),
                                   internal=False)
        path.append(phasepoint)
        _add_path_data(path_data, phasepoint)
    return path, path_data


def create_test_paths(npath=5):
    """Create some paths we can use.

    Parameters
    ----------
    npath : integer
        The number of paths to create.

    """
    paths = []
    correct_data = []
    for i in range(npath):
        length = randint(10, 100)
        path = Path(rgen=MockRandomGenerator(seed=0))
        data_copy = {'order': [], 'pos': [], 'vel': [],
                     'ekin': [], 'vpot': []}
        for j in range(length):
            phasepoint = set_up_system([(i + 1) * j, j, j*10],
                                       np.ones((9, 3)) * i * j,
                                       np.ones((9, 3)) * i * j * -1,
                                       vpot=(i + j),
                                       ekin=(i + j + 1))
            path.append(phasepoint)
            data_copy['order'].append(phasepoint.order)
            for key in ('pos', 'vel', 'ekin', 'vpot'):
                data_copy[key].append(getattr(phasepoint.particles, key))
        paths.append(path)
        correct_data.append(data_copy)
    return paths, correct_data


def remove_file(filename):
    """Remove the given file."""
    try:
        os.remove(filename)
    except OSError:
        pass
