# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test methods from pyretis.core.common"""
import os
import sys
import inspect
import logging
import pytest
import numpy as np
from io import StringIO
from pyretis.core.common import (
    big_fat_comparer,
    crossing_finder,
    generic_factory,
    import_from,
    inspect_function,
    numpy_allclose,
    segments_counter,
    select_and_trim_a_segment,
    trim_path_between_interfaces,
    wirefence_weight_and_pick,
    _pick_out_arg_kwargs,
)
from pyretis.core.random_gen import RandomGenerator
from unittest.mock import patch
from .test_path import PATHTEST0, PATHTEST1, PATHTEST2

LOCAL_DIR = os.path.abspath(os.path.dirname(__file__))


# Define some functions used for testing:
# pylint: disable=unused-argument
def function1(arg1, arg2, arg3, arg4):
    """To test positional arguments."""
# pylint: enable=unused-argument


# pylint: disable=unused-argument
def function2(arg1, arg2, arg3, arg4=10):
    """To test positional and keyword arguments."""
# pylint: enable=unused-argument


# pylint: disable=unused-argument
def function3(arg1, arg2, arg3, arg4=100, arg5=10):
    """To test positional and keyword arguments."""
# pylint: enable=unused-argument


# pylint: disable=unused-argument
def function4(*args, **kwargs):
    """To test positional and keyword arguments."""
# pylint: enable=unused-argument


# pylint: disable=unused-argument
def function5(arg1, arg2, *args, arg3=100, arg4=100):
    """To test positional and keyword arguments."""
# pylint: enable=unused-argument


# pylint: disable=unused-argument
def function6(arg1, arg2, arg3=100, *, arg4=10):
    """To test positional and keyword arguments."""
# pylint: enable=unused-argument


# pylint: disable=unused-argument, keyword-arg-before-vararg
def function7(arg1, arg2, arg3=100, *args, arg4, arg5=10):
    """To test positional and keyword arguments."""
# pylint: enable=unused-argument, keyword-arg-before-vararg


# pylint: disable=unused-argument
def function8(arg1, arg2=100, self='something'):
    """To test redefinition of __init__."""
# pylint: enable=unused-argument


CORRECT = [
    {'args': ['arg1', 'arg2', 'arg3', 'arg4'],
     'varargs': [], 'kwargs': [], 'keywords': []},
    {'args': ['arg1', 'arg2', 'arg3'],
     'varargs': [], 'kwargs': ['arg4'], 'keywords': []},
    {'args': ['arg1', 'arg2', 'arg3'], 'varargs': [],
     'kwargs': ['arg4', 'arg5'], 'keywords': []},
    {'args': [], 'varargs': ['args'], 'kwargs': [],
     'keywords': ['kwargs']},
    {'args': ['arg1', 'arg2'], 'kwargs': ['arg3', 'arg4'],
     'varargs': ['args'], 'keywords': []},
    {'args': ['arg1', 'arg2'], 'kwargs': ['arg3', 'arg4'],
     'varargs': [], 'keywords': []},
    {'args': ['arg1', 'arg2'], 'kwargs': ['arg3', 'arg4', 'arg5'],
     'varargs': ['args'], 'keywords': []},
]


def test_import():
    """Test that we can import."""
    module = os.path.join(LOCAL_DIR, 'fooengine.py')
    klass = 'FooEngine'
    imp = import_from(module, klass)
    sys.path.insert(0, LOCAL_DIR)
    from fooengine import FooEngine
    del sys.path[0]
    assert inspect.getsource(imp) == inspect.getsource(FooEngine)


def test_import_importerror():
    """Test that we can import."""
    module = os.path.join(LOCAL_DIR, 'fooengine.py')
    klass = 'DoesNotExist'
    with pytest.raises(ValueError):
        import_from(module, klass)
    module = 'path-that-should-not-be'
    with pytest.raises(ValueError):
        import_from(module, klass)
    module = os.path.join(LOCAL_DIR, 'isnothere.py')
    with pytest.raises(ValueError):
        import_from(module, klass)


def test_inspect():
    """Test the inspect function method."""
    # Define some test functions:
    functions = [function1, function2, function3, function4,
                 function5, function6, function7]

    for i, func in enumerate(functions):
        args = inspect_function(func)
        assert args == CORRECT[i]


def test_arg_kind():
    """Test a function with only positional arguments."""
    args = inspect_function(range.__eq__)
    assert not args['keywords']
    assert not args['varargs']
    assert not args['kwargs']
    for i in ('self', 'value'):
        assert i in args['args']


def test_pick_out_kwargs():
    """Test pick out of "self" for kwargs."""
    settings = {'arg1': 10, 'arg2': 100, 'self': 'text'}

    # pylint: disable=too-few-public-methods
    class Abomination:
        """Just to allow redefination of __init__."""
    # pylint: enable=too-few-public-methods

    abo = Abomination()
    abo.__init__ = function8

    args, kwargs = _pick_out_arg_kwargs(abo, settings)
    assert 'self' not in args
    assert 'self' not in kwargs

    with pytest.raises(ValueError):
        _pick_out_arg_kwargs(abo, {})


class Klass1:
    """A class used for factory testing."""

    def __init__(self):
        self.stuff = 10

    def method1(self):
        """Return stuff."""
        return self.stuff

    def method2(self, add):
        """Add to stuff."""
        self.stuff += add


def test_factory():
    """Test that the factory works as intended."""
    factory_map = {'klass1': {'cls': Klass1}}

    settings = {'class': 'Klass1'}
    cls = generic_factory(settings, factory_map, name='Testing')
    assert isinstance(cls, Klass1)

    settings = {'clAsS': 'Klass1'}
    with patch('pyretis.core.common.logger.critical') as mock_critical:
        cls = generic_factory(settings, factory_map, name='Testing')
        assert cls is None
        assert mock_critical.called

    settings = {'Klass': 'Klass1'}
    with patch('pyretis.core.common.logger.critical') as mock_critical:
        cls = generic_factory(settings, factory_map, name='Testing')
        assert cls is None
        assert mock_critical.called


def test_numpy_allclose():
    """Test that the numpy allclose comparison works as intended."""
    val1 = np.array([1.2345, 2.34567, 9.87654321])
    val2 = np.array([1.2345, 2.34567, 9.87654321])
    val3 = np.array([1.2345, 2.34567, 10.])
    # Both are valid:
    assert numpy_allclose(val1, val2)
    assert not numpy_allclose(val1, val3)
    # Both are None:
    assert numpy_allclose(None, None)
    # One is None:
    assert not numpy_allclose(val1, None)
    assert not numpy_allclose(None, val1)
    # One is something else:
    assert not numpy_allclose(val1, 'the roof is on fire')
    # Both is something else:
    assert not numpy_allclose('the roof is on fire',
                              'the roof is on fire')


def test_fat_part_1():
    """Test the big_fat_comparer."""
    assert not big_fat_comparer(np.array([1]), np.array([1, 2, 3]))
    with pytest.raises(ValueError):
        big_fat_comparer(np.array([1]), np.array([1, 2, 3]), hard=True)

    assert not big_fat_comparer(np.array([1, 2, 4]),
                                np.array([1, 2, 3]))
    with pytest.raises(ValueError):
        big_fat_comparer(np.array([1, 2, 4]),
                         np.array([1, 2, 3]), hard=True)

    assert not big_fat_comparer(1, 2)
    with pytest.raises(ValueError):
        big_fat_comparer(1, 2, hard=True)

    assert not big_fat_comparer(1, [1])
    with pytest.raises(ValueError):
        big_fat_comparer(1, [1], hard=True)

    assert not big_fat_comparer((1, 3, 4), (1, 3))
    with pytest.raises(ValueError):
        big_fat_comparer((1, 3, 4), (1, 3), hard=True)
    assert not big_fat_comparer((1, 3), (1, 3, 4))
    with pytest.raises(ValueError):
        big_fat_comparer((1, 3), (1, 3, 4), hard=True)
    assert not big_fat_comparer((1, 2, 3), (2, 3, 4))
    with pytest.raises(ValueError):
        big_fat_comparer((1, 2, 3), (2, 3, 4), hard=True)

    assert not big_fat_comparer('Ciao', 'Mamma')
    with pytest.raises(ValueError):
        big_fat_comparer('Ciao', 'Mamma', hard=True)

    assert not big_fat_comparer([1, 2], [0])
    with pytest.raises(ValueError):
        big_fat_comparer([1, 2], [0], hard=True)

    a, b = np.array([0, 2]), np.array([0, 1])
    with pytest.raises(ValueError):
        big_fat_comparer(a, b, hard=True)

    a, b = [0, 2], [0, {'b': 'c'}]
    with pytest.raises(ValueError):
        big_fat_comparer(a, b, hard=True)


def test_fat_part_2():
    """Test the big_fat_comparer."""
    try1 = {'correct': CORRECT,
            'kkk': [1, 2, 3],
            'giorgio': {'santo': 'subito',
                        'yes': [12, 'b']}}

    try2 = {'correct': CORRECT,
            'kkk': [1, 2, 3],
            'giorgio': {'santo': 'subito',
                        'yes': [12, 'b']}}

    try3 = {'correct': CORRECT,
            'kkk': [1, 2, 3],
            'giorgio': {'santo': 'subit',
                        'yes': [12, 'b']}}
    try4 = {'correct': CORRECT,
            'kkk': [1, 2, 3],
            'giorgio': {'santo': 'subit',
                        'no': 'interest',
                        'yes': {'b': 12}}}

    assert not big_fat_comparer(try1, try3)
    with pytest.raises(ValueError):
        big_fat_comparer(try1, try3, hard=True)

    assert not big_fat_comparer(try1, try4)
    with pytest.raises(ValueError):
        big_fat_comparer(try1, try4, hard=True)

    assert not big_fat_comparer(try3, try4)
    with pytest.raises(ValueError):
        big_fat_comparer(try3, try4, hard=True)

    assert not big_fat_comparer(try2, try4)
    with pytest.raises(ValueError):
        big_fat_comparer(try2, try4, hard=True)

    assert not big_fat_comparer(
        try1, {**{'Ihate': 'Cinnamon'}, **try2}
    )
    with pytest.raises(ValueError):
        big_fat_comparer(try1, {**{'Ihate': 'Cinnamon'},
                                **try2}, hard=True)

    assert not big_fat_comparer({**{'Iate': 'Cinnamon'}, **try2},
                                try1)
    with pytest.raises(ValueError):
        big_fat_comparer({**{'Iate': 'Cinnamon'}, **try2},
                         try1, hard=True)

    assert not big_fat_comparer(try1, 'TestMcTest')
    with pytest.raises(ValueError):
        big_fat_comparer(try1, 'TestMcTest', hard=True)

    assert not big_fat_comparer('TestMcTest', try1)
    with pytest.raises(ValueError):
        big_fat_comparer('TestMcTest', try1, hard=True)

    assert big_fat_comparer(try1, try2)
    assert big_fat_comparer(try1, try2, hard=True)

    try2['giorgio']['yes'][0] = 1
    assert not big_fat_comparer(try1, try2)
    with pytest.raises(ValueError):
        big_fat_comparer(try1, try2, hard=True)

    try2['giorgio']['yes'][0] = 'c'
    assert not big_fat_comparer(try1, try2)
    with pytest.raises(ValueError):
        big_fat_comparer(try1, try2, hard=True)

    try2['giorgio']['yes'][0] = 12
    try1['giorgio']['yes'].append(12)
    assert not big_fat_comparer(try1, try2)
    with pytest.raises(ValueError):
        big_fat_comparer(try1, try2, hard=True)

    del try1['giorgio']['yes'][2]
    assert big_fat_comparer(try1, try2)
    assert big_fat_comparer(try1, try2, hard=True)


def test_x1_crossing_finder():
    path = PATHTEST0.copy()
    # Test downward crossing (at seed=1, index 1/3/5 is picked)
    path.rgen = RandomGenerator(seed=1)
    ph1, ph2 = crossing_finder(path=path, interface=2.9)
    assert ph1.order[0] == 15
    assert ph2.order[0] == 0

    # Test upward crossing (at seed=0, index 4 is picked)
    path.rgen = RandomGenerator(seed=0)
    ph1, ph2 = crossing_finder(path=path, interface=2.9)
    assert ph1.order[0] == 2
    assert ph2.order[0] == 3


def test_segments_counter():
    path = PATHTEST0.copy()
    assert 4 == segments_counter(path=path,
                                 interface_l=2.9, interface_r=4.1)
    assert 4 == segments_counter(path=path,
                                 interface_l=2, interface_r=9)


def test_trim_path_between_interfaces():
    """Test for trimming trajectories withing a range in path space."""
    path0 = trim_path_between_interfaces(
        path=PATHTEST0.copy(), interface_l=3.0,
        interface_r=7.0
    )
    path2 = PATHTEST2.copy()
    path2.generated = 'ct'

    with patch('sys.stdout', new=StringIO()):
        assert path0.length == 12
        assert path2 == path0
        assert PATHTEST0 != path0


def test_wirefence_weight_and_pick():
    """Test for wirefernce function in a pathological case."""
    path0 = trim_path_between_interfaces(
        path=PATHTEST0.copy(), interface_l=3.0,
        interface_r=7.0
    )
    # Create a path with L-L, L-R and R-R sub-paths between 1.5 and 2.5:
    for i, j in enumerate([1]+[2]*3+[1]+[2]*1+[3]+[2]*4+[3]):
        path0.phasepoints[i].order[0] = j

    q_tot, segment = wirefence_weight_and_pick(path=path0,
                                               intf_l=1.5, intf_r=2.5)
    # Only the L-L and L-R phase points to be counted
    assert q_tot == 4

    # Create a path with phase points i < intf_l and i+1 > intf_r,
    # and also a R-L sub-path:
    for i, j in enumerate([1]+[3]*6+[2]*4+[1]):
        path0.phasepoints[i].order[0] = j
    q_tot, segment = wirefence_weight_and_pick(path=path0,
                                               intf_l=1.5, intf_r=2.5,)
    # Only the R-L phase points to be counted
    assert q_tot == 4


def test_select_and_trim_a_segment():
    """Test to trim a selected segment keeping the crossing."""
    path = PATHTEST0.copy()
    path.rgen = RandomGenerator(seed=0)
    path_rnd = select_and_trim_a_segment(path=path, interface_l=4.1,
                                         interface_r=5.9)
    path3 = select_and_trim_a_segment(path=path.copy(), interface_l=4.1,
                                      interface_r=5.9, segment_to_pick=0)
    path1 = PATHTEST1.copy()
    path1.generated = 'sg'

    with patch('sys.stdout', new=StringIO()):
        assert path3.length == 3
        assert path1 == path3
        assert PATHTEST1 != path3

    path2 = path1.empty_path()
    path2.append(path1.phasepoints[0])
    path2.append(path1.phasepoints[1])

    path3 = select_and_trim_a_segment(path=path2, interface_l=4.1,
                                      interface_r=5.9, segment_to_pick=0)
    with patch('sys.stdout', new=StringIO()):
        assert path3 != path_rnd
        assert path3.length == 2
        assert PATHTEST1 != path3
