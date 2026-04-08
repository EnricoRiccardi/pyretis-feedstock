# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the classes and methods from pyretis.core.box"""
import pytest
import numpy as np
from pyretis.core.box import (
    create_box,
    RectangularBox,
    TriclinicBox,
    box_matrix_to_list,
    box_vector_angles,
    angles_from_box_matrix,
    box_from_restart
)


def test_create_empty_box():
    """Test the creation of boxes with no arguments."""
    box = create_box()
    assert isinstance(box, RectangularBox)
    assert box.low == [-float('inf')]
    assert box.high == [float('inf')]
    assert box.length == [float('inf')]
    assert box.periodic == [False]


def test_create_missing_periodic():
    """Test default behaviour of creation without periodic arguments."""
    box = create_box(cell=[10, 10, 10], periodic=None)
    assert isinstance(box, RectangularBox)
    assert box.periodic == [True, True, True]
    assert np.allclose(box.low, [0., 0., 0.])
    assert np.allclose(box.high, [10., 10., 10.])
    assert np.allclose(box.length, [10., 10., 10.])

    box = create_box(cell=[10, 10, 10], periodic=[False])
    assert box.periodic == [False, True, True]


def test_create_box_fail():
    """Test some corner cases when we fail to create a box."""
    # Test that we fail for too many periodic values given:
    with pytest.raises(ValueError):
        create_box(cell=[10, 10, 10], periodic=[False]*4)
    # Test that we fail when low == high:
    with pytest.raises(ValueError):
        create_box(cell=[10, 10], low=[0, 10], high=[10, 10])
    # Test that we fail for inconsistent lengths:
    with pytest.raises(ValueError):
        create_box(cell=[10, 10], low=[0, 0], high=[10, 20])


def test_box_size():
    """Test that giving a size works as expected."""
    test_in = (
        {'length': [10, 10]},
        {'low': [0., -10.], 'high': [10., 20.]},
        {'low': [1, 1, 1], 'length': [10, 11, 12]},
        {'high': [1, 1, 1], 'length': [10, 11, 12]},
        {'high': [10, 9, 8]},
        {'low': [-10, 10]},
    )
    correct = (
        {'low': [0., 0.], 'high': [10., 10.], 'length': [10., 10.]},
        {'low': [0., -10.], 'high': [10., 20.], 'length': [10., 30.]},
        {'low': [1, 1, 1], 'high': [11., 12., 13.],
         'length': [10., 11., 12.]},
        {'low': [-9, -10, -11], 'high': [1., 1., 1.],
         'length': [10., 11., 12.]},
        {'low': [0., 0., 0.], 'high': [10, 9, 8], 'length': [10, 9, 8]},
        {'low': [-10, 10], 'high': [float('inf'), float('inf')],
         'length': [float('inf'), float('inf')]},
    )
    for case, corr in zip(test_in, correct):
        box = create_box(low=case.get('low', None),
                         high=case.get('high', None),
                         cell=case.get('length', None),
                         periodic=case.get('periodic', None))
        assert isinstance(box, RectangularBox)
        assert np.allclose(box.length, corr['length'])
        assert np.allclose(box.low, corr['low'])
        assert np.allclose(box.high, corr['high'])


def test_volume_calculate():
    """Test calculation of volume."""
    box = create_box(cell=[10])
    assert isinstance(box, RectangularBox)
    assert box.calculate_volume() == pytest.approx(10.)
    box2 = create_box(cell=[10, 10])
    assert isinstance(box2, RectangularBox)
    assert box2.calculate_volume() == pytest.approx(100.)
    box3 = create_box(cell=[10, 10, 10])
    assert isinstance(box3, RectangularBox)
    assert box3.calculate_volume() == pytest.approx(1000.)
    box4 = create_box(low=[-10, 0, 11], high=[10, 5, 19])
    assert isinstance(box4, RectangularBox)
    assert box4.calculate_volume() == pytest.approx(20.*5.*8.)


def test_faulty_input():
    """Test that the initialisation fails as we expect."""
    with pytest.raises(ValueError):
        create_box(cell=[10, -10, 10])
    with pytest.raises(TypeError):
        create_box(cell=10)
    with pytest.raises(ValueError):
        create_box(cell=[10, (-10, 10), (0, -15)])
    with pytest.raises(ValueError):
        create_box(low=[0, 0], high=[10, 10], cell=[11, 11])
    with pytest.raises(ValueError):
        create_box(cell=['crash', 15])
    with pytest.raises(ValueError):
        create_box(low=[10, 10], high=[10, 10])


def test_update_box():
    """Test update of box size."""
    box = create_box(cell=[10, 10, 10])
    new_length = [10, 11, 12]
    box.update_size(new_length)
    for i, j in zip(box.length, new_length):
        assert i == pytest.approx(j)
    new_length2 = [13, 12, 11, 10]
    box.update_size(new_length2)  # This should NOT update.
    for i, j in zip(box.length, new_length):
        assert i == pytest.approx(j)
    new_length3 = [3, 3]
    box.update_size(new_length3)  # This should NOT update.
    for i, j in zip(box.length, new_length):
        assert i == pytest.approx(j)
    new_length4 = None
    box.update_size(new_length4)
    for i, j in zip(box.length, new_length):
        assert i == pytest.approx(j)


def test_bounds():
    """Test the bounds method."""
    box = create_box(cell=[10, 11, 12])
    correct = [[0., 10.], [0., 11.], [0., 12.]]
    bounds = box.bounds()
    for bound, corr in zip(bounds, correct):
        for i, j in zip(bound, corr):
            assert i == pytest.approx(j)


def test_print_length():
    """Test that we print out cell parameters correctly."""
    lengths = [
        [1],
        [1, 2],
        [1, 2, 3],
        [1, 2, 3, 4, 5, 6],
        [1, 2, 3, 4, 5, 6, 7, 8, 9],
    ]
    fmt = '{:6.2f}'
    for i in lengths:
        box = create_box(cell=i)
        correct = ' '.join([fmt.format(j) for j in i])
        assert correct == box.print_length(fmt=fmt)


def test_restart():
    """Test that we create and read restart info for a box."""
    box_settings = [
        {'low': [10, -1, 101], 'high': [12, 5, 102],
         'periodic': [True, True, False], 'cell': [2., 6., 1.]},
        {'low': [9, 8, 7], 'periodic': [False, True, False],
         'cell': [1, 2, 3, 4, 5, 6]},
    ]
    keys = ('cell', 'periodic', 'low', 'high')
    for setting in box_settings:
        box = create_box(**setting)
        restart = box.restart_info()
        assert all(k in restart for k in keys)
        for key, val in restart.items():
            if key not in setting:
                continue
            if key == 'periodic':
                assert val == setting[key]
            else:
                assert np.allclose(val, setting[key])
        info = {'box': restart}
        box2 = box_from_restart(info)
        # todo The line under or above should be removed once the restart
        # strategy is unified.
        box2.load_restart_info(restart)
        assert box == box2
        assert not box != box2


def test_pbc_coordinate_dim():
    """Test pbc for a specific dimension."""
    length = [10, 11, 12]
    box = create_box(cell=length, periodic=[False, True, True])
    pos = [11, 10, 14]
    pbc_pos = box.pbc_coordinate_dim(pos[2], 2)
    assert pbc_pos == pytest.approx(2.0)
    pbc_pos = box.pbc_coordinate_dim(pos[0], 0)
    assert pbc_pos == pytest.approx(pos[0])


def test_pbc_wrap():
    """Test pbc wrap for coordinates."""
    length = [10, 11, 12]
    box = create_box(cell=length, periodic=[False, True, True])
    pos = np.array([[11, 10, 14], ])
    correct = np.array([[11, 10, 2], ])
    pbc_pos = box.pbc_wrap(pos)
    assert np.allclose(correct, pbc_pos)
    assert not np.allclose(correct, pos)


def test_pbc_dist_matrix():
    """Test pbc wrap for a distance vector."""
    box = create_box(cell=[10, 10, 10], periodic=[False, True, True])
    dist = np.array([[8., 7., 9.], ])
    pbc_dist = box.pbc_dist_matrix(dist)
    correct = np.array([[8., -3., -1.], ])
    assert np.allclose(correct, pbc_dist)
    assert np.allclose(correct, dist)


def test_triclinic_create_box():
    """Test creation of TriclinicBox."""
    box1 = create_box(
        cell=[17.5092040633036, 7.58170825120892, 6.95903807579504,
              4.37730063346742, 0.0, 0.0]
    )
    assert isinstance(box1, TriclinicBox)
    box2 = create_box(cell=[10, 10, 10, 0.0, 0.0, 0.0])
    assert isinstance(box2, TriclinicBox)
    box3 = create_box(cell=[1, 2, 3, 4, 5, 6, 7, 8, 9])
    assert isinstance(box3, TriclinicBox)


def test_triclinic_volume_calculate():
    """Test calculation of volume."""
    box1 = create_box(cell=[10, 1, 1, 0, 0, 0])
    assert box1.calculate_volume() == pytest.approx(10.)
    box2 = create_box(cell=[10, 11, 12, 0, 0, 0])
    # wait box 1 was used?
    assert box1.calculate_volume() == pytest.approx(10.)
    # Original test was:
    # box2 = create_box(cell=[10, 11, 12, 0, 0, 0])
    # assert box2.calculate_volume() == 10.*11.*12.
    assert box2.calculate_volume() == pytest.approx(10.*11.*12.)
    length = [17.5092040633036, 7.58170825120892, 6.95903807579504,
              4.37730063346742, 0.0, 0.0]
    box3 = create_box(cell=length)
    assert box3.calculate_volume() == pytest.approx(923.810056228)
    length = [17.5092040633036, 7.58170825120892, 6.95903807579504,
              4.37730063346742, 0.0, 0.0, 0.0, 0.0, 0.0]
    box4 = create_box(cell=length)
    assert box4.calculate_volume() == pytest.approx(923.810056228)


def test_triclinic_update_box():
    """Test update for triclinic box."""
    box = create_box(cell=[1, 2, 3, 0, 0, 0])
    new_size = [1, 2, 3, 4, 5, 6]
    box.update_size(new_size)
    correct_size = np.array([[1., 4., 5.], [0., 2., 6.], [0., 0, 3.]])
    assert np.allclose(box.box_matrix, correct_size)
    new_size = [-1, -1, -1]
    box.update_size(new_size)  # This should NOT update size.
    assert np.allclose(box.box_matrix, correct_size)
    new_size = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    box.update_size(new_size)
    correct_size = np.array([[1., 4., 5.], [6., 2., 7.], [8., 9., 3.]])
    assert np.allclose(box.box_matrix, correct_size)


def test_box_matrix_to_list():
    """Test that we return the box matrix as a list as expected."""
    new_size = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    box = create_box(cell=new_size)
    out = box_matrix_to_list(box.box_matrix)
    for i, j in zip(new_size, out):
        assert i == pytest.approx(j)
    assert box_matrix_to_list(None) is None


def test_box_vector_angles():
    """Test conversion from a,b,c, alpha, beta, gamma."""
    test_data = [
        {'length': [3., 3., 3.],
         'alpha': 100., 'beta': 80., 'gamma': 75.,
         'box': np.array([[3.0, 0.77646, 0.52094],
                          [0.0, 2.89778, -0.67891],
                          [0.0, 0.0, 2.87536]])},
        {'length': [3., 5.1, 1.9],
         'alpha': 100., 'beta': 85., 'gamma': 66.,
         'box': np.array([[3.0, 2.07436, 0.16559],
                          [0.0, 4.65908, -0.43488],
                          [0.0, 0.0, 1.84213]])},
    ]
    for i in test_data:
        box_matrix = box_vector_angles(i['length'], i['alpha'],
                                       i['beta'], i['gamma'])
        assert np.allclose(box_matrix, i['box'], atol=1e-4)


def test_angles_from_box_matrix():
    """Test that we can get box angles and lengths from a box matrix."""
    test_data = [
        {'length': [3., 3., 3.],
         'alpha': 100., 'beta': 80., 'gamma': 75.,
         'box': np.array([[3.0, 0.77646, 0.52094],
                          [0.0, 2.89778, -0.67891],
                          [0.0, 0.0, 2.87536]])},
        {'length': [3., 5.1, 1.9],
         'alpha': 100., 'beta': 85., 'gamma': 66.,
         'box': np.array([[3.0, 2.07436, 0.16559],
                          [0.0, 4.65908, -0.43488],
                          [0.0, 0.0, 1.84213]])},
    ]
    for i in test_data:
        length, alpha, beta, gamma = angles_from_box_matrix(i['box'])
        assert np.allclose(length, i['length'], atol=1e-5)
        assert alpha == pytest.approx(i['alpha'], abs=1e-3)
        assert beta == pytest.approx(i['beta'], abs=1e-3)
        assert gamma == pytest.approx(i['gamma'], abs=1e-3)


def test_box_equality():
    """Test that we can compare boxes."""
    box1 = create_box(cell=[1, 2, 3, 4, 5, 6, 7, 8, 9])
    box2 = create_box(cell=[1, 2, 3, 4, 5, 6, 7, 8, 9])
    assert box1 == box2
    # Test failure for different periodic settings:
    box1.periodic = [False, True, True]
    assert box1 != box2
    box1.periodic = [False]
    assert box1 != box2
    box1.periodic = [True, True, True]
    assert box1 == box2
    box1.length = None
    assert box1 != box2
    box2.length = None
    assert box1 == box2
    if hasattr(box2, 'length'):
        del box2.length
    assert box1 != box2
    box2.length2 = 100
    assert box1 != box2
    box3 = create_box(cell=[1, 2])
    assert box1 != box3
    box4 = create_box(cell=[1, 2])
    assert box3 == box4
    box3.dim = 100
    assert box3 != box4


def test_box_copy():
    """Test that we can copy a box."""
    box1 = create_box(cell=[1, 2, 3, 4, 5, 6, 7, 8, 9])
    box2 = box1.copy()
    assert box1 is not box2
    assert box1 == box2
    box1 = create_box(cell=[1, 2, 3])
    box2 = box1.copy()
    assert box1 is not box2
    assert box1 == box2
