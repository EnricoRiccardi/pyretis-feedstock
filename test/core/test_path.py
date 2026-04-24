# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test functionality for the path classes from pyretis.core.path."""
import logging
import pytest
import numpy as np
from pyretis.core.path import (Path, paste_paths, check_crossing)
from pyretis.core.random_gen import RandomGenerator
from .help import (NegativeOrder, SameOrder, set_up_system,
                   PATHTEST0, PATHTEST1, PATHTEST2)


def fill_forward_backward(pathf, pathb, npoints=20):
    """Fill in data for forward and backward paths."""
    for i in range(npoints):
        pathf.append(
            set_up_system([1.0 * i], np.zeros(3) * i, np.zeros(3) * i, i, i)
        )
        pathb.append(
            set_up_system([-1.0 * i], np.zeros(3) * i * -1.0,
                          np.zeros(3) * i * -1.0, -1.0 * i, -1.0 * i)
        )  # Pretty exotic kinetic energy.


def test_paste_paths1():
    """Test that we can paste paths together."""
    rgen = RandomGenerator(seed=0)
    pathf = Path(rgen, maxlen=1000)
    pathb = Path(rgen, maxlen=1000)
    fill_forward_backward(pathf, pathb, npoints=20)
    path = paste_paths(pathb, pathf, overlap=False, maxlen=None)
    assert path.length == 40
    for i, phasepoint in enumerate(path.phasepoints):
        if i <= 19:
            assert phasepoint.order[0] == pytest.approx(i - 19.)
        else:
            assert phasepoint.order[0] == pytest.approx(i - 20.)
    path = paste_paths(pathb, pathf, overlap=True, maxlen=None)
    assert path.length == 39
    for i, phasepoint in enumerate(path.phasepoints):
        assert phasepoint.order[0] == pytest.approx(i - 19.)


def test_paste_paths2():
    """Test that we can paste paths together when we truncate."""
    rgen = RandomGenerator(seed=0)
    pathf = Path(rgen, maxlen=30)
    pathb = Path(rgen, maxlen=30)
    fill_forward_backward(pathf, pathb, npoints=20)
    path = paste_paths(pathb, pathf, overlap=True, maxlen=None)
    assert path.length == 30
    for i, phasepoint in enumerate(path.phasepoints):
        assert phasepoint.order[0] == pytest.approx(i - 19.)
    pathf = Path(rgen, maxlen=32)
    pathb = Path(rgen, maxlen=31)
    fill_forward_backward(pathf, pathb, npoints=20)
    path = paste_paths(pathb, pathf, overlap=True, maxlen=None)
    assert path.length == 32
    for i, phasepoint in enumerate(path.phasepoints):
        assert phasepoint.order[0] == pytest.approx(i - 19.)
    path = paste_paths(pathb, pathf, overlap=True, maxlen=10)
    assert path.length == 10
    for i, phasepoint in enumerate(path.phasepoints):
        assert phasepoint.order[0] == pytest.approx(i - 19.)


def test_paste_paths3():
    """Test that the backwards phasepoint survive paste-pathing."""
    rgen = RandomGenerator(seed=0)
    pathf = Path(rgen, maxlen=1000)
    pathb = Path(rgen, maxlen=1000)
    fill_forward_backward(pathf, pathb, npoints=20)
    pathb.phasepoints[0].order[0] = -11.11
    pathb.phasepoints[-1].order[0] = -22.22
    pathf.phasepoints[0].order[0] = 11.11
    pathf.phasepoints[-1].order[0] = 22.22
    path = paste_paths(pathb, pathf, overlap=True, maxlen=None)
    # The bacwkards path is the one that survives, while the
    # forwards path loses a phasepoint. Note that the backwards
    # path is being flipped in time.
    assert path.phasepoints[0].order[0] == pathb.phasepoints[-1].order[0]
    assert path.phasepoints[19].order[0] == pathb.phasepoints[0].order[0]
    assert path.phasepoints[20].order[0] == pathf.phasepoints[1].order[0]
    assert path.phasepoints[-1].order[0] == pathf.phasepoints[-1].order[0]


def test_check_crossing():
    """Test the check crossing method."""
    leftof, cross = check_crossing(0, 1.0, [-1.0, 0.0, 1.1], None)
    assert leftof == [False, False, True]
    assert not cross
    leftof, cross = check_crossing(1, 1.2, [-1.0, 0.0, 1.1], leftof)
    assert leftof == [False, False, False]
    assert cross[0] == (1, 2, '+')
    leftof, cross = check_crossing(10, -2, [-1.0, 0.0, 1.1], leftof)
    assert leftof == [True, True, True]
    assert cross[0] == (10, 0, '-')
    assert cross[1] == (10, 1, '-')
    assert cross[2] == (10, 2, '-')


def test_length():
    """Test for calculate the length of a path."""
    assert PATHTEST0.length == 64
    assert PATHTEST2.length == 12


def test_minvalue():
    """Test for calculate the min value of a path."""
    assert PATHTEST0.ordermin[0] == 0
    assert PATHTEST2.ordermin[0] == 4


def test_maxvalue():
    """Test for calculate the max value of a path."""
    assert PATHTEST0.ordermax[0] == 15
    assert PATHTEST2.ordermax[0] == 6


def test___eq__():
    """Test for equal properties for paths."""
    assert PATHTEST0 != PATHTEST2
    # Compare a path with something that is not a path:
    assert PATHTEST0 != 'This_is_not_a_path'
    # Compare where one path has "lost" one attribute:
    path2 = PATHTEST0.copy()
    assert PATHTEST0 == path2
    del path2.status
    assert PATHTEST0 != path2
    # Test when a phase point is different:
    path2 = PATHTEST0.copy()
    path2.phasepoints[-1].particles.pos = np.ones(3)
    assert PATHTEST0 != path2
    # Test with different time origin as an example of differing
    # attributes:
    path2 = PATHTEST0.copy()
    path2.time_origin = -1
    assert PATHTEST0 != path2


def test___ne__():
    """Test for non equal properties for paths."""
    assert PATHTEST0 != PATHTEST2


def test_path_reverse():
    """Test if we reverse correctly for class Path."""
    rgen = RandomGenerator(seed=0)
    path = Path(rgen)
    for _ in range(50):
        path.append(
            set_up_system([rgen.rand()[0]], np.ones(3),
                          np.ones(1), 0.0, 0.0)
        )
    path_rev = path.reverse(order_function=None)
    path_rev2 = path.reverse(order_function=SameOrder())
    assert path_rev == path_rev2
    for original, rev in zip(reversed(path.phasepoints),
                             path_rev.phasepoints):
        assert original.order[0] == pytest.approx(rev.order[0])
        assert np.allclose(original.particles.get_vel(),
                           rev.particles.get_vel() * -1)
        assert np.allclose(original.particles.get_pos(),
                           rev.particles.get_pos())
    # Test for cases when the order parameter depends on the
    # velocity:
    path_rev = path.reverse(NegativeOrder())
    for original, rev in zip(reversed(path.phasepoints),
                             path_rev.phasepoints):
        assert original.order[0] == pytest.approx(rev.order[0] * -1)
        assert np.allclose(original.particles.get_vel(),
                           rev.particles.get_vel() * -1)
        assert np.allclose(original.particles.get_pos(),
                           rev.particles.get_pos())


def test_path_exceed_maxlen():
    """Test that we stop adding points if we exceed the path max-length."""
    rgen = RandomGenerator(seed=0)
    path = Path(rgen, maxlen=10)
    for _ in range(path.maxlen):
        add = path.append(
            set_up_system([rgen.rand()[0]], np.zeros(3),
                          np.zeros(3), 0.0, 0.0)
        )
        assert add
    for _ in range(path.maxlen):
        add = path.append(
            set_up_system([rgen.rand()[0]], np.zeros(3),
                          np.zeros(3), 0.0, 0.0)
        )
        assert not add


def test_empty_path_creation():
    """Test that empty paths are created with correct type/settings."""
    rgen = RandomGenerator(seed=0)
    maxlen = 10
    path = Path(rgen, maxlen=maxlen)
    for _ in range(maxlen + 5):
        path.append(
            set_up_system([rgen.rand()[0]], np.zeros(3),
                          np.zeros(3), 0.0, 0.0)
        )

    path2 = path.empty_path(maxlen=maxlen)

    assert isinstance(path2, Path)
    assert path.maxlen == path2.maxlen
    assert path.rgen == path2.rgen


def test_get_min_max():
    """Test the getting of min/max order parameter."""
    rgen = RandomGenerator(seed=0)
    path = Path(rgen, maxlen=100)
    all_order = []
    for i in range(20):
        order = -1.0 * i
        if i == 10:
            order = 100.
        elif i == 15:
            order = -100.
        all_order.append(order)
        path.append(
            set_up_system([order], np.zeros(3), np.zeros(3), 0.0, 0.0)
        )
    ordermin, ordermax = path.ordermin, path.ordermax
    assert min(all_order) == pytest.approx(ordermin[0])
    assert max(all_order) == pytest.approx(ordermax[0])
    assert 15 == pytest.approx(ordermin[1])
    assert 10 == pytest.approx(ordermax[1])


def test_get_shooting_point_exp_clamps_upper_slab():
    """Points on the right interface must not overflow the slab index."""
    rgen = RandomGenerator(seed=0)
    path = Path(rgen, maxlen=100)
    for order in [0.0, 0.25, 0.5, 0.75, 1.0]:
        path.append(
            set_up_system([order], np.zeros(3), np.zeros(3), 0.0, 0.0)
        )

    phasepoint, idx = path.get_shooting_point(
        criteria='exp', interfaces=[0.0, 0.5, 1.0]
    )

    assert phasepoint is path.phasepoints[idx]
    assert 0 <= idx < path.length


def test_check_interfaces():
    """Test the check interfaces method."""
    path = Path(None, maxlen=100)
    ret = path.check_interfaces([1.0, 4.0, 5.0])
    assert all((i is None for i in ret))
    for i in range(5):
        path.append(
            set_up_system([i], np.zeros(3), np.zeros(3), 0.0, 0.0)
        )
    ret = path.check_interfaces([1.0, 3.0, 5.0])
    assert ret[0] == 'L'
    assert ret[1] is None
    assert ret[2] == 'M'
    assert ret[3] == [True, True, False]


def test_start_end():
    """Test the get start/end points method."""
    path = Path(None, maxlen=100)
    for i in range(5):
        path.append(
            set_up_system([i * -1.0], np.zeros(3), np.zeros(3), 0.0, 0.0)
        )
    assert path.get_start_point(0) == 'L'
    assert path.get_end_point(1) == 'L'
    assert path.get_end_point(0, 1) == 'L'
    assert path.get_end_point(-10, -6) == 'R'
    assert path.get_end_point(-100, -1) is None
    assert path.get_start_point(0) == 'L'
    assert path.get_start_point(0, 1) == 'L'
    assert path.get_start_point(-4, -3) == 'R'
    assert path.get_start_point(-2, 1) is None


def test_get_path_data():
    """Test the get_path_data and set_move methods."""
    path = Path(None, maxlen=100)
    path.set_move('fake')
    for i in range(5):
        path.append(
            set_up_system([i], np.zeros(3), np.zeros(3), 0.0, 0.0)
        )
    path_info = path.get_path_data('ACC', [1.0, 2.0, 3.0])
    correct = {'length': 5, 'ordermax': (4, 4), 'ordermin': (0, 0),
               'generated': ('fake', 0, 0, 0),
               'status': 'ACC',
               'interface': ('L', 'M', 'R')}
    for key, val in correct.items():
        assert val == path_info[key]


def test_success():
    """Test the success method."""
    path = Path(None, maxlen=100)
    for i in range(5):
        path.append(
            set_up_system([i], np.zeros(3), np.zeros(3), 0.0, 0.0)
        )
    assert path.success(3.0)
    assert not path.success(4.0)


def test_update_energies():
    """Test the update energies method."""
    path = Path(None, maxlen=100)
    for i in range(5):
        path.append(
            set_up_system([i], np.zeros(3), np.zeros(3), 1.0, 2.0)
        )
    # Test if ekin and vpot have correct length:
    ekin = [10] * path.length
    vpot = [10] * path.length
    path.update_energies(ekin, vpot)
    ekin2 = [i.particles.ekin for i in path.phasepoints]
    vpot2 = [i.particles.vpot for i in path.phasepoints]
    assert ekin == ekin2
    assert vpot == vpot2
    # Test for one is longer than the other:
    vpot = [11] * path.length + [12]
    path.update_energies(ekin, vpot)
    vpot2 = [i.particles.vpot for i in path.phasepoints]
    assert ekin == ekin2
    assert vpot[:path.length] == vpot2
    # Test when energies are too short:
    ekin = [10] * (path.length - 2)
    vpot = [10] * (path.length - 1)
    path.update_energies(ekin, vpot)
    ekin2 = [i.particles.ekin for i in path.phasepoints]
    vpot2 = [i.particles.vpot for i in path.phasepoints]
    assert ekin + [None, None] == ekin2
    assert vpot + [None] == vpot2


def test_add():
    """Test the __iadd__ method."""
    path = Path(None, maxlen=10)
    path2 = Path(None, maxlen=10)
    for i in range(5):
        path.append(
            set_up_system([i], np.zeros(3), np.zeros(3), 0.0, 0.0)
        )
        path2.append(
            set_up_system([i * 10], np.zeros(3), np.zeros(3), 0.0, 0.0)
        )
    path += path2
    assert path.length == 10
    for i, phasepoint in enumerate(path.phasepoints):
        if i <= 4:
            assert phasepoint.order[0] == i
        else:
            assert phasepoint.order[0] == (i - 5) * 10
    # Try to add some more points (we are now > maxlen):
    path += path2
    assert path.length == 10


def test_delete():
    """Test the delete method."""
    path = Path(None, maxlen=10)
    for i in range(5):
        path.append(
            set_up_system([i], np.zeros(3), np.zeros(3), 0.0, 0.0)
        )
    assert path.length == 5
    for i in range(path.length):
        path.delete(0)
    # Check that all points were deleted:
    assert path.length == 0


def test_sorting():
    """Test the sorting method."""
    data = [
        ([1, 300, 1, 2, -123], ('a', 'b'), ('a', 'b'), 1, 1),
        ([3, -4], ('e', 'f'), ('e', 'f'), 5, 4),
        ([2, 0, 0], ('c', 'd'), ('c', 'd'), 3, 2),
        ([5], ('i', 'j'), ('i', 'j'), 8, 7),
        ([4, 700000], ('g', 'h'), ('g', 'h'), 8, 7),
    ]
    path = Path(None, maxlen=10)
    for datai in data:
        path.append(
            set_up_system(datai[0], datai[1], datai[2],
                          datai[3], datai[4], internal=False)
        )
    correct_order = [1, 2, 3, 4, 5]
    correct_ekin = [1, 2, 4, 7, 7]
    correct_vpot = [1, 3, 5, 8, 8]
    # Sort path by order:
    path.sorting('order', reverse=False)
    sort = [i.order[0] for i in path.phasepoints]
    assert sort == correct_order
    # Sort path by order, reversed:
    path.sorting('order', reverse=True)
    sort = [i.order[0] for i in path.phasepoints]
    assert sort[::-1] == correct_order
    # Sort path by kinetic energy:
    path.sorting('ekin', reverse=False)
    sort = [i.particles.ekin for i in path.phasepoints]
    assert sort == correct_ekin
    # Sort path by kinetic energy, reversed:
    path.sorting('ekin', reverse=True)
    sort = [i.particles.ekin for i in path.phasepoints]
    assert sort[::-1] == correct_ekin
    # Sort path by potential energy:
    path.sorting('vpot', reverse=False)
    sort = [i.particles.vpot for i in path.phasepoints]
    assert sort == correct_vpot
    # Sort path by potential energy, reversed:
    path.sorting('vpot', reverse=True)
    sort = [i.particles.vpot for i in path.phasepoints]
    assert sort[::-1] == correct_vpot
    with pytest.raises(AttributeError):
        path.sorting('barbara')
    with pytest.raises(AttributeError):
        path.sorting('gigio')


def test_external_reverse():
    """Test that we can reverse external paths."""
    path = Path(None)
    path.append(
        set_up_system([0.0, None], ('initial.g96', None),
                      False, None, None, internal=False)
    )
    for i in range(5):
        path.append(
            set_up_system([(i + 1) * 10, None], ('trajB.trr', i),
                          True, None, None, internal=False)
        )
    for i in range(5):
        path.append(
            set_up_system([(i + 1) * 20, None], ('trajF.trr', i),
                          False, None, None, internal=False)
        )
    rev = path.reverse(SameOrder())
    for point1, point2 in zip(rev.phasepoints,
                              reversed(path.phasepoints)):
        assert point1.particles.vel_rev != point2.particles.vel_rev
        assert point1.particles.get_pos() == point2.particles.get_pos()
        assert point1.order[0] == point2.order[0]


def test_external_restart():
    """Tests the restart feature for external paths."""
    path = Path(None)
    path.append(
        set_up_system([0.0, None], ('initial.g96', None),
                      False, None, None, internal=False)
    )
    for i in range(5):
        path.append(
            set_up_system([(i + 1) * 10, None], ('trajB.trr', i),
                          True, None, None, internal=False)
        )
    for i in range(5):
        path.append(
            set_up_system([(i + 1) * 20, None], ('trajF.trr', i),
                          False, None, None, internal=False)
        )
    dic = {
        'generated': ('fake', 0, 0, 0),
        'status': 'ACC',
        'weight': 5,
    }
    for attr, value in dic.items():
        setattr(path, attr, value)
    info = path.restart_info()
    path2 = path.empty_path()
    path2.load_restart_info(info)

    for i, j in enumerate(path.rgen.get_state()['state']):
        if isinstance(j, np.ndarray):
            path2_l = path2.rgen.get_state()['state'][i].tolist()
            assert j.tolist() == path2_l
        else:
            assert j == path2.rgen.get_state()['state'][i]
    assert path.rgen.get_state()['seed'] == path2.rgen.get_state()['seed']
    assert path.generated == path2.generated
    assert path.maxlen == path2.maxlen
    assert path.time_origin == path2.time_origin
    assert path.status == path2.status
    assert path.weight == path2.weight
    for i, j in enumerate(path.phasepoints):
        assert j.order[0] == path2.phasepoints[i].order[0]
