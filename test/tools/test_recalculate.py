# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the recalculate orderparameter tools in pyretis.tools"""
from contextlib import contextmanager
from io import StringIO
import logging
import os
import tarfile
from tempfile import NamedTemporaryFile
import pytest
from unittest.mock import patch
import numpy as np
from pyretis.tools.recalculate_order import (
    recalculate_from_frame,
    recalculate_from_trj,
    recalculate_order,
    get_traj_files,
)

from pyretis.orderparameter import OrderParameter


logging.disable(logging.CRITICAL)


@contextmanager
def turn_on_logging():
    """Turn on logging so that tests can detect it."""
    logging.disable(logging.INFO)
    try:
        yield
    finally:
        logging.disable(logging.CRITICAL)


CORRECT_TRI = [
    [1.2000000476837158, 1.100000023841858, 0.9526299834251404,
     0.0, 0.0, 0.0, 0.0, 0.550000011920929, 0.0],
    [1.3190346956253052, 1.1973484754562378, 1.0101481676101685,
     0.0, 0.0, 0.0, 0.0, 0.6045578122138977, 0.0],
    [1.445665717124939, 1.2725739479064941, 1.043269157409668,
     0.0, 0.0, 0.0, 0.0, 0.6625972390174866, 0.0],
    [1.549871802330017, 1.3433609008789062, 1.0741169452667236,
     0.0, 0.0, 0.0, 0.0, 0.7103575468063354, 0.0],
    [1.6492804288864136, 1.3989951610565186, 1.0993438959121704,
     0.0, 0.0, 0.0, 0.0, 0.7559200525283813, 0.0],
    [1.7235912084579468, 1.4501303434371948, 1.1253700256347656,
     0.0, 0.0, 0.0, 0.0, 0.7899801731109619, 0.0],
    [1.8234614133834839, 1.5086286067962646, 1.1493875980377197,
     0.0, 0.0, 0.0, 0.0, 0.8357531428337097, 0.0],
    [1.9387280941009521, 1.5735180377960205, 1.1708111763000488,
     0.0, 0.0, 0.0, 0.0, 0.888584315776825, 0.0],
    [2.0414657592773438, 1.6291977167129517, 1.1705245971679688,
     0.0, 0.0, 0.0, 0.0, 0.9356740117073059, 0.0],
    [2.138443946838379, 1.6822302341461182, 1.191916584968567,
     0.0, 0.0, 0.0, 0.0, 0.9801228642463684, 0.0],
    [2.234323024749756, 1.7414289712905884, 1.219726800918579,
     0.0, 0.0, 0.0, 0.0, 1.0240683555603027, 0.0],
]


CORRECT_XYZ = [
    [-1.0912E-05, 1.85E-07, -2.310E-06],
    [-2.1178E-05, 6.2E-08, -4.785E-06],
    [-3.1073E-05, -3.8E-07, -7.123E-06],
    [-4.1081E-05, -1.026E-06, -9.074E-06],
    [-5.1297E-05, -1.657E-06, -1.0489E-05],
    [-6.1028E-05, -2.03E-06, -1.1326E-05],
]


HERE = os.path.abspath(os.path.dirname(__file__))


class OrderBox(OrderParameter):
    """Order parameter returning the box cell as a numpy array.

    Intentionally returns a raw numpy array to exercise the
    ``wrap_orderparameter`` normalisation path.
    """

    def __init__(self):
        super().__init__(description='Box order parameter')

    def calculate(self, system):
        return system.box.cell


class Xpos1partExt(OrderParameter):
    """Order parameter returning the X position of particle 0 as a scalar.

    Intentionally returns a raw scalar to exercise the
    ``wrap_orderparameter`` normalisation path.
    """

    def __init__(self):
        super().__init__(description='X pos particle 1')

    def calculate(self, system):
        return system.particles.pos[0][0]


class Velocity(OrderParameter):
    """Order parameter returning the per-axis velocity sum as a numpy array.

    Intentionally returns a raw numpy array to exercise the
    ``wrap_orderparameter`` normalisation path.
    """

    def __init__(self):
        super().__init__(description='Box order parameter')

    def calculate(self, system):
        return system.particles.vel.sum(axis=0)


class TestRecalculateOrder:
    """Test the calculation of order parameters from trajectory files."""

    def test_rectangular_box(self):
        """Test recalculate for a rectangular box."""
        orderf = OrderBox()
        box = np.array([2.7200000286102295]*3)
        trr = os.path.join(HERE, '4water.trr')
        with patch('sys.stdout', new=StringIO()):
            for i in recalculate_order(orderf, trr, {}):
                assert np.allclose(i, box)

    def test_novel(self, caplog):
        """Test recalculate for a rectangular box."""
        trr = os.path.join(HERE, 'traj-novel.trr')
        with turn_on_logging(), patch('sys.stdout', new=StringIO()):
            with caplog.at_level(logging.WARNING,
                                 logger='pyretis.tools.recalculate_order'):

                for order in recalculate_order(Velocity(), trr, {}):
                    for i in order:
                        assert i == pytest.approx(0.0)
                assert any("No velocities found in" in rec.message
                           for rec in caplog.records)

    def test_tri_box(self):
        """Test recalculate for a triclinic box."""
        orderf = OrderBox()
        trr = os.path.join(HERE, '4water_tri_dyn.trr')
        orderp = []
        with patch('sys.stdout', new=StringIO()):
            for i, j in zip(recalculate_order(orderf, trr, {}), CORRECT_TRI):
                orderp.append(i)
                assert np.allclose(i, j)
            assert len(orderp) == len(CORRECT_TRI)
            for i, j in zip(
                    recalculate_order(orderf, trr, {'reverse': True}),
                    orderp[::-1]):
                assert np.allclose(i, j)

    def test_sub_selection(self):
        """Test when we specify a min/max index to consider."""
        orderf = OrderBox()
        trr = os.path.join(HERE, '4water_tri_dyn.trr')
        idx = 3
        with patch('sys.stdout', new=StringIO()):
            orderp = [i for i in recalculate_order(orderf, trr,
                                                   {'maxidx': idx})]
            for i, j in zip(orderp, CORRECT_TRI[:idx+1]):
                assert np.allclose(i, j)
        with patch('sys.stdout', new=StringIO()):
            orderp = [i for i in recalculate_order(orderf, trr,
                                                   {'minidx': idx})]
            for i, j in zip(orderp, CORRECT_TRI[idx:]):
                assert np.allclose(i, j)

    def test_recalculate_gro(self):
        """Recalculate from .gro file."""
        orderf = OrderBox()
        tar = tarfile.open(
            os.path.join(HERE, '4water_tri_dyn.gro.tgz'), 'r:gz'
        )
        members = sorted([i.name for i in tar.getmembers()])
        i = 0
        with patch('sys.stdout', new=StringIO()):
            for member in members:
                gro = tar.extractfile(member)
                if gro is not None:
                    grofile = NamedTemporaryFile(suffix='.gro')
                    with open(grofile.name, 'wb') as output:
                        output.write(gro.read())
                    box = recalculate_order(orderf, grofile.name, {})
                    for j, k in zip(box[0], CORRECT_TRI[i]):
                        assert j == pytest.approx(k, abs=10**-(5))
                    i += 1
        # Do some additional tests:
        with patch('sys.stdout', new=StringIO()):
            gro = tar.extractfile(members[1])
            grofile = NamedTemporaryFile(suffix='.gro')
            with open(grofile.name, 'wb') as output:
                output.write(gro.read())
            vel1, = recalculate_from_frame(Velocity(), grofile.name,
                                           {'ext': '.gro', 'reverse': False})
            vel2, = recalculate_from_frame(Velocity(), grofile.name,
                                           {'ext': '.gro', 'reverse': True})
            # vel1 and vel2 are normalised lists; use np.negative for negation
            assert np.allclose(vel1, np.negative(vel2))
        with patch('sys.stdout', new=StringIO()):
            with pytest.raises(ValueError):
                recalculate_from_frame(Velocity(), 'something.txt',
                                       {'ext': '.txt'})

    def test_recalculate_g96(self):
        """Recalculate from .g96."""
        orderf = OrderBox()
        tar = tarfile.open(
            os.path.join(HERE, '4water_tri_dyn.g96.tgz'), 'r:gz'
        )
        members = sorted([i.name for i in tar.getmembers()])
        i = 0
        with patch('sys.stdout', new=StringIO()):
            for member in members:
                g96 = tar.extractfile(member)
                if g96 is not None:
                    g96file = NamedTemporaryFile(suffix='.g96')
                    with open(g96file.name, 'wb') as output:
                        output.write(g96.read())
                    box = recalculate_order(orderf, g96file.name, {})
                    for j, k in zip(box[0], CORRECT_TRI[i]):
                        assert j == pytest.approx(k, abs=10**-(5))
                    i += 1

    def test_recalculate_xyz(self):
        """Recalculate from .xyz."""
        orderf = Velocity()
        filename = os.path.join(HERE, 'traj.xyz')
        with patch('sys.stdout', new=StringIO()):
            order = [i for i in recalculate_order(orderf, filename, {})]
            order_r = [
                i for i in recalculate_order(
                    orderf, filename, {'reverse': True}
                )
            ]
            assert len(order) == len(CORRECT_XYZ)
            assert len(order) == len(order_r)
            for vel1, vel2, vel3 in zip(order, order_r[::-1], CORRECT_XYZ):
                for i, j, k in zip(vel1, vel2, vel3):
                    assert i == pytest.approx(-j, abs=10**-(7))
                    assert i == pytest.approx(k, abs=10**-(7))
        # Test with subsection:
        idx = 2
        with patch('sys.stdout', new=StringIO()):
            order = [
                i for i in recalculate_order(orderf, filename, {'maxidx': idx})
            ]
            for i, j in zip(order, CORRECT_XYZ[:idx+1]):
                assert np.allclose(i, j)
        with patch('sys.stdout', new=StringIO()):
            order = [
                i for i in recalculate_order(orderf, filename, {'minidx': idx})
            ]
            for i, j in zip(order, CORRECT_XYZ[idx:]):
                assert np.allclose(i, j)
        # Test with missing box:
        filename = os.path.join(HERE, 'traj-no-box.xyz')
        orderf = OrderBox()
        with patch('sys.stdout', new=StringIO()):
            order = [
                i for i in recalculate_order(orderf, filename,
                                             {'box': [6., 6., 6.]})
            ]
            for orderi in order:
                assert len(orderi) == 3
                for i in orderi:
                    assert i == pytest.approx(6.)

    def test_get_traj_files(self):
        """Test the method for getting traj files."""
        traj = os.path.join(HERE, 'traj-ext.txt')
        correct = {
            'trajB.trr': {'minidx': 0, 'maxidx': 6, 'reverse': True},
            'trajF.trr': {'minidx': 0, 'maxidx': 7, 'reverse': False},
        }
        files = get_traj_files(traj)
        for key, val in correct.items():
            assert key in files
            assert val == files[key]

    def test_use_exteral(self):
        """Test the method for load and recalculate op from ext trjs."""
        orderf = Xpos1partExt()
        sim_file = os.path.join(HERE, '2water')

        ordercheck = [0.1259, 0.1257, 0.1255, 0.1253, 0.1251, 0.1249]

        with patch('sys.stdout', new=StringIO()):
            order = [i for i in recalculate_order(
                order_parameter=orderf,
                traj_file=sim_file+'.trr',
                options={'top': sim_file+'.gro'})]

        for orderi, orderj in zip(order, ordercheck):
            assert orderi[0] == pytest.approx(orderj, 3, 0.1)

        with patch('sys.stdout', new=StringIO()):
            order = [i for i in recalculate_from_trj(
                orderf, sim_file + '.trr',
                {'top': sim_file + '.gro'})]

        for orderi, orderj in zip(order, ordercheck):
            assert orderi[0] == pytest.approx(orderj, 3, 0.1)

        one_op = recalculate_from_trj(orderf, sim_file + '.trr',
                                      {'top': sim_file + '.gro',
                                       'idx': 3})

        with patch('sys.stdout', new=StringIO()):
            assert next(one_op)[0] == pytest.approx(ordercheck[3], 3, 0.1)
