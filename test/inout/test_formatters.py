# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A test module for the formatters."""
import logging
import os
import pytest
import numpy as np
from pyretis.core import create_box, Particles, System
from pyretis.core.path import Path
from pyretis.core.pathensemble import PathEnsemble
from pyretis.core.random_gen import MockRandomGenerator
from pyretis.inout.formats.formatter import (
    OutputFormatter,
    apply_format,
    format_number,
    get_log_formatter,
    PyretisLogFormatter,
    LOG_FMT,
    LOG_DEBUG_FMT,
)
from pyretis.inout.formats.cross import CrossFormatter
from pyretis.inout.formats.energy import EnergyFormatter, EnergyPathFormatter
from pyretis.inout.formats.order import OrderFormatter, OrderPathFormatter
from pyretis.inout.formats.path import PathExtFormatter, PathIntFormatter
from pyretis.inout.formats.snapshot import (
    SnapshotFormatter,
    adjust_coordinate,
)
from pyretis.inout.formats.pathensemble import PathEnsembleFormatter
from .help import (
    assert_equal_path_dict,
    create_external_path,
    turn_on_logging,
    set_up_system,
    create_test_paths,
    DATA_CROSS,
    DATA_CROSS2,
    CORRECT_CROSS,
    DATA_ENERGY,
    CORRECT_ENERGY,
    DATA_ORDER,
    CORRECT_ORDER,
    CORRECT_PATH_EXT,
    PATH_EXT_RAW,
)


logging.disable(logging.CRITICAL)
HERE = os.path.abspath(os.path.dirname(__file__))


class TestOutputFormatter:
    """Test the generic output formatter."""

    def test_initiation(self, caplog):
        """Test initiation and the header."""
        formatter = OutputFormatter('test-formatter', header=None)
        assert formatter.header is None

        with pytest.raises(AttributeError):
            OutputFormatter('test-file', header='Just a text header')

        txt = 'This is the header, it is very long indeed'
        formatter = OutputFormatter('test-file', header={'text': txt})
        assert formatter.header == txt

        formatter = OutputFormatter('test-file', header={'text': txt,
                                                         'width': 10})
        assert formatter.header == txt

        formatter = OutputFormatter(
            'test-file', header={'width': (10, 10),
                                 'labels': ('Label1', 'Label2')},
        )
        assert formatter.header == '#   Label1     Label2'
        formatter.header = txt
        assert formatter.header == txt

    def test_parse(self, caplog):
        """Test that the generic parser work for simple data."""
        formatter = OutputFormatter('test-formatter', header=None)
        raw = '1 2.0 3.0     4.123'
        data = formatter.parse(raw)
        correct = [1.0, 2.0, 3.0, 4.123]
        assert len(data) == len(correct)
        for i, j in zip(data, correct):
            assert i == pytest.approx(j)


class TestCrossFormatter:
    """Test the CrossFormatter formatter."""

    def test_cross_formatter(self, caplog):
        """Test that we can format crossing data."""
        formatter = CrossFormatter()
        assert formatter.header == '#     Step  Int Dir'
        for i, datai in enumerate(DATA_CROSS):
            for line in formatter.format(i, datai):
                assert line == CORRECT_CROSS[i]

    def test_parse(self, caplog):
        """Test that we can parse crossing data."""
        formatter = CrossFormatter()
        for i, j in zip(CORRECT_CROSS, DATA_CROSS2):
            data = formatter.parse(i)
            assert data == j


class TestEnergyFormatter:
    """Test the EnergFMT formatter."""

    def test_energy_file_formatter(self, caplog):
        """Test that we can format energy data."""
        formatter = EnergyFormatter()
        fields = ['vpot', 'ekin', 'etot', 'temp']
        assert formatter.header == (
            '#     Time      Potential        Kinetic'
            '          Total    Temperature')
        for step, i in enumerate(DATA_ENERGY):
            data = {key: j for (key, j) in zip(fields, i)}
            for line in formatter.format(step, data):
                assert line == CORRECT_ENERGY[step]

    def test_energy_missing_data(self, caplog):
        """Test the energy formatting when some data is missing."""
        formatter = EnergyFormatter()
        data = {'vpot': 1.0, 'ekin': 2.0, 'etot': 3.0}
        for line in formatter.format(0, data):
            assert line == ('         0       1.000000       2.000000'
                            '       3.000000            nan')


class TestEnergyPathFormatter:
    """Test the EnergPathFMT formatter."""

    def test_energy_path_formatter(self, caplog):
        """Test that we can format energy data from paths."""
        formatter = EnergyPathFormatter()
        paths, _ = create_test_paths(npath=1)
        paths[0].generated = 'Cha-Cha'
        correct = ['# Cycle: 0, status: ACC, move: Cha-Cha',
                   '#     Time      Potential        Kinetic']
        for i, line in enumerate(formatter.format(0, [paths[0], 'ACC'])):
            if i < 2:
                assert line == correct[i]
            else:
                point = paths[0].phasepoints[i-2]
                txt = '{:>10d} {:>14.6f} {:>14.6f}'.format(
                    i-2, point.particles.vpot, point.particles.ekin
                )
                assert line == txt
        # Test handling of empty paths:
        with pytest.raises(StopIteration):
            next(formatter.format(0, [None, 'ACC']))
        # Just check what happens if we are missing one of the attributes:
        del paths[0].phasepoints[0].particles.vpot
        for i, line in enumerate(formatter.format(0, [paths[0], 'ACC'])):
            if i == 2:
                split = line.split()
                assert 'nan' == split[1]


class TestOrderFormatter:
    """Test the OrderFormatter formatter."""

    def test_order_parameter_formatter(self, caplog):
        """Test that we can format order parameter data."""
        formatter = OrderFormatter()
        assert formatter.header == '#     Time       Orderp'
        for i, datai in enumerate(DATA_ORDER):
            for line in formatter.format(i, datai):
                assert line == CORRECT_ORDER[i]


class TestOrderPathFormatter:
    """Test the OrderPathFormatter formatter."""

    def test_order_path_formatter(self, caplog):
        """Test that we can format order parameter data for paths."""
        formatter = OrderPathFormatter()
        paths, _ = create_test_paths(npath=1)
        paths[0].generated = 'ld'
        correct = ['# Cycle: 0, status: ACC, move: ld',
                   '#     Time       Orderp']
        for i, line in enumerate(formatter.format(0, [paths[0], 'ACC'])):
            if i < 2:
                assert line == correct[i]
            else:
                point = paths[0].phasepoints[i-2]
                fmt = ' '.join(['{:>10d}'] + len(point.order)*['{:>12.6f}'])
                txt = fmt.format(i-2, *point.order)
                assert line == txt
        # Test handling of empty paths:
        with pytest.raises(StopIteration):
            next(formatter.format(0, [None, 'ACC']))


def generate_system_snapshots(ndim, length=10, npart=10):
    """Generate a set of system snapshots.

    Parameters
    ----------
    ndim : int
        The number of dimensions to consider.
    length : int, optional
        The number of snapshots to generate.
    npart : int, optional
        The number of particles to add.

    """
    box = create_box(low=[0.]*ndim, high=[10.]*ndim)
    system = System(box=box)
    system.particles = Particles(dim=ndim)
    for i in range(npart):
        name = 'Kr'
        if i % 2 == 0:
            name = 'Ar'
        system.add_particle(np.zeros(ndim), mass=1.0, name=name, ptype=0)
    for i in range(length):
        for j in range(npart):
            for k in range(ndim):
                system.particles.pos[j, k] = (i + 1) * (j + 1) * (k + 1)
                system.particles.vel[j, k] = (i + 0.5) * (j + 0.5) * (k + 0.5)
                if i % 2 == 0:
                    system.particles.vel *= -1
        yield system


def compare_traj_fmt_file(formatter, correct_file, ndim):
    """Compare formatted trajectory output with a file.

    Parameters
    ----------
    formatter : object like :py:class:`.SnapshotFormatter`
        Formatter used to format the trajectory/systems.
    correct_file : string
        Path to the file with the correctly formatted output.
    ndim : int
        The number of dimensions to use (1, 2 or 3).

    Returns
    -------
    out : boolean
        True, if the data is correctly formatted.

    """
    correct = []
    with open(correct_file, 'r', encoding='utf-8') as infile:
        for line in infile:
            correct.append(line.rstrip())
    j = 0
    for i, system in enumerate(generate_system_snapshots(ndim,
                                                         length=5)):
        for line in formatter.format(i, system):
            assert line == correct[j]
            j += 1


class TestSnapshotFormatter:
    """Test the SnapshotFormatter formatter."""

    def test_traj_formatter_full(self, caplog):
        """Test that we can format trajectories."""
        formatter1 = SnapshotFormatter(write_vel=False, fmt='full')
        formatter2 = SnapshotFormatter(write_vel=True, fmt='full')
        cases = [
            (1, 'trajfmt_full_novel_1D.txt', formatter1),
            (2, 'trajfmt_full_novel_2D.txt', formatter1),
            (3, 'trajfmt_full_novel_3D.txt', formatter1),
            (1, 'trajfmt_full_vel_1D.txt', formatter2),
            (2, 'trajfmt_full_vel_2D.txt', formatter2),
            (3, 'trajfmt_full_vel_3D.txt', formatter2),
        ]
        for (ndim, filename, formatter) in cases:
            correct_file = os.path.join(HERE, filename)
            compare_traj_fmt_file(formatter, correct_file, ndim)

    def test_traj_formatter(self, caplog):
        """Test that we can format trajectories."""
        formatter1 = SnapshotFormatter(write_vel=False)
        with turn_on_logging():
            with caplog.at_level(logging.WARNING,
                                 logger='pyretis.inout.formats.snapshot'):
                formatter1.parse('some text')
        formatter2 = SnapshotFormatter(write_vel=True)
        cases = [
            (1, 'trajfmt_novel_1D.txt', formatter1),
            (2, 'trajfmt_novel_2D.txt', formatter1),
            (3, 'trajfmt_novel_3D.txt', formatter1),
            (1, 'trajfmt_vel_1D.txt', formatter2),
            (2, 'trajfmt_vel_2D.txt', formatter2),
            (3, 'trajfmt_vel_3D.txt', formatter2),
        ]
        for (ndim, filename, formatter) in cases:
            correct_file = os.path.join(HERE, filename)
            compare_traj_fmt_file(formatter, correct_file, ndim)

    def test_adjust_coordinate(self, caplog):
        """Test that we can adjust coordinates."""
        # 1D:
        coord = np.array([10, ])
        coord_ = adjust_coordinate(coord)
        assert np.allclose(coord_, [[10.0, 0.0, 0.0]])
        # 1 particle, 1D:
        particles = Particles(dim=1)
        particles.add_particle(np.array([1.0]),
                               np.zeros(1),
                               np.zeros(1))
        pos = adjust_coordinate(particles.pos)
        assert np.allclose(pos, np.array([1.0, 0.0, 0.0]))
        # 1 particle, 2D:
        particles = Particles(dim=2)
        particles.add_particle(np.array([1.0, 1.0]),
                               np.zeros(2),
                               np.zeros(2))
        pos = adjust_coordinate(particles.pos)
        assert np.allclose(pos, np.array([1.0, 1.0, 0.0]))
        # 1 particle, 3D:
        particles = Particles(dim=3)
        particles.add_particle(np.array([1.0, 1.0, 1.0]),
                               np.zeros(3),
                               np.zeros(3))
        pos = adjust_coordinate(particles.pos)
        assert np.allclose(pos, np.array([1.0, 1.0, 1.0]))
        # 2 particles, 1D:
        particles = Particles(dim=1)
        particles.add_particle(np.array([1.0]),
                               np.zeros(1),
                               np.zeros(1))
        particles.add_particle(np.array([-1.0]),
                               np.zeros(1),
                               np.zeros(1))
        pos = adjust_coordinate(particles.pos)
        assert np.allclose(pos, np.array([[1.0, 0.0, 0.0],
                                          [-1.0, 0.0, 0.0]]))
        # 2 particles, 2D:
        particles = Particles(dim=2)
        particles.add_particle(np.array([1.0, -1.0]),
                               np.zeros(2),
                               np.zeros(2))
        particles.add_particle(np.array([-1.0, 1.0]),
                               np.zeros(2),
                               np.zeros(2))
        pos = adjust_coordinate(particles.pos)
        assert np.allclose(pos, np.array([[1.0, -1.0, 0.0],
                                          [-1.0, 1.0, 0.0]]))

        # 3 particles, 3D:
        particles = Particles(dim=3)
        particles.add_particle(np.array([1.0, -1.0, 0.5]),
                               np.zeros(3),
                               np.zeros(3))
        particles.add_particle(np.array([-1.0, 1.0, -0.5]),
                               np.zeros(3),
                               np.zeros(3))
        pos = adjust_coordinate(particles.pos)
        assert np.allclose(pos, np.array([[1.0, -1.0, 0.5],
                                          [-1.0, 1.0, -0.5]]))


class TestPathIntFormatter:
    """Test that the PathFMT formatter work as intended."""

    def test_pathint_formatter(self, caplog):
        """Test the formatter for internal paths."""
        formatter = PathIntFormatter()
        paths, _ = create_test_paths(npath=2)
        k = -1
        for i, line in enumerate(formatter.format(0, [paths[1], 'ACC'])):
            if i == 0:
                assert '# Cycle: 0, status: ACC' == line
            else:
                j = i - 1
                if j % 10 == 0:
                    k += 1
                    assert f'Snapshot: {k}' == line
                else:
                    # Assume 3D:
                    k_float = float(k)
                    pos = [k_float for _ in range(3)]
                    vel = [-k_float for _ in range(3)]
                    fmt = ' '.join(['{:15.9f}'] * 6)
                    txt = fmt.format(*pos, *vel)
                    assert txt == line
        # Test handling of empty paths:
        with pytest.raises(StopIteration):
            next(formatter.format(0, [None, 'ACC']))

    def test_parse(self, caplog):
        """Test that we can parse trajectory data."""
        formatter = PathIntFormatter()
        raw = '1.0 2.0 3.0 4.0 5.0 6.0'
        pos, vel = formatter.parse(raw)
        correct_pos = np.array([1., 2., 3.])
        correct_vel = np.array([4., 5., 6.])
        assert np.allclose(pos, correct_pos)
        assert np.allclose(vel, correct_vel)
        with pytest.raises(ValueError):
            raw = '1.0'
            formatter.parse(raw)
        with pytest.raises(ValueError):
            raw = '1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0'
            formatter.parse(raw)


class TestPathExtFormatter:
    """Test that the external path formatter work as intended."""

    def test_pathext_formatter(self, caplog):
        """Test the formatter for external paths."""
        path, _ = create_external_path()
        formatter = PathExtFormatter()
        for corr, snap in zip(CORRECT_PATH_EXT,
                              formatter.format(0, [path, 'ACC'])):
            assert corr == snap
        # Test handling of empty paths:
        with pytest.raises(StopIteration):
            next(formatter.format(0, [None, 'ACC']))

    def test_parse(self, caplog):
        """Test parsing of path data."""
        formatter = PathExtFormatter()
        k = 0
        for line in CORRECT_PATH_EXT:
            if line.strip().startswith('#'):
                continue
            data = formatter.parse(line)
            assert list(data) == PATH_EXT_RAW[k]
            k += 1


class TestPathEnsembleFormatter:
    """Test the formatting of path ensemble data."""

    def test_formatting(self, caplog):
        """Test the formatting of path data."""
        formatter = PathEnsembleFormatter()
        rgen = MockRandomGenerator(seed=0)
        path = Path(rgen)
        for i in range(10):
            path.append(
                set_up_system([i], np.ones((5, 3)) * i, np.ones((5, 3)) * i,
                              vpot=0.0, ekin=0.0, internal=True)
            )
        path.generated = ('sh', 0, 0, 0)
        ens = PathEnsemble(1, [0.5, 5.0, 8.0])
        ens.add_path_data(path, 'ACC')
        corrects = [
            ('         0          1          1 L M R      10 ACC sh'
             '  0.000000000e+00  9.000000000e+00       0       9'
             '  0.000000000e+00       0       0  1.000000000e+00'),
        ]
        for line, correct in zip(formatter.format(0, ens), corrects):
            assert line == correct

    def test_parse(self, caplog):
        """Test parsing of formatted data."""
        formatter = PathEnsembleFormatter()
        raw = [
            ('         2          1          0 L * L    1019 '
             'NCR sh -9.020623174e-01 -3.052974085e-01    1018'
             '     554 -4.283683762e-01     721     292  2.000000000e+00'),
            ('         7          4          1 L M R    1140 '
             'ACC sh -9.018613917e-01  1.001603711e+00       0'
             '    1139 -2.765514723e-01     498     322  2.000000000e+00'),
        ]
        correct = [
            {
                'cycle': 2, 'generated': ['sh', -0.4283683762, 721, 292],
                'interface': ('L', '*', 'L'), 'length': 1019,
                'ordermax': (-0.3052974085, 554),
                'ordermin': (-0.9020623174, 1018), 'status': 'NCR',
                'weight': 2.
            },
            {
                'cycle': 7, 'generated': ['sh', -0.2765514723, 498, 322],
                'interface': ('L', 'M', 'R'), 'length': 1140,
                'ordermax': (1.001603711, 1139),
                'ordermin': (-0.9018613917, 0),
                'status': 'ACC',
                'weight': 2.
            },
        ]

        for i, j in zip(raw, correct):
            assert_equal_path_dict(formatter.parse(i), j)
        with turn_on_logging():
            with caplog.at_level(
                    logging.WARNING,
                    logger='pyretis.inout.formats.pathensemble'):
                formatter.parse('1 2 3')
        # Also test a line with a comment:
        out = formatter.parse('1 2 3 A B C 4 ACC sh 5 6 7 8 9 8 7 2 # comment')
        cor = {'cycle': 1, 'generated': ['sh', 9.0, 8, 7],
               'interface': ('A', 'B', 'C'), 'length': 4, 'weight': 2.,
               'ordermax': (6.0, 8), 'ordermin': (5.0, 7), 'status': 'ACC'}

        assert_equal_path_dict(out, cor)

        # Also test a shorter (old version) line:
        out = formatter.parse('1 2 3 A B C 4 ACC sh 5 6 7 8 9 8 7')
        cor = {'cycle': 1, 'generated': ['sh', 9.0, 8, 7],
               'interface': ('A', 'B', 'C'), 'length': 4, 'weight': 1.,
               'ordermax': (6.0, 8), 'ordermin': (5.0, 7), 'status': 'ACC'}
        assert_equal_path_dict(out, cor)

        # Let it fail just because:
        out = formatter.parse('1 2 3 A B C 4 ACC sh 5 6 7 8 9 8 7')
        cor = {'cycle': 1, 'generated': ['sh', 9.0, 8, 7],
               'interface': ('A', 'B', 'C'), 'length': 4, 'weight': 2.,
               'ordermax': (6.0, 8), 'ordermin': (5.0, 7), 'status': 'ACC'}

        with pytest.raises(AssertionError) as err:
            assert_equal_path_dict(out, cor)
        assert str(err.value) == 'Different weight'


class TestMethods:
    """Test some of the methods defined in the formatters module."""

    def test_apply_format(self, caplog):
        """Test that we can apply a format."""
        txt = apply_format(12345.7, '{:7.2f}')
        assert txt == '1.2e+04'
        txt = apply_format(-1234568.9, '{:7.2f}')
        assert txt == ' -1e+06'
        txt = apply_format(-1234568.9, '{:8.2f}')
        assert txt == '-1.2e+06'
        txt = apply_format(-1234568.9, '{:9.2f}')
        assert txt == '-1.23e+06'
        txt = apply_format(123.45, '{:>10.2f}')
        assert txt == '    123.45'

    def test_format_number(self, caplog):
        """Test that we can format numbers as expected."""
        txt = format_number(99.9, 0.0, 100, fmtf='{0:<4.2f}',
                            fmte='{0:<4.2e}')
        assert txt == '99.90'
        txt = format_number(100.1, 0.0, 100, fmtf='{0:<4.2f}',
                            fmte='{0:<4.2e}')
        assert txt == '1.00e+02'
        txt = format_number(-100.1, 0.0, 100, fmtf='{0:<4.2f}',
                            fmte='{0:<4.2e}')
        assert txt == '-1.00e+02'

    def test_get_formatter(self, caplog):
        """Test that we can select the log formatter."""
        formatter = get_log_formatter(0)
        assert isinstance(formatter, PyretisLogFormatter)
        # pylint: disable=protected-access
        assert formatter._fmt == LOG_DEBUG_FMT

        formatter = get_log_formatter(100)
        assert isinstance(formatter, PyretisLogFormatter)
        # pylint: disable=protected-access
        assert formatter._fmt == LOG_FMT
