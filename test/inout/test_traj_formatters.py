# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the trajectory writers."""
import logging
import pytest
import tempfile
import os
import numpy as np
from numpy.random import rand, random
from pyretis.core import create_box, System, Particles, Path
from pyretis.tools.lattice import generate_lattice
from pyretis.inout.formats.snapshot import (
    SnapshotFile,
    SnapshotFormatter,
)
from pyretis.inout.formats.path import (
    PathIntFormatter,
    PathExtFormatter,
    PathExtFile,
    PathIntFile,
)
from .help import create_external_path, CORRECT_PATH_EXT, set_up_system
logging.disable(logging.CRITICAL)


HERE = os.path.abspath(os.path.dirname(__file__))


def create_test_system():
    """Create a system we can use for testing."""
    xyz, size = generate_lattice('fcc', [3, 3, 3], density=0.9)
    low, high = [], []
    for i in size:
        low.append(i[0])
        high.append(i[1])
    box = create_box(low=low, high=high)
    system = System(units='lj', box=box)
    system.particles = Particles(dim=3)
    for xyzi in xyz:
        system.add_particle(name='Ar', pos=xyzi, vel=np.zeros_like(xyzi))
    system.particles.vel[-1] = np.array([-1.0, 0.123, 1.0])
    return system


def create_path():
    """Setup a simple path for a test."""
    box = create_box(cell=rand(3))
    path = Path(None)
    phasepoints = []
    for _ in range(10):
        pos = rand(10, 3)
        vel = rand(*pos.shape)
        vpot = random()
        ekin = random()
        phasepoint = set_up_system(
            [pos[0][0], pos[1][0]], pos, vel, vpot=vpot, ekin=ekin
        )
        phasepoint.box = box
        phasepoints.append(phasepoint)
        path.append(phasepoint)
    return phasepoints, path


class TestTraj:
    """Test trajectory writing work as intended."""

    def test_txt_writer(self):
        """Test the SnapshotFormatter."""
        system = create_test_system()
        txt_writer = SnapshotFormatter(write_vel=True)
        snapshot = txt_writer.format_snapshot(0, system)
        correct = os.path.join(HERE, 'generated.txt')
        with open(correct, 'r', encoding='utf-8') as fileh:
            for lines1, lines2 in zip(fileh, snapshot):
                assert lines1.rstrip() == lines2.rstrip()

    def test_traj_writer_novel(self):
        """Test the SnapShotFormatter class when we exclude velocities."""
        writer = SnapshotFormatter(write_vel=False, fmt='full')
        phasepoints, path = create_path()
        with tempfile.NamedTemporaryFile() as tmp:
            for step, snapshot in enumerate(path.phasepoints):
                for line in writer.format_snapshot(step, snapshot):
                    string = f'{line}\n'
                    tmp.write(string.encode('utf-8'))
            tmp.flush()
            del writer
            reader = SnapshotFile(tmp.name, 'r')
            for block, snapshot in zip(reader.load(), phasepoints):
                assert np.allclose(block['box'], snapshot.box.length)
                xyz = np.transpose(np.vstack((block['x'],
                                              block['y'],
                                              block['z'])))
                assert np.allclose(xyz, snapshot.particles.get_pos())

    def test_traj_writer_vel(self):
        """Test the SnapshotFormatter class when we include velocities."""
        writer = SnapshotFormatter(write_vel=True, fmt='full')
        phasepoints, path = create_path()
        with tempfile.NamedTemporaryFile() as tmp:
            for step, snapshot in enumerate(path.phasepoints):
                for line in writer.format_snapshot(step, snapshot):
                    string = f'{line}\n'
                    tmp.write(string.encode('utf-8'))
            tmp.flush()
            del writer
            reader = SnapshotFile(tmp.name, 'r')
            for block, snapshot in zip(reader.load(), phasepoints):
                assert np.allclose(block['box'], snapshot.box.length)
                xyz = np.transpose(np.vstack((block['x'],
                                              block['y'],
                                              block['z'])))
                vel = np.transpose(np.vstack((block['vx'],
                                              block['vy'],
                                              block['vz'])))
                assert np.allclose(xyz, snapshot.particles.get_pos())
                assert np.allclose(vel, snapshot.particles.get_vel())

    def test_path_int_writer(self):
        """Test the path internal writer."""
        phasepoints, path = create_path()
        writer = PathIntFormatter()
        fmt = ' '.join([writer._FMT] * 6)
        idxs = 0
        idx = 0
        for i, lines in enumerate(writer.format(0, (path, 'ACC'))):
            if i == 0:
                assert '# Cycle: 0, status: ACC', lines
            else:
                if lines.startswith('Snapshot'):
                    assert f'Snapshot: {idxs}' == lines
                    idxs += 1
                    idx = 0
                else:
                    posvel = fmt.format(
                        *phasepoints[idxs - 1].particles.get_pos()[idx],
                        *phasepoints[idxs - 1].particles.get_vel()[idx]
                    )
                    assert lines == posvel
                    idx += 1

    def test_pathint_read_write(self):
        """Test that we can read write with PathIntFormatter."""
        writer = PathIntFormatter()
        statuses = ('ACC', 'BWI', 'FTL')
        path_phasepoints = []
        with tempfile.NamedTemporaryFile() as tmp:
            for i, status in zip(range(3), statuses):
                phasepoints, path = create_path()
                path_phasepoints.append(phasepoints)
                for line in writer.format(i, (path, status)):
                    tmp.write(f'{line}\n'.encode('utf-8'))
            tmp.flush()
            del writer
            reader = PathIntFile(tmp.name, 'r')
            for status, path, phasepoints in zip(statuses,
                                                 reader.load(),
                                                 path_phasepoints):
                assert status == path['comment'][0].split()[-1]
                for snapshot, phasepoint in zip(path['data'], phasepoints):
                    particles = phasepoint.particles
                    assert np.allclose(snapshot['pos'], particles.get_pos())
                    assert np.allclose(snapshot['vel'], particles.get_vel())

    def test_pathint_read_error(self):
        """Test what happens in we read faulty input with PathIntFile."""
        filename = os.path.join(HERE, 'traj-error1.txt')
        reader = PathIntFile(filename, 'r').load()
        next(reader)
        with pytest.raises(ValueError):
            next(reader)
        filename = os.path.join(HERE, 'traj-error2.txt')
        reader = PathIntFile(filename, 'r').load()
        with pytest.raises(ValueError):
            next(reader)

    def test_path_ext_writer(self):
        """Test the path external writer."""
        path, _ = create_external_path()
        writer = PathExtFormatter()
        for corr, snap in zip(CORRECT_PATH_EXT,
                              writer.format(0, (path, 'ACC'))):
            assert corr == snap

    def test_path_ext_read_write(self):
        """Test the read/write for the PathExtFormatter."""
        writer = PathExtFormatter()
        statuses = ('ACC', 'BWI', 'FTL')
        all_path_data = []
        with tempfile.NamedTemporaryFile() as tmp:
            for i, status in zip(range(3), statuses):
                path, data = create_external_path(random_length=True)
                all_path_data.append(data)
                for line in writer.format(i, (path, status)):
                    string = f'{line}\n'
                    tmp.write(string.encode('utf-8'))
            tmp.flush()
            del writer
            reader = PathExtFile(tmp.name, 'r')
            for status, path, data in zip(statuses, reader.load(),
                                          all_path_data):
                assert status == path['comment'][0].split()[-1]
                for snapshot, datai in zip(path['data'], data):
                    assert int(snapshot[0]) == datai[0]
                    assert snapshot[1] == datai[1]
                    if datai[2] is None:
                        assert int(snapshot[2]) == 0
                    else:
                        assert int(snapshot[2]) == datai[2]
                    vel = snapshot[3] == '-1'
                    assert vel == datai[3]

    def test_path_eq_none(self):
        """
        If nullmoves=False, the first and last ensembles are sometimes not
        updated (in first and last path ensembles when they are excluded from
        a replica exchange move) and nothing is supposed to be written to file.
        This is tested here.
        """
        writer = PathExtFormatter()
        assert list(writer.format(4, (None, 'ACC'))) == []
        writer = PathIntFormatter()
        assert list(writer.format(13, (None, 'ACC'))) == []
