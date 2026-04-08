# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A test module for the xyz in/out module."""
import logging
import pytest
import tempfile
import filecmp
import os
import gzip
import numpy as np
from pyretis.inout.formats.xyz import (
    format_xyz_data,
    read_xyz_file,
    convert_snapshot,
    write_xyz_file,
    write_xyz_trajectory,
    reverse_xyz_file,
    xyz_merge,
    txt_to_xyz,
)


logging.disable(logging.CRITICAL)
HERE = os.path.abspath(os.path.dirname(__file__))


class TestXYZIO:
    """Test the xyz format reading & writing."""

    def test_format_data(self):
        """Test that we format data as expected."""
        correct_xyz = [
            ['2', 'PyRETIS XYZ writer',
             'X         1.000000000     2.000000000     3.000000000',
             'X         4.000000000     5.000000000     6.000000000'],
            ['2', 'Custom header',
             'X         1.000000000     2.000000000     3.000000000',
             'X         4.000000000     5.000000000     6.000000000'],
            ['2', 'PyRETIS XYZ writer',
             'A         1.000000000     2.000000000     3.000000000',
             'B         4.000000000     5.000000000     6.000000000'],
            ['2', 'PyRETIS XYZ writer',
             'X          1.0      2.0      3.0',
             'X          4.0      5.0      6.0'],
            ['2', 'PyRETIS XYZ writer',
             ('X         1.000000000     2.000000000     3.000000000'
              '     9.000000000     8.000000000     7.000000000'),
             ('X         4.000000000     5.000000000     6.000000000'
              '     6.000000000     5.000000000     4.000000000')],
        ]
        pos = np.array([[1., 2., 3.], [4., 5., 6.]])
        vel = np.array([[9., 8., 7.], [6., 5., 4.]])
        cases = [
            {'pos': pos},
            {'pos': pos, 'header': 'Custom header'},
            {'pos': pos, 'names': ['A', 'B']},
            {'pos': pos, 'fmt': '{0:5s} {1:8.1f} {2:8.1f} {3:8.1f}'},
            {'pos': pos, 'vel': vel},
        ]
        for case, correct in zip(cases, correct_xyz):
            formatted = format_xyz_data(
                case.get('pos', None),
                vel=case.get('vel', None),
                names=case.get('names', None),
                header=case.get('header', None),
                fmt=case.get('fmt', None),
            )
            for i, j in zip(formatted, correct):
                assert i == j

    def test_read_file(self):
        """Test that we can read a single snapshot."""
        correct_data = [
            {'box': np.array([1., 2., 3.]),
             'xyz': np.array([[0., 0., 0.], [0.5, 0.5, 0.5], [0.5, 0.5, 0.],
                              [0.5, 0., 0.5], [0., 0.5, 0.5]]),
             'vel': np.zeros((5, 3)),
             'names': ['Ba', 'Hf', 'O', 'O', 'O']},
            {'box': np.array([1., 2., 3., 4., 5., 6.]),
             'xyz': np.array([[1., 2., 3.], [7., 8., 9.], [0.1, 0.2, 0.3]]),
             'vel': np.array([[4., 5., 6.], [10., 11., 12.], [0.6, 0.7, 0.8]]),
             'names': ['X', 'Y', 'Z']},
        ]
        infiles = ('config_box.xyz', 'config_vel.xyz')
        for files, correct in zip(infiles, correct_data):
            infile = os.path.join(HERE, files)
            snapshot = read_xyz_file(infile)
            box, xyz, vel, names = convert_snapshot(next(snapshot))
            data = (box, xyz, vel)
            var = ('box', 'xyz', 'vel')
            for i, j in zip(data, var):
                assert np.allclose(i, correct[j])
            for i, j in zip(names, correct['names']):
                assert i == j

    def test_xyz_file(self):
        """Test that we can write an xyz-file."""
        infile = os.path.join(HERE, 'config_read.xyz')
        snapshot = next(read_xyz_file(infile))
        _, xyz, vel, names = convert_snapshot(snapshot)
        header = snapshot['header']
        with tempfile.NamedTemporaryFile() as temp:
            write_xyz_file(temp.name, xyz, vel=vel, names=names,
                           header=header)
            temp.flush()
            compare = filecmp.cmp(temp.name, infile)
            assert compare

    def test_write_trajectory(self):
        """Test that we can write a xyz-trajectory."""
        posin = [np.ones((2, 3))*1.5, np.ones((2, 3))*2.]
        velin = [np.ones((2, 3)), np.ones((2, 3))*2.]
        namesin = ['A', 'B']
        boxin = [[3., 3., 3.], [5., 5., 5., 5., 5., 5.]]
        step = 0
        headers = [
            '# Step: 0 Box:    3.0000    3.0000    3.0000',
            ('# Box:    5.0000    5.0000    5.0000'
             '    5.0000    5.0000    5.0000'),
        ]
        with tempfile.NamedTemporaryFile() as temp:
            for i, j, k in zip(posin, velin, boxin):
                write_xyz_trajectory(temp.name, i, j, namesin,
                                     k, step=step, append=True)
                step = None
            temp.flush()
            for i, snapshot in enumerate(read_xyz_file(temp.name)):
                box, xyz, vel, names = convert_snapshot(snapshot)
                assert np.allclose(xyz, posin[i])
                assert np.allclose(vel, velin[i])
                assert np.allclose(box, boxin[i])
                assert snapshot['header'] == headers[i]
                for j, k in zip(namesin, names):
                    assert j == k

    def test_reverse_xyz_file(self):
        """Test that we can reverse an xyz-file."""
        infile = os.path.join(HERE, 'backward.xyz')
        with tempfile.NamedTemporaryFile() as temp:
            reverse_xyz_file(infile, temp.name)
            temp.flush()
            correct = reversed([snap for snap in read_xyz_file(infile)])
            for snapshot, corr in zip(read_xyz_file(temp.name), correct):
                assert snapshot['header'] == corr['header']
                for key, val in snapshot.items():
                    if key in ('header', 'atomname'):
                        assert val == corr[key]
                    else:
                        assert np.allclose(val, corr[key])

    def test_reverse_larger_file(self):
        """Test that we can reverse a larger xyz-file."""
        infile = os.path.join(HERE, 'traj.xyz')
        with tempfile.NamedTemporaryFile() as temp:
            reverse_xyz_file(infile, temp.name)
            temp.flush()
            for i, snapshot in enumerate(read_xyz_file(temp.name)):
                step = int(snapshot['header'].split()[2])
                assert i == step - 1
                posx = np.ones(3) * (i + 1)
                assert np.allclose(snapshot['x'], posx)

    def test_merge(self):
        """Test that we can merge forward and backward xyz-files."""
        backward = os.path.join(HERE, 'backward.xyz')
        forward = os.path.join(HERE, 'forward.xyz')
        merged_correct = os.path.join(HERE, 'merged.xyz')
        with tempfile.NamedTemporaryFile() as temp:
            xyz_merge(backward, forward, temp.name)
            temp.flush()
            compare = filecmp.cmp(temp.name, merged_correct)
            assert compare

    def test_txt_convert(self):
        """Test that we can convert from the internal txt format."""
        infile = os.path.join(HERE, 'traj.txt.gz')
        output_name = os.path.join(HERE, 'DELETE')
        out = []
        correct = ['correct-000-ACC.xyz.gz', 'correct-000-ACC.xyz.gz']
        with tempfile.NamedTemporaryFile() as temp:
            with open(temp.name, 'wb') as output:
                with gzip.open(infile, 'rb') as gfile:
                    output.write(gfile.read())
            temp.flush()
            out = txt_to_xyz(temp.name, output_name, None,
                             selection='ACC', nzero=3)
        for files, correctf in zip(out, correct):
            correct_name = os.path.join(HERE, correctf)
            with tempfile.NamedTemporaryFile() as temp:
                with open(temp.name, 'wb') as output:
                    with gzip.open(correct_name, 'rb') as gfile:
                        output.write(gfile.read())
                temp.flush()
                compare = filecmp.cmp(temp.name, files)
                assert compare
            os.remove(files)
