# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test common methods that are used in testing."""
import os
import pytest
from pyretis.testing.compare import (
    compare_files,
    compare_files_lines,
    compare_files_columns,
    compare_files_numerical,
    compare_traj_archive,
    compare_pathensemble_files,
)


HERE = os.path.abspath(os.path.dirname(__file__))
TEST_PATH = os.path.join(HERE, 'files_for_comparison')
# Output messages to expect from the comparisons:
MESSAGES = (
    'Files are equal',
    'Different items in block data',
    'Block comment differs',
    'Block terms differ',
    'Numerical data differ',
    'The number of lines in the files differ',
    'Line 8 differs: # Cycle: 100, status: ACC != # NOPE',
    'Files are not equal',
    'Files differ on line 2, column 6',
    'Files differ on line 0, column 9',
)
# Some example files for testing numerical comparison:
ENERGY1 = os.path.join(TEST_PATH, 'energy1.txt')
ENERGY2 = os.path.join(TEST_PATH, 'energy2.txt')
ENERGY3 = os.path.join(TEST_PATH, 'energy3.txt')
ENERGY4 = os.path.join(TEST_PATH, 'energy4.txt')
ENERGY5 = os.path.join(TEST_PATH, 'energy5.txt')
ORDER1 = os.path.join(TEST_PATH, 'order1.txt')
ORDER2 = os.path.join(TEST_PATH, 'order2.txt')
PATHENSEMBLE1 = os.path.join(TEST_PATH, 'pathensemble1.txt')
PATHENSEMBLE2 = os.path.join(TEST_PATH, 'pathensemble2.txt')
PATHENSEMBLE3 = os.path.join(TEST_PATH, 'pathensemble3.txt')
PATHENSEMBLE4 = os.path.join(TEST_PATH, 'pathensemble4.txt')
PATHENSEMBLE5 = os.path.join(TEST_PATH, 'pathensemble5.txt')
# Some directories for testing archive comparison:
ARCHIVE1 = os.path.join(TEST_PATH, 'traj1')
ARCHIVE2 = os.path.join(TEST_PATH, 'traj1')
ARCHIVE3 = os.path.join(TEST_PATH, 'traj3')
ARCHIVE4 = os.path.join(TEST_PATH, 'traj4')
ARCHIVE5 = os.path.join(TEST_PATH, 'traj5')


class TestCompareFiles:
    """Test that we can compare files."""

    def test_compare_numerical_columns(self):
        """Compare energy.txt and order.txt type of files."""
        # Cases:
        cases = [
            # First, two equal files:
            {
                'args': (ENERGY1, ENERGY2, 'energy'),
                'kwargs': {'skip': None},
                'result': (True, MESSAGES[0]),
            },
            # Now, two files with different column data:
            {
                'args': (ENERGY1, ENERGY3, 'energy'),
                'kwargs': {'skip': None},
                'result': (False, MESSAGES[1]),
            },
            # Now, two files with different block comments:
            {
                'args': (ENERGY1, ENERGY4, 'energy'),
                'kwargs': {'skip': None},
                'result': (False, MESSAGES[2]),
            },
            # Now, two files with different data in a single column,
            # and we request skipping the comparison for this column:
            {
                'args': (ENERGY1, ENERGY5, 'energy'),
                'kwargs': {'skip': ['ekin']},
                'result': (True, MESSAGES[0]),
            },
            # Repeat for these two files, but do not skip:
            {
                'args': (ENERGY1, ENERGY5, 'energy'),
                'kwargs': {'skip': None},
                'result': (False, MESSAGES[3]),
            },
        ]
        for case in cases:
            result = compare_files_columns(
                *case['args'], **case['kwargs']
            )
            assert result == case['result']

    def test_compare_numerical(self):
        """Compare energy.txt and order.txt type of files."""
        # Cases:
        cases = [
            # First, two equal files:
            {
                'args': (ENERGY1, ENERGY2),
                'result': (True, MESSAGES[0]),
            },
            # Two different files:
            {
                'args': (ENERGY1, ENERGY5),
                'result': (False, MESSAGES[4]),
            },
            # Order parameter file:
            {
                'args': (ORDER1, ORDER2),
                'result': (True, MESSAGES[0]),
            },
        ]
        for case in cases:
            result = compare_files_numerical(*case['args'])
            assert result == case['result']

    def test_compare_archive(self):
        """Test the method for comparing archive output."""
        # Cases:
        cases = [
            # Two archives that are equal:
            {
                'args': (ARCHIVE1, ARCHIVE2),
                'result': [],
            },
            # Two archives that differs for a traj.txt file:
            {
                'args': (ARCHIVE1, ARCHIVE3),
                'result': [('traj.txt', 'traj.txt')],
            },
            # Two archives that differs in the number of files:
            {
                'args': (ARCHIVE1, ARCHIVE4),
                'result': [('traj1', 'traj4')],
            },
            # Two archives that have the same number of files but the
            # expected file names does not appear in one of them:
            {
                'args': (ARCHIVE1, ARCHIVE5),
                'result': [
                    ('traj.txt', 'file3.txt'),
                    ('energy.txt', 'file1.txt'),
                    ('order.txt', 'file2.txt'),
                ],
            },
        ]
        for case in cases:
            result = compare_traj_archive(*case['args'])
            result_base = []
            for (i, j) in result:
                result_base.append(
                    (os.path.basename(i), os.path.basename(j))
                )
            assert len(result_base) == len(case['result'])
            for i in case['result']:
                assert i in result_base

    def test_compare_lines(self):
        """Test the method for comparing files line-by-line."""
        # Cases:
        cases = [
            # Two equal files:
            {
                'args': (ENERGY1, ENERGY2),
                'kwargs': {'skip': None},
                'result': (True, MESSAGES[0]),
            },
            # Two unequal files:
            {
                'args': (ENERGY1, ENERGY4),
                'kwargs': {'skip': None},
                'result': (False, MESSAGES[6]),
            },
            # Two unequal files, but we skip the unequal line:
            {
                'args': (ENERGY1, ENERGY4),
                'kwargs': {'skip': [8]},
                'result': (True, MESSAGES[0]),
            },
            # Files with different data, and number of lines:
            {
                'args': (ENERGY1, ORDER1),
                'kwargs': {'skip': [8]},
                'result': (False, MESSAGES[5]),
            },
        ]
        for case in cases:
            result = compare_files_lines(*case['args'], **case['kwargs'])
            assert result == case['result']

    def test_compare_files(self):
        """Test the method for comparing files."""
        cases = [
            # Equal files, numerical comparison:
            {
                'args': (ENERGY1, ENERGY2),
                'kwargs': {'skip': None, 'mode': 'numerical'},
                'result': (True, MESSAGES[0])
            },
            # Equal files, line-by-line:
            {
                'args': (ENERGY1, ENERGY2),
                'kwargs': {'skip': None, 'mode': 'line'},
                'result': (True, MESSAGES[0])
            },
            # Equal files, using filecmp.cmp:
            {
                'args': (ENERGY1, ENERGY2),
                'kwargs': {'skip': None, 'mode': 'cmp'},
                'result': (True, MESSAGES[0])
            },
            # Uneaqual files, numerical.
            # Note: By a numerical comparison, these two files are
            # equal. They only differ in comments and not in the
            # actual numerical data in the file.
            {
                'args': (ENERGY1, ENERGY4),
                'kwargs': {'skip': None, 'mode': 'numerical'},
                'result': (True, MESSAGES[0])
            },
            # Unequal files, line-by-line. This includes comparing
            # comments:
            {
                'args': (ENERGY1, ENERGY4),
                'kwargs': {'skip': None, 'mode': 'line'},
                'result': (False, MESSAGES[6])
            },
            # Unequal files, using filecmp.cmp:
            {
                'args': (ENERGY1, ENERGY4),
                'kwargs': {'skip': None, 'mode': 'cmp'},
                'result': (False, MESSAGES[7])
            },
        ]
        for case in cases:
            result = compare_files(*case['args'], **case['kwargs'])
            assert result == case['result']

    def test_compare_pathensemble_files(self):
        """Test comparison for path ensemble files."""
        cases = [
            # Equal numerical data. Differs in comments:
            {
                'args': (PATHENSEMBLE1, PATHENSEMBLE2),
                'result': (True, MESSAGES[0]),
                'kwargs': {},
            },
            # Different number of lines:
            {
                'args': (PATHENSEMBLE1, PATHENSEMBLE3),
                'result': (False, MESSAGES[5]),
                'kwargs': {},
            },
            # Small numerical difference between files:
            {
                'args': (PATHENSEMBLE1, PATHENSEMBLE4),
                'result': (True, MESSAGES[0]),
                'kwargs': {},
            },
            # Small numerical difference between files,
            # test usage of relative tolerance:
            {
                'args': (PATHENSEMBLE1, PATHENSEMBLE4),
                'result': (False, MESSAGES[9]),
                'kwargs': {'rel_tol': 1e-20},
            },
            # Large numerical difference between files:
            {
                'args': (PATHENSEMBLE1, PATHENSEMBLE5),
                'result': (False, MESSAGES[8]),
                'kwargs': {},
            },
            # Large numerical difference between files,
            # skipping the column with differences:
            {
                'args': (PATHENSEMBLE1, PATHENSEMBLE5),
                'result': (True, MESSAGES[0]),
                'kwargs': {'skip': [6]}
            },
        ]
        for case in cases:
            result = compare_pathensemble_files(
                *case['args'], **case['kwargs']
            )
            assert result == case['result']
        # Check that the cases fail for a line-by-line comparison:
        for case in cases:
            equal, _ = compare_files(*case['args'], skip=None, mode='line')
            assert not equal
