# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the common methods in pyretis.pyvisa.common."""
import colorama
import os
import logging
import subprocess
import pytest
import numpy as np
import tempfile
import shutil
from unittest import mock
from io import StringIO
from pyretis.pyvisa.common import (
    try_data_shift,
    where_from_to,
    read_traj_txt_file,
    find_rst_file,
    _get_trjs,
    find_data,
    get_cv_names,
    recalculate_all,
)
logging.disable(logging.CRITICAL)

HERE = os.path.abspath(os.path.dirname(__file__))
color = colorama.Fore.CYAN


def del_list_ind(data_list, to_delete):
    """Function deleting indices in list."""
    for i in reversed(to_delete):
        del data_list[i]


class TestMethods:
    """Test some of the methods from pyretis.pyvisa.common."""

    def test_try_data_shift(self):
        """Test to try_data_shift and shift_data functions."""
        n = 100  # number of data points
        # The Correct data, normal distr around origo
        # with a slope m=-1/3
        x, y = [], []
        np.random.seed(n)
        for i in range(n):
            r1 = np.random.normal()
            r2 = 3*np.random.normal()
            x.append(r1-r2)
            y.append(r1+r2)
        # The shifted data, to be corrected
        avx, avy = np.average(x), np.average(y)
        sx, sy = [], []
        for i in range(n):
            if x[i] > avx:
                sx.append(-90+x[i])
            else:
                sx.append(x[i])
            if y[i] > avy:
                sy.append(-90+y[i])
            else:
                sy.append(y[i])
        # Evaluate correct data, should do nothing:
        nx, ny = try_data_shift(x, y, False)
        assert x == nx
        assert y == ny
        # Evaluate x-shifted data:
        nx, ny = try_data_shift(sx, y, False)
        assert ny == y
        assert not np.average(sx == np.average(nx))
        # Evaluate y-shifted data:
        nx, ny = try_data_shift(x, sy, False)
        assert x == nx
        assert not np.average(sy == np.average(ny))
        # Evaluate x&y-shifted data:
        nx, ny = try_data_shift(sx, sy, False)
        assert not np.average(sx == np.average(nx))
        assert not np.average(sy == np.average(ny))
        # Evaluate x&y-shifted data, but with op1 flag to hold:
        nx, ny = try_data_shift(sx, sy, True)
        assert sx == nx
        assert not np.average(sy == np.average(ny))

    def test_where_to_from(self):
        """Test to get the right L and R and * trj start and end labels."""
        trj = [1, 2, 3, 4, 5, 6, 4]
        int_as, int_bs = [0, 2, 2, 0, 5, 3], [7, 7, 3, 0, 5, 2]
        results = [('*', '*'), ('L', '*'), ('L', 'R'),
                   ('R', 'R'), ('L', 'L'), ('L', 'R')]

        # For a certain trj, we test different interfaces.
        for int_a, int_b, result in zip(int_as, int_bs,  results):
            start, end = where_from_to(trj, int_a, int_b)
            assert (start, end) == result

        # Let's check the zero ensemble.
        trjs = [[0, 2], [6, 5], [1, 6]]
        results = [('*', '*'), ('R', 'R'), ('*', 'R')]

        for trj, result in zip(trjs, results):
            start, end = where_from_to(trj, int_a=3)
            assert (start, end) == result

    def test_find_rst_file(self):
        """Test for finding the rst file."""
        start_path = os.path.join(HERE, 'test_simulation_dir/000')
        os.chdir(start_path)
        rst_path = os.path.join(HERE, 'test_simulation_dir/input.rst')
        final_path = find_rst_file(start_path)
        assert rst_path == final_path
        final_path = find_rst_file('/')
        assert '/' == '/'

    def test_read_traj_txt_file(self):
        """Test for reading a traj.txt file."""
        txt_file = os.path.join(HERE, 'test_simulation_dir/traj.txt')
        files = read_traj_txt_file(txt_file)
        assert len(files.keys()) == 2

    def test_get_trjs(self):
        """Test for _get_trjs."""
        dir_file = os.path.join(HERE, 'test_simulation_dir')
        trj_list = _get_trjs(dir_file)
        found = sum([1 for trj in trj_list if 'traj.' in trj])
        assert found > 0

    def test_find_data(self):
        """Test for find_data."""
        ass = find_data(HERE, data=HERE + '/test_simulation_dir/traj.txt')
        assert 'traj.' in ass['000']['traj']['0']['traj'][0]

        source = os.path.join(HERE, 'test_simulation_dir/003/order.txt')
        tmp = os.path.join(HERE, 'test_simulation_dir/003/banana')
        os.rename(source, tmp)
        ass = find_data(os.path.join(HERE, 'test_simulation_dir'))
        target = os.path.join(HERE, 'test_simulation_dir',
                              '003', 'traj', 'traj-acc',
                              '0', 'traj', 'traj.xyz')
        assert ass['003']['traj']['0']['traj'][0] == target
        os.rename(tmp, source)

        ass = find_data(HERE, data=HERE + '/test_simulation_dir/traj.txt')
        assert ('test_simulation_dir/traj.txt' in
                ass['000']['traj']['0']['traj'][0])

    def test_recalculate_all(self):
        """Test for recalculate_all."""
        with tempfile.TemporaryDirectory() as tempdir:
            tmp_dir = os.path.join(tempdir, 'test_simulation_dir')
            shutil.copytree(HERE + '/test_simulation_dir', tmp_dir)
            with mock.patch('sys.stdout', new=StringIO()):
                assert not recalculate_all(tmp_dir, 'input.rst')
                file_new = os.path.join(tmp_dir, 'input.rst')
                subprocess.run(['sed', '-i', "s/class = Position/ /",
                                file_new], check=True)
                assert recalculate_all(tmp_dir, 'input.rst')
                assert not recalculate_all(tmp_dir, 'input.rst',
                                           ensemble_names=[])


class TestGetCvNames:
    """Test the get_cv_names label-generation helper."""

    @staticmethod
    def _settings(op_name=None, cv_specs=()):
        """Build a minimal settings dict for the given names."""
        settings = {'orderparameter': {}}
        if op_name is not None:
            settings['orderparameter']['name'] = op_name
        if cv_specs:
            settings['collective-variable'] = []
            for cv_name in cv_specs:
                block = {}
                if cv_name is not None:
                    block['name'] = cv_name
                settings['collective-variable'].append(block)
        return settings

    def test_single_op_default(self):
        """Single OP with no name and no column count uses indexed default."""
        labels = get_cv_names(self._settings())
        assert labels == ['op_1']

    def test_single_op_default_multi_columns(self):
        """A single block expands to fill the column count."""
        labels = get_cv_names(self._settings(), num_columns=4)
        assert labels == ['op_1', 'op_2', 'op_3', 'op_4']

    def test_single_op_string_name_expanded(self):
        """A string name on a single multi-element OP gets indexed."""
        labels = get_cv_names(self._settings(op_name='pippo'),
                              num_columns=3)
        assert labels == ['pippo_1', 'pippo_2', 'pippo_3']

    def test_single_op_string_name_single_column(self):
        """A string name on a single single-element OP keeps the bare name."""
        labels = get_cv_names(self._settings(op_name='pippo'),
                              num_columns=1)
        assert labels == ['pippo']

    def test_single_op_list_name(self):
        """A list name on a single OP must match the column count."""
        labels = get_cv_names(self._settings(op_name=['a', 'b', 'c']),
                              num_columns=3)
        assert labels == ['a', 'b', 'c']

    def test_single_op_list_name_mismatch(self):
        """A list name with the wrong length raises a ValueError."""
        with pytest.raises(ValueError):
            get_cv_names(self._settings(op_name=['a', 'b']),
                         num_columns=3)

    def test_op_and_cv_default(self):
        """Defaults distinguish the main OP from the collective variables."""
        labels = get_cv_names(
            self._settings(cv_specs=[None, None]),
            num_columns=3,
        )
        assert labels == ['op_1', 'cv1_1', 'cv2_1']

    def test_op_and_cv_named(self):
        """Named OP + CVs combine cleanly when sizes match."""
        labels = get_cv_names(
            self._settings(op_name='main', cv_specs=['extra']),
            num_columns=2,
        )
        assert labels == ['main', 'extra']

    def test_op_and_cv_with_list_name(self):
        """A list ``name`` lets a multi-valued CV take several columns."""
        labels = get_cv_names(
            self._settings(op_name='main',
                           cv_specs=[['x', 'y', 'z']]),
            num_columns=4,
        )
        assert labels == ['main', 'x', 'y', 'z']

    def test_block_count_mismatch_raises(self):
        """A column count that cannot be reconciled with the blocks errors."""
        with pytest.raises(ValueError):
            get_cv_names(
                self._settings(op_name='main', cv_specs=['extra']),
                num_columns=5,
            )

    def test_no_blocks(self):
        """Empty settings yield an empty list."""
        assert get_cv_names({}) == []
