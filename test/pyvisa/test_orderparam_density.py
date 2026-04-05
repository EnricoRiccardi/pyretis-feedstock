# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the oder_parameter methods in pyretis.pyvisa.orderparam_density."""
import logging
import os
import warnings
import pytest
from io import StringIO
import numpy as np
import pandas as pd
from subprocess import DEVNULL, STDOUT, check_call
from pyretis.pyvisa.orderparam_density import (PathDensity,
                                               PathVisualize,
                                               remove_nan,
                                               pyvisa_compress)
from unittest import mock

pytest.importorskip('tables', exc_type=ImportError)
logging.disable(logging.CRITICAL)
warnings.filterwarnings('ignore', category=pd.io.pytables.PerformanceWarning)

BASE = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                    'test_simulation_dir')
BASE_ZIP = os.path.join(BASE, 'pyvisa_compressed_data.hdf5.zip')
INPUTFILE = os.path.join(BASE, 'input.rst')

CORRECT = {
    'path': ['000', '001', '002', '003', '004', '005'],
    'op_labels': ['op1', 'op2'],
    'energy_labels': ['time', 'cycle', 'potE', 'kinE', 'totE'],
    'interfaces': [-0.9, -0.8, -0.7, -0.5, -0.4, -0.3, 1.0],
    'num_op': 2
}

CORRECT_LENGTHS_SLOW = {
    '000 A': [968, 968, 968, 968, 3872, 3872, 3872, 3872],
    '000 R': [0, 0, 0, 0, 1500, 1500, 1500, 1500],
    '001 A': [247, 247, 247, 247, 1235, 1235, 1235, 1235],
    '001 R': [454, 454, 454, 454, 2719, 2719, 2719, 2719],
    '002 A': [247, 247, 247, 247, 2223, 2223, 2223, 2223],
    '002 R': [454, 454, 454, 454, 2719, 2719, 2719, 2719],
    '003 A': [247, 247, 247, 247, 2223, 2223, 2223, 2223],
    '003 R': [454, 454, 454, 454, 2719, 2719, 2719, 2719],
}

LENGTHS = {'000': {0: [484], 1: [500], 2: [500], 3: [500],
                   4: [484], 5: [484],
                   6: [484], 7: [484], 8: [484], 9: [484], 10: [484]},
           '001': {0: [247], 1: [500], 2: [492], 3: [247],
                   4: [454], 5: [247], 6: [247],
                   7: [247], 8: [454], 9: [365], 10: [454]},
           '002': {0: [247], 1: [500], 2: [492, 500, 492],
                   3: [247, 247], 4: [454, 454],
                   5: [247], 6: [247],
                   7: [247], 8: [454], 9: [365], 10: [454]},
           '005': {0: [17], 1: [19], 2: [19]}}

STATUS = {'005': {0: 'ACC', 1: 'ACC', 2: 'NCR',
                  3: 'NCR', 4: 'ACC', 5: 'ACC', 6: 'ACC'}}
MOVES = {'005': {0: 'ki', 1: 'tr', 2: 's+',
                 3: 's+', 4: 'tr', 5: 's-', 6: 'tr'}}


def init_pathdensity(base=BASE, ifile=INPUTFILE):
    """Initialize PathDensity with input file."""
    with mock.patch('sys.stdout', new=StringIO()):
        data_dict = PathDensity(base, ifile)
    return data_dict


def init_pathvisualize(base=BASE, pfile=INPUTFILE):
    """Initialize PathVisualize."""
    with mock.patch('sys.stdout', new=StringIO()):
        visual = PathVisualize(base, pfile=pfile)
    return visual


def compare_keys(dict1, correct):
    """Test if given dictionary contains the correct set of keys."""
    check = True
    for key in correct:
        check = True if key in dict1 else False
    return check


class TestMethods:
    """Testing class of pyretis.pyvisa.orderparam_density."""

    def test_PathVisualize_blank(self):
        """Testing initiate PathVisualize."""
        mori = PathVisualize()
        mori.pfile = os.path.join(BASE, 'traj.txt')
        with pytest.raises(ValueError) as err:
            mori.load_whatever()
            assert 'not recognised' in err.value

    def test_PathDensity_blank(self):
        """Testing initiate PathDensity with a compressed input file."""
        with pytest.raises(FileNotFoundError) as err:
            init_pathdensity(ifile='blank.hdf5')
            assert 'File blank.hdf5 does not exist.' in err.value

    def test_PathDensity(self):
        """Testing the PathDensity class creation from file."""
        ipd = init_pathdensity(BASE, INPUTFILE)
        assert ipd.iofile == INPUTFILE
        assert compare_keys(ipd.infos, CORRECT.keys())
        assert compare_keys(ipd.infos['ensemble_names'],
                            CORRECT['path'])
        assert compare_keys(ipd.infos['op_labels'],
                            CORRECT['op_labels'])
        assert compare_keys(ipd.infos['energy_labels'],
                            CORRECT['energy_labels'])
        assert compare_keys(ipd.infos['interfaces'],
                            CORRECT['interfaces'])
        assert ipd.infos['num_op'] == CORRECT['num_op']

        re_d = PathDensity(basepath='.', iofile=ipd)
        assert isinstance(re_d, PathDensity)
        # test minor functionalities
        re_d.get_traj_energy('no_energy')

    def test_PathDensity_instance(self):
        """Test comparing PathDensity object."""
        pe_d = init_pathdensity()
        pe_e = init_pathdensity()
        assert pe_e.infos == pe_d.infos
        assert pe_e.traj_dict == pe_d.traj_dict

    def test_walk_Dirs(self):
        """Testing the walk_Dirs function of PathDensity."""
        pe_d = init_pathdensity()
        with mock.patch('sys.stdout', new=StringIO()):
            pe_d.walk_dirs()

    def test_pyvisa_compress(self):
        """Testing the compression of input files."""
        pyvisa_dict_hdf5 = {'only_order': True}
        with mock.patch('sys.stdout', new=StringIO()):
            pyvisa_compress(BASE, INPUTFILE,
                            pyvisa_dict_hdf5)
        assert os.path.isfile('pyvisa_compressed_data.hdf5.zip')

        with pytest.raises(ImportError) as err:
            pyvisa_compress('.', 'pyvisa_compressed_data.hdf5.zip', {})
            assert 'Cannot compress an already' in err.value
        check_call(['unzip', 'pyvisa_compressed_data.hdf5.zip'],
                   stdout=DEVNULL, stderr=STDOUT)

        with mock.patch('sys.stdout', new=StringIO()):
            re_d = PathDensity('.', 'pyvisa_compressed_data.hdf5')
        assert isinstance(re_d, PathDensity)
        os.remove('pyvisa_compressed_data.hdf5.zip')
        os.remove('pyvisa_compressed_data.hdf5')

    def test_hdf5(self):
        """Testing the saving of a hdf5 file."""
        pe_d = init_pathdensity()
        with mock.patch('sys.stdout', new=StringIO()):
            pe_d.walk_dirs()
            # This is to remove the Deprecation Warning due to pandas
            with mock.patch('sys.stderr', new=StringIO()):
                pe_d.hdf5_data()
        assert os.path.isfile('pyvisa_compressed_data.hdf5.zip')
        os.remove('pyvisa_compressed_data.hdf5.zip')

    def test_hdf5ing_and_loading(self):
        """Test for hdf5ing data to file and load with PathVisualize."""
        re_d = init_pathdensity()
        with mock.patch('sys.stdout', new=StringIO()):
            re_d.walk_dirs()
            # This is to remove the Deprecation Warning due to pandas.
            with mock.patch('sys.stderr', new=StringIO()):
                re_d.hdf5_data()
        assert os.path.isfile('pyvisa_compressed_data.hdf5.zip')
        init_pathvisualize('pyvisa_compressed_data.hdf5.zip')
        os.rename('pyvisa_compressed_data.hdf5.zip', 'myhdf.hdf6.zip')
        with mock.patch('sys.stdout', new=StringIO()):
            with pytest.raises(TypeError) as err:
                init_pathvisualize('myhdf.hdf6.zip')
                assert 'recognised' in err
        os.remove('myhdf.hdf6.zip')

    def test_get_trajectory_op_length(self):
        """Test for correct lengths of collected trajectories."""
        re_d = init_pathdensity()
        with mock.patch('sys.stdout', new=StringIO()):
            re_d.walk_dirs()
        for ens in LENGTHS.keys():
            for cyc in LENGTHS[ens]:
                for idx in range(0, len(LENGTHS[ens][cyc])):
                    assert (LENGTHS[ens][cyc][0] ==
                            re_d.traj_dict[ens][cyc].info['length'])

    def test_load_trajectory(self):
        """Test for loading data."""
        return
        re_d = init_pathdensity()
        with mock.patch('sys.stdout', new=StringIO()):
            re_d.walk_dirs()
            with mock.patch('sys.stderr', new=StringIO()):
                re_d.hdf5_data()
        re_w = init_pathvisualize('pyvisa_compressed_data.hdf5.zip')
        os.remove(BASE_ZIP)

        criteria = {'x': 'op1', 'y': 'op2', 'z': 'time',
                    'cycles': re_d.infos['cycles'],
                    'weight': False}

        for ens in re_d.infos['ensemble_names']:
            x, y, z = [], [], []

            for cyc in range(criteria['cycles'][0], criteria['cycles'][1] + 1):
                if cyc not in re_d.traj_dict[ens].keys():
                    continue
                else:
                    traj = re_d.traj_dict[ens][cyc]
                    for xyz in {'x': x, 'y': y, 'z': z}.items():
                        if criteria[xyz[0]] == 'time':
                            xyz[1].extend(range(0, traj.info['length']))
                        elif criteria[xyz[0]] == 'cycle':
                            xyz[1].extend([traj.info['cycle']] *
                                          traj.info['length'])
                        elif criteria[xyz[0]] == 'None':
                            xyz[1].extend([1] * traj.info['length'])
                        elif criteria[xyz[0]] in traj.frames.columns:
                            xyz[1].extend(traj.frames[criteria[xyz[0]]])

            criteria['ensemble_name'] = ens
            x2, y2, z2, _ = re_w.load_traj(criteria)

            assert x == x2
            assert y == y2
            assert z == z2

    def test_move_status(self):
        """Test for move status."""
        pe_d = init_pathdensity()
        with mock.patch('sys.stdout', new=StringIO()):
            pe_d.walk_dirs()
        for cyc in pe_d.traj_dict['005'].keys():
            traj = pe_d.traj_dict['005'][cyc]
            assert traj.info['status'] == STATUS['005'][cyc]
            assert traj.info['MC-move'] == MOVES['005'][cyc]

    def test_load_traj_acc_rej_both_hdf5(self):
        """Load ACC, REJ and BOTH trajectories from hdf5."""
        pe_d = init_pathdensity()
        with mock.patch('sys.stdout', new=StringIO()):
            pe_d.walk_dirs()
            pe_d.hdf5_data()
        pe_w = init_pathvisualize('.', 'pyvisa_compressed_data.hdf5.zip')
        os.remove('pyvisa_compressed_data.hdf5.zip')
        c_acc = {'status': 'ACC', 'x': 'op1', 'y': 'potE',
                 'z': 'time', 'cycles': pe_d.infos['cycles'],
                 'ensemble_name': '005', 'weight': False}
        c_rej = {'status': 'REJ', 'x': 'op2', 'y': 'totE',
                 'z': 'cycle', 'cycles': pe_d.infos['cycles'],
                 'ensemble_name': '005', 'weight': False}
        c_both = {'status': 'BOTH', 'x': 'op1', 'y': 'kinE',
                  'z': 'None', 'cycles': pe_d.infos['cycles'],
                  'ensemble_name': '005', 'weight': False}
        x_a, y_a, z_a = [], [], []
        x_r, y_r, z_r = [], [], []
        x_b, y_b, z_b = [], [], []
        for cycle in pe_d.traj_dict['005'].keys():
            traj = pe_d.traj_dict['005'][cycle]
            if traj.info['status'] == 'ACC':
                x_a.extend(traj.frames['op1'])
                y_a.extend(traj.frames['potE'])
                z_a.extend(range(0, traj.info['length']))
            else:
                x_r.extend(traj.frames['op2'])
                y_r.extend(traj.frames['totE'])
                z_r.extend([traj.info['cycle']] *
                           traj.info['length'])

            x_b.extend(traj.frames['op1'])
            y_b.extend(traj.frames['kinE'])
            z_b.extend([1] * traj.info['length'])

        x2_a, y2_a, z2_a, _ = pe_w.load_traj(c_acc)
        # Use assert_array_equal which correctly handles NaN == NaN
        np.testing.assert_array_equal(x_a, x2_a)
        np.testing.assert_array_equal(y_a, y2_a)
        np.testing.assert_array_equal(z_a, z2_a)

        x2_r, y2_r, z2_r, _ = pe_w.load_traj(c_rej)
        np.testing.assert_array_equal(x_r, x2_r)
        np.testing.assert_array_equal(y_r, y2_r)
        np.testing.assert_array_equal(z_r, z2_r)

        x2_b, y2_b, z2_b, _ = pe_w.load_traj(c_both)
        np.testing.assert_array_equal(x_b, x2_b)
        np.testing.assert_array_equal(y_b, y2_b)
        np.testing.assert_array_equal(z_b, z2_b)

    def test_no_op(self):
        """Remove Trajectory with no OP values."""
        pe_d = init_pathdensity()
        with mock.patch('sys.stdout', new=StringIO()):
            pe_d.walk_dirs()
        assert 3 not in pe_d.traj_dict['005'].keys()
        assert 0 in pe_d.traj_dict['005'].keys()

    def test_more_ops(self):
        """Remove Trajectory with no OP values."""
        pe_d = init_pathdensity()
        with mock.patch('sys.stdout', new=StringIO()):
            pe_d.walk_dirs()
        assert 5 == len(pe_d.traj_dict['005'][4].frames.columns)
        assert 5 == len(pe_d.traj_dict['005'][5].frames.columns)
        assert 2 == len(pe_d.traj_dict['004'][1].frames.columns)

    def test_remove_nan(self):
        """Test for stability in the most dramatic case.

        Note that the remove_nan function alters the input in place.

        """
        frames = {'0': {'op1': float('nan'), 'op2': 0},
                  '1': {'op1': float('nan'), 'op2': float('nan')},
                  '2': {'op1': float('nan'), 'op2': float('nan')},
                  '3': {'op1': float('nan'), 'op2': float('nan')}}
        data = pd.DataFrame.from_dict(frames, orient='index')
        data['op3'] = [float('nan'), float('nan'), float('nan'), 4]
        data['op4'] = [0, 2, float('nan'), 4]

        remove_nan(data)
        # Nothing can be fixed here.
        assert not data['op1'].iloc[0] * 0 == 0
        # Nothing can be fixed here and the start shall remain the same.
        assert data['op2'].iloc[0] == 0
        assert not data['op2'].iloc[3] * 0 == 0
        # Test nan at the beginning, it should be fixed.
        assert data['op3'].iloc[0] == 4
        assert data['op3'].iloc[3] == 4
        # Test nan disappearing (most common case).
        assert data['op4'].iloc[2] == 4

        # Test that the function works with normal lists also.
        bad = [0, float('nan'), 1]
        not_bad = [0, 1, 1]
        remove_nan([bad])
        assert bad == not_bad
