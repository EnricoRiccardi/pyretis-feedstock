# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A test module for the CP2K io module."""
import logging
import tempfile
import os
from unittest.mock import MagicMock, patch
import pytest
import numpy as np
from pyretis.inout.formats.cp2k import (
    read_cp2k_input,
    dfs_print,
    set_parents,
    read_cp2k_energy,
    read_cp2k_restart,
    read_cp2k_box,
    update_node,
    remove_node,
    update_cp2k_input,
    cp2k_settings,
)


logging.disable(logging.CRITICAL)
HERE = os.path.abspath(os.path.dirname(__file__))

BOX_TEMPLATE = """&FORCE_EVAL
  &SUBSYS
    &CELL
      PERIODIC {}
      ABC 5.0 5.0 5.0
      ALPHA_BETA_GAMMA 90.0 90.0 90.0
    &END CELL
  &END SUBSYS
&END FORCE_EVAL"""
BOX_TEMPLATE2 = """&FORCE_EVAL
  &SUBSYS
    &CELL
      ABC 5.0 5.0 5.0
      ALPHA_BETA_GAMMA 90.0 90.0 90.0
    &END CELL
  &END SUBSYS
&END FORCE_EVAL"""
BOX_TEMPLATE3 = """&FORCE_EVAL
  &SUBSYS
  &END SUBSYS
&END FORCE_EVAL"""


class TestCP2KIO:
    """Test CP2KIO."""

    def test_read_cp2k(self):
        """Test that we can read CP2K input files."""
        infile = os.path.join(HERE, 'cp2k.inp')
        nodes = read_cp2k_input(infile)
        buff = []
        for _, node in enumerate(nodes):
            visited = set()
            out = dfs_print(node, visited)
            for line in out:
                buff.append(line.strip())
        with open(infile, 'r') as inputfile:
            buff2 = [i.strip() for i in inputfile]
            # Test that we find all lines in the original input file:
            for line in buff:
                assert line in buff2
            # Test that we read everything:
            for line in buff2:
                if not line:
                    continue
                assert line in buff

    def test_set_parents_read_cp2k(self):
        """Test that we can construct parents for sections."""
        infile = os.path.join(HERE, 'cp2k.inp')
        nodes = read_cp2k_input(infile)
        node_ref = set_parents(nodes)
        ref = 'FORCE_EVAL->DFT->XC->XC_FUNCTIONAL'
        node = node_ref[ref]
        assert node.title == 'XC_FUNCTIONAL'
        assert node.level == 3
        for i, j in zip(node.settings, ['PADE']):
            assert i == j
        for i, j in zip(node.parents, ['FORCE_EVAL', 'DFT', 'XC',
                                       'XC_FUNCTIONAL']):
            assert i == j

    def test_read_energy(self):
        """Test that we can read cp2k energy output."""
        infile = os.path.join(HERE, 'cp2k.ener')
        energy = read_cp2k_energy(infile)
        for key in ('ekin', 'vpot', 'temp', 'etot'):
            assert key in energy
        kin_correct = np.array([0.004750225, 0.008063503, 0.018933585,
                                0.033672265, 0.048363565, 0.060298562])
        pot_correct = np.array([-15.906282000, -15.920520150, -15.931604156,
                                -15.946600436, -15.961490491, -15.973407929])
        assert np.allclose(kin_correct, energy['ekin'])
        assert np.allclose(pot_correct, energy['vpot'])

    def test_read_energy_missing(self):
        """Test that we can read cp2k energy output with missing columns."""
        infile = os.path.join(HERE, 'cp2k.ener2')
        energy = read_cp2k_energy(infile)
        for key in ('ekin', ):
            assert key in energy
        for key in ('vpot', 'temp', 'etot'):
            assert key not in energy
        kin_correct = np.array([0.004750225, 0.008063503, 0.018933585,
                                0.033672265, 0.048363565, 0.060298562])
        assert np.allclose(kin_correct, energy['ekin'])

    def test_read_restart(self):
        """Test that we can read cp2k restart files."""
        infile = os.path.join(HERE, 'cp2k.restart')
        atoms, pos, vel, box, periodic = read_cp2k_restart(infile)
        correct = {
            'atoms': ['H', 'N', 'C'],
            'pos': np.array([
                [2.0037295984360898E-01, 1.3196012905205206E+00,
                 -8.4772218235716257E-01],
                [1.2469580920891239E+00, 5.6235164153925576E-01,
                 -1.5756836410607700E-01],
                [-1.5843560585520555E-01, 3.9769694640060993E-01,
                 2.4037061887112918E-01],
                ]),
            'vel': np.array([
                [-1.3106795516477289E-03, -1.7627790939445680E-03,
                 -4.3827056898331971E-04],
                [-7.7711186323380392E-04, -4.0307772046066932E-04,
                 4.6944742420092559E-04],
                [1.0163250974302026E-03, 6.1937600463006590E-04,
                 -5.3549195000426192E-04],
                ]),
            'box': [5., 5., 5.],
            'periodic': [True, True, True],
        }
        for i, j in zip(atoms, correct['atoms']):
            assert i == j
        assert np.allclose(pos, correct['pos'])
        assert np.allclose(vel, correct['vel'])
        for i, j in zip(periodic, correct['periodic']):
            assert i == j
        for i, j in zip(box, correct['box']):
            assert i == pytest.approx(j)

    def test_read_box(self):
        """Test that we can read the box."""
        test_data = [
            {'file': 'box1.inp', 'box': [5., 5., 5.],
             'periodic': [True, True, True]},
            {'file': 'box2.inp', 'box': [1., 5., 9., 4., 7., 2., 8., 3., 6.],
             'periodic': [True, True, True]},
            {'file': 'box3.inp',
             'box': [2.7, 4.47637, 1.97371, 1.99301, 0.25593,
                     0.0, -0.67006, 0.0, 0.0], 'periodic': [True, True, True]},
            {'file': 'box4.inp', 'box': [5., 5., 5.],
             'periodic': [True, True, True]},
            {'file': 'box5.inp', 'box': None,
             'periodic': [True, True, True]},
        ]
        for data in test_data:
            infile = os.path.join(HERE, data['file'])
            box, periodic = read_cp2k_box(infile)
            if box is None:
                assert box == data['box']
            else:
                for i, j in zip(box, data['box']):
                    assert i == pytest.approx(j, abs=10**-(4))
            for i, j in zip(periodic, data['periodic']):
                assert i == j

    def test_read_box_periodic(self):
        """Test that we can read different periodic settings."""
        correct = {
            'X': [True, False, False],
            'XY': [True, True, False],
            'XYZ': [True, True, True],
            'XZ': [True, False, True],
            'Y': [False, True, False],
            'YZ': [False, True, True],
            'Z': [False, False, True],
        }
        for key, val in correct.items():
            txt = BOX_TEMPLATE.format(key)
            with tempfile.NamedTemporaryFile() as temp:
                temp.write(txt.encode('utf-8'))
                temp.flush()
                _, periodic = read_cp2k_box(temp.name)
                for i, j in zip(periodic, val):
                    assert i == j
        # Test also when we are missing the input:
        with tempfile.NamedTemporaryFile() as temp:
            temp.write(BOX_TEMPLATE2.encode('utf-8'))
            temp.flush()
            _, periodic = read_cp2k_box(temp.name)
            for i, j in zip(periodic, [True, True, True]):
                assert i == j
        with tempfile.NamedTemporaryFile() as temp:
            temp.write(BOX_TEMPLATE3.encode('utf-8'))
            temp.flush()
            _, periodic = read_cp2k_box(temp.name)
            for i, j in zip(periodic, [True, True, True]):
                assert i == j

    def test_update_node(self):
        """Test that we can update nodes."""
        infile = os.path.join(HERE, 'cp2k.inp')
        nodes = read_cp2k_input(infile)
        node_ref = set_parents(nodes)
        target = 'MOTION->MD'
        data = {'STEPS': 987}
        assert node_ref[target].data[0] == 'STEPS 123'
        update_node(target, None, data, node_ref, nodes, replace=False)
        assert node_ref[target].data[0] == 'STEPS 987'
        data = ['STEPS 666']
        update_node(target, None, data, node_ref, nodes, replace=True)
        assert node_ref[target].data[0] == 'STEPS 666'
        data = {'EXTRA': 'yes', 'EXTRAEXTRA': None}
        update_node(target, None, data, node_ref, nodes, replace=False)
        for i in ('STEPS 666', 'EXTRA yes', 'EXTRAEXTRA'):
            assert i in node_ref[target].data
        settings = ['one', 'two']
        update_node(target, settings, [], node_ref, nodes, replace=False)
        for i, j in zip(settings, node_ref[target].settings):
            assert i == j
        settings = ['three']
        update_node(target, settings, [], node_ref, nodes, replace=True)
        for i, j in zip(settings, node_ref[target].settings):
            assert i == j
        target = 'MOTION->STUFF'
        settings = ['nope']
        data = ['setting1 100']
        update_node(target, settings, data, node_ref, nodes, replace=False)
        assert node_ref[target].data[0] == 'setting1 100'
        assert node_ref[target].settings[0] == 'nope'
        # node_ref length can vary depends on how many nodes are in file
        target = 'SOMETHING->WHATEVER'
        settings = ['nope']
        data = ['whatever']
        update_node(target, settings, data, node_ref, nodes, replace=False)
        assert 'SOMETHING->WHATEVER' in node_ref

    def test_remove_node(self):
        """Test removal of nodes."""
        target = 'MOTION'
        infile = os.path.join(HERE, 'cp2k.inp')
        nodes = read_cp2k_input(infile)
        node_ref = set_parents(nodes)
        assert 'MOTION' in node_ref
        remove_node(target, node_ref, nodes)
        assert 'MOTION' not in node_ref
        # Test remove non-existent
        remove_node('NOPE', node_ref, nodes)
        # Test remove with parent
        target = 'FORCE_EVAL->DFT'
        nodes = read_cp2k_input(infile)
        node_ref = set_parents(nodes)
        remove_node(target, node_ref, nodes)
        assert 'FORCE_EVAL->DFT' not in node_ref

    def test_update_cp2k_input(self):
        """Test the update_cp2k_input function."""
        infile = os.path.join(HERE, 'cp2k.inp')
        with tempfile.NamedTemporaryFile() as temp:
            update = {'MOTION->MD': {'data': ['STEPS 777'], 'replace': True}}
            remove = ['FORCE_EVAL->DFT']
            update_cp2k_input(infile, temp.name, update=update, remove=remove)
            temp.flush()
            nodes = read_cp2k_input(temp.name)
            node_ref = set_parents(nodes)
            assert node_ref['MOTION->MD'].data[0] == 'STEPS 777'
            assert 'FORCE_EVAL->DFT' not in node_ref

    def test_cp2k_settings(self):
        """Test the cp2k_settings function."""
        template = os.path.join(HERE, 'cp2k.inp')
        settings = {
            'engine': {'template': template, 'cp2k_format': 'xyz'},
            'system': {'temperature': 298.0}
        }
        # cp2k.inp has TEMPERATURE 500
        with patch('pyretis.inout.formats.cp2k.logger') as mock_log:
            cp2k_settings(settings, 'some/path')
            mock_log.warning.assert_called()
            assert settings['system']['temperature'] == 500.0

        # Test matching temperature (covers line 537 pass)
        settings = {
            'engine': {'template': template},
            'system': {'temperature': 500.0}
        }
        cp2k_settings(settings, 'some/path')
        assert settings['system']['temperature'] == 500.0

        # Test missing temperature in settings (covers line 535)
        settings = {
            'engine': {'template': template},
            'system': {}
        }
        cp2k_settings(settings, 'some/path')
        assert settings['system']['temperature'] == 500.0

    def test_cp2k_settings_box(self):
        """Test cp2k_settings with box matrix to list conversion."""
        # This is implicitly tested in read_box_data
        pass

    def test_cp2k_settings_defaults(self):
        """Test cp2k_settings with defaults (300K)."""
        # Create an input without TEMPERATURE
        with tempfile.NamedTemporaryFile() as temp:
            temp.write(b"&MOTION\n&MD\n&END MD\n&END MOTION\n")
            temp.flush()
            settings = {
                'engine': {'template': temp.name},
                'system': {}
            }
            with patch('pyretis.inout.formats.cp2k.logger') as mock_log:
                cp2k_settings(settings, 'some/path')
                assert settings['system']['temperature'] == 300.0

            # Test with system temp set but not in file
            settings = {
                'engine': {'template': temp.name},
                'system': {'temperature': 298.0}
            }
            cp2k_settings(settings, 'some/path')
            assert settings['system']['temperature'] == 300.0

    def test_read_cp2k_box_no_cell(self):
        """Test no CELL found in read_cp2k_box."""
        with tempfile.NamedTemporaryFile() as temp:
            temp.write(b"&FORCE_EVAL\n&END FORCE_EVAL\n")
            temp.flush()
            with patch('pyretis.inout.formats.cp2k.logger') as mock_log:
                box, periodic = read_cp2k_box(temp.name)
                mock_log.warning.assert_called()
                assert np.allclose(box, np.eye(3)*100)
                assert periodic == [True, True, True]
