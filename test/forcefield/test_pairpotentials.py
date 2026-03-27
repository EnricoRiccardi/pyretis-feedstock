# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A test of the pairpotential module."""
import logging
import pytest
import numpy as np
from pyretis.forcefield.potentials.pairpotentials.pairpotential import (
    mixing_parameters,
    generate_pair_interactions,
    _check_pair_parameters,
)
logging.disable(logging.CRITICAL)


class TestMixingParameters:
    """Test that we can set up for mixing."""

    def test_mixing_geometric(self):
        """Test the geometric mixing."""
        epsilon_i, epsilon_j = 2.0, 3.0
        sigma_i, sigma_j = 5.0, 6.0
        rcut_i, rcut_j = 10.0, 11.0
        mixed = mixing_parameters([epsilon_i, epsilon_j],
                                  [sigma_i, sigma_j],
                                  [rcut_i, rcut_j],
                                  mixing='geometric')
        for i, j in zip(mixed, ((epsilon_i * epsilon_j),
                                (sigma_i * sigma_j),
                                (rcut_i * rcut_j))):
            assert i == np.sqrt(j)

    def test_mixing_arithmetic(self):
        """Test arithmetic mixing."""
        epsilon_i, epsilon_j = 2.0, 3.0
        sigma_i, sigma_j = 5.0, 6.0
        rcut_i, rcut_j = 10.0, 11.0
        correct = (np.sqrt(epsilon_i * epsilon_j),
                   0.5 * (sigma_i + sigma_j),
                   0.5 * (rcut_i + rcut_j))
        mixed = mixing_parameters([epsilon_i, epsilon_j],
                                  [sigma_i, sigma_j],
                                  [rcut_i, rcut_j],
                                  mixing='arithmetic')
        for i, j in zip(mixed, correct):
            assert i == j

    def test_mixing_sixthpower(self):
        """Test 6th power mixing."""
        epsilon_i, epsilon_j = 2.0, 3.0
        sigma_i, sigma_j = 5.0, 6.0
        rcut_i, rcut_j = 10.0, 11.0
        si3 = sigma_i**3
        si6 = si3**2
        sj3 = sigma_j**3
        sj6 = sj3**2
        avgs6 = 0.5 * (si6 + sj6)
        mixed = mixing_parameters([epsilon_i, epsilon_j],
                                  [sigma_i, sigma_j],
                                  [rcut_i, rcut_j],
                                  mixing='sixthpower')
        correct = (np.sqrt(epsilon_i * epsilon_j) * si3 * sj3 / avgs6,
                   avgs6**(1.0 / 6.0),
                   (0.5 * (rcut_i**6 + rcut_j**6))**(1.0 / 6.0))
        for i, j in zip(mixed, correct):
            assert i == j

    def test_mixing_missing(self, caplog):
        """Test mixing when we supply something not recognized."""
        logging.disable(logging.INFO)
        mod = 'pyretis.forcefield.potentials.pairpotentials.pairpotential'
        with caplog.at_level(logging.WARNING, logger=mod):
            mixed = mixing_parameters([11., 14], [12., 15], [13., 16.],
                                      mixing='Life On Mars?')
        logging.disable(logging.CRITICAL)
        for i in mixed:
            assert i == 0.0


class TestGeneratePairInteraction:
    """Test generation of pair parameters."""

    def test_generate_full(self):
        """Test generation of full parameters."""
        parameters = {0: {'sigma': 1.0, 'epsilon': 1.5, 'rcut': 2.},
                      1: {'sigma': 2.0, 'epsilon': 2.5, 'rcut': 3.},
                      2: {'sigma': 3.0, 'epsilon': 3.5, 'rcut': 4.}}
        mix = generate_pair_interactions(parameters, 'arithmetic')
        correct_mix = {
            (0, 0): {'rcut': 2.0, 'epsilon': 1.5, 'sigma': 1.0},
            (2, 2): {'rcut': 4.0, 'epsilon': 3.5, 'sigma': 3.0},
            (1, 1): {'rcut': 3.0, 'epsilon': 2.5, 'sigma': 2.0},
            (0, 1): {'rcut': 2.5, 'epsilon': 1.9364916731037085, 'sigma': 1.5},
            (0, 2): {'rcut': 3.0, 'epsilon': 2.2912878474779199, 'sigma': 2.0},
            (1, 2): {'rcut': 3.5, 'epsilon': 2.9580398915498081, 'sigma': 2.5},
        }
        for key, val in correct_mix.items():
            for key2, val2 in val.items():
                assert val2 == mix[key][key2]
        for key1, key2 in zip(((1, 0), (2, 0), (2, 1)),
                              ((0, 1), (0, 2), (1, 2))):
            for key, val in mix[key1].items():
                assert val == mix[key2][key]

    def test_generate_pair(self):
        """Testing generate where we specify some cross interactions."""
        parameters = [
            {0: {'sigma': 1.0, 'epsilon': 1.5, 'rcut': 2.},
             1: {'sigma': 2.0, 'epsilon': 2.5, 'rcut': 3.},
             (0, 1): {'sigma': 1.0, 'epsilon': 1.0, 'rcut': 100.}},
            {0: {'sigma': 1.0, 'epsilon': 1.5, 'rcut': 2.},
             1: {'sigma': 2.0, 'epsilon': 2.5, 'rcut': 3.},
             (1, 0): {'sigma': 1.0, 'epsilon': 1.0, 'rcut': 100.}},
        ]
        for param in parameters:
            mix = generate_pair_interactions(param, 'arithmetic')
            correct_mix = {
                (0, 0): param[0],
                (1, 1): param[1],
            }
            if (0, 1) in param:
                correct_mix[(0, 1)] = param[(0, 1)]
                correct_mix[(1, 0)] = param[(0, 1)]
            elif (1, 0) in param:
                correct_mix[(0, 1)] = param[(1, 0)]
                correct_mix[(1, 0)] = param[(1, 0)]
            for key, val in correct_mix.items():
                for key2, val2 in val.items():
                    assert val2 == mix[key][key2]

    def test_check_pair_params(self):
        """Test the check pair params."""
        parameters = {
            0: {'sigma': 1.0, 'epsilon': 1.5, 'rcut': 2.},
            1: {'sigma': 2.0, 'rcut': 3.}
        }
        assert 'epsilon' not in parameters[1]
        _check_pair_parameters(parameters)
        assert 'epsilon' in parameters[1]
        assert parameters[1]['epsilon'] == 0.0
        # Missing rcut:
        parameters = {0: {'sigma': 1.0, 'epsilon': 1.5, 'factor': 123.}}
        assert 'rcut' not in parameters[0]
        _check_pair_parameters(parameters)
        assert 'rcut' in parameters[0]
        assert parameters[0]['rcut'] == 123.
