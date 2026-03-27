# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A test module for the energy analysis."""
import logging
import pytest
import os
import pickle
import numpy as np
from pyretis.analysis.energy_analysis import analyse_energies
from pyretis.inout.settings import SECTIONS
from pyretis.inout.formats.energy import EnergyFile

logging.disable(logging.CRITICAL)
HERE = os.path.abspath(os.path.dirname(__file__))


class EnergyTest:
    """Test that we can analyse energies."""

    def test_energy_analysis(self):
        """Test the energy analysis."""
        filename = os.path.join(HERE, 'energy.txt')
        data = None
        with EnergyFile(filename, 'r') as efile:
            data = efile.load()
        settings = {
            'particles': {'npart': 1},
            'system': {'dimensions': 1, 'beta': 1, 'temperature': 1},
        }
        settings['analysis'] = SECTIONS['analysis']

        correct_file = os.path.join(HERE, 'energy-results.dat')
        for i in data:
            results = analyse_energies(i['data'], settings)
            break
        correct_data = {}
        with open(correct_file, 'rb') as infile:
            correct_data = pickle.load(infile)
        for key0 in ('time', 'temp', 'etot', 'vpot', 'ekin'):
            assert key0 in correct_data
            assert key0 in results
            for key, val in correct_data[key0].items():
                for i, j in zip(val, results[key0][key]):
                    assert np.allclose(i, j)
