# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A test module for flux analysis."""
import logging
import pytest
import warnings
import os
import pickle
import sys
import numpy as np
from pyretis.analysis import analyse_md_flux
from pyretis.analysis.flux_analysis import analyse_flux, find_crossings
from pyretis.inout.formats.cross import CrossFile
from pyretis.inout.formats.order import OrderFile
from pyretis.inout.settings import SECTIONS

logging.disable(logging.CRITICAL)
HERE = os.path.abspath(os.path.dirname(__file__))


class TestFluxAnalysis:
    """Test that we can analyse for initial flux."""

    def test_flux_analysis(self):
        """Test the flux analysis."""
        filename = os.path.join(HERE, 'cross.txt')
        data = None
        with CrossFile(filename, 'r') as crossfile:
            data = crossfile.load()
        settings = {
            'simulation': {'endcycle': 250000, 'interfaces': [-0.9]},
            'engine': {'timestep': 0.002},
            'particles': {'npart': 1},
            'system': {'dimensions': 1,
                       'temperature': 666,
                       'beta': 1},
        }
        settings['analysis'] = SECTIONS['analysis']
        # What if empty?
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            results_all = analyse_md_flux(crossdata=[],
                                          energydata=[],
                                          orderdata=np.zeros(shape=(5, 2)),
                                          settings=settings)

        results = results_all['flux']
        assert not results['eff_cross']

        correct_file = os.path.join(HERE, 'flux-results.dat')
        for i in data:
            results = analyse_flux(i['data'], settings)
            break
        correct_data = {}
        with open(correct_file, 'rb') as infile:
            correct_data = pickle.load(infile)
        for key0 in ('eff_cross', 'ncross', 'neffcross',
                     'flux', 'runflux', 'interfaces', 'cross_time',
                     'neffc/nc', 'pMD', '1-p', 'teffMD', 'corrMD'):
            assert key0 in correct_data
            assert key0 in results
            for i, j in zip(correct_data[key0], results[key0]):
                try:
                    assert np.allclose(np.asanyarray(i), np.asanyarray(j))
                except Exception as e:
                    sys.stderr.write(f"\nFAILED on key: {key0}\n")
                    sys.stderr.write(f"i: {i} ({type(i)})\n")
                    sys.stderr.write(f"j: {j} ({type(j)})\n")
                    raise e
        assert correct_data['totalcycle'] == results['totalcycle']
        for key, val in correct_data['times'].items():
            assert val == pytest.approx(results['times'][key])
        for data1, data2 in zip(correct_data['errflux'],
                                results['errflux']):
            for i, j in zip(data1, data2):
                assert np.allclose(np.asanyarray(i), np.asanyarray(j))

    def test_flux_analysis_B(self):
        """Test the flux analysis considers also flux from state B."""
        fluxdata = [(592, 2, -1), (1171, 2, 1), (1642, 1, -1), (2279, 1, 1),
                    (2851, 3, -1), (3497, 2, 1), (3937, 2, -1), (4838, 2, -1),
                    (4967, 3, -1), (5288, 3, 1), (6666, 2, -1), (6727, 2, -1),
                    (6828, 3, -1), (6869, 3, 1), (6971, 2, -1), (9294, 3, 1)]
        settings = {
            'simulation': {'endcycle': 9999, 'interfaces': [-0.9, -0.8, -0.7]},
            'engine': {'timestep': 0.002},
            'particles': {'npart': 1},
            'system': {'dimensions': 1,
                       'temperature': 666,
                       'beta': 1},
        }
        settings['analysis'] = SECTIONS['analysis']
        results_all = analyse_flux(fluxdata, settings)
        assert results_all['eff_cross'][0][0] == (637, 2279)
        assert results_all['teffMD'][0] == 9999.0
        assert results_all['corrMD'][0] == 3.334111370456819

    def test_crossing_calculation(self):
        """Test calculation of crossings from order parameter data."""
        filename = os.path.join(HERE, 'order-calculate-crossings.txt')
        data = None
        with OrderFile(filename, 'r') as orderfile:
            data = orderfile.load()
        orderp = next(data)['data'][:, 1]
        cross = find_crossings(orderp, [-0.9, -0.85, -0.70])

        filename_cross = os.path.join(HERE, 'cross-correct.txt')
        cross_correct = None
        with CrossFile(filename_cross, 'r') as crossfile:
            cross_correct = next(crossfile.load())['data']
        assert len(cross) == len(cross_correct)
        for i, j in zip(cross, cross_correct):
            assert i[0] == j[0]
            assert i[1] == j[1] - 1
            assert i[2] in ('+', '-')
            if i[2] == '-':
                assert j[2] == -1
            elif i[2] == '+':
                assert j[2] == 1

    def test_flux_analysis_line186(self):
        """Test the logic path for line 186 (simulation ends in A)."""
        fluxdata = [(10, 1, -1)]
        settings = {
            'simulation': {'endcycle': 100, 'interfaces': [-0.9, -0.8, -0.7]},
            'engine': {'timestep': 0.1, 'subcycles': 1},
            'analysis': {'skipcross': 1, 'maxblock': 1, 'blockskip': 1}
        }
        results = analyse_flux(fluxdata, settings)
        assert results['times']['A'] == 90  # 100 - 10
