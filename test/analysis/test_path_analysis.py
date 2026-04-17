# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A test module for path ensemble analysis."""
import copy
import logging
import os
import pytest
import numpy as np
from pyretis.analysis.path_analysis import (
    analyse_path_ensemble,
    analyse_repptis_ensemble,
    match_probabilities,
    perm_calculations,
    retis_flux,
    retis_rate,
    skip_paths
)
from pyretis.core.path import Path
from pyretis.core.system import System
from pyretis.inout.analysisio.analysisio import (
    repptis_running_pcross_analysis,
)
from pyretis.inout.formats.pathensemble import PathEnsembleFile
from pyretis.inout.settings import SECTIONS
from .help import turn_on_logging
logging.disable(logging.CRITICAL)


HERE = os.path.abspath(os.path.dirname(__file__))


def add_some_paths(ensemble, raw_data):
    """Add some fake path to the given ensemble."""
    interfaces = ensemble.interfaces
    for _, path in enumerate(raw_data.get_paths()):
        dummy_path = Path(None, maxlen=None, time_origin=0)
        dummy_path.generated = path['generated']
        # Add fake points for the interfaces:
        start, _, end = path['interface']
        order = ([((path['ordermax'][0] + path['ordermin'][0]) * 0.5, 0)] *
                 path['length'])
        for i in (start, end):
            if i == 'L':
                order.append([interfaces[0] - 0.1])
            elif start == 'R':
                order.append([interfaces[2] + 0.1])
            else:
                order.append([0.5*(interfaces[0] + interfaces[2])])
        order[path['ordermax'][1]] = path['ordermax']
        order[path['ordermin'][1]] = path['ordermin']
        for i in order:
            phasepoint = System()
            phasepoint.order = i
            dummy_path.append(phasepoint)
        ensemble.add_path_data(dummy_path, path['status'])


class TestAnalysePathEnsemble:
    """Test that we run the analysis for PathEnsemble results."""

    def test_no_accept_ens_1(self, caplog):
        """Test when no accepted paths are available in 0+."""
        ens = {
            'ensemble_number': 1,
            'interfaces': [-float('inf'), -0.9, -0.9],
            'file': 'pathensemble_no_acc.txt',
            'test': {'fluxlength': (0, 1541.7987987987988)},
            }
        settings = {'analysis': SECTIONS['analysis']}
        settings['tis'] = {'detect': -0.9}
        filename = os.path.join(HERE, ens['file'])
        raw_data = PathEnsembleFile(filename, 'r', ensemble_settings=ens)
        with caplog.at_level(logging.WARNING):
            result = analyse_path_ensemble(raw_data, settings)
        assert 'pathlength' not in result
        assert 'No accepted paths' in caplog.text

    def test_no_accept_ens_0(self, caplog):
        """Test when no accepted paths are available in 0-."""
        ens = {
            'ensemble_number': 0,
            'interfaces': [-float('inf'), -0.9, -0.9],
            'file': 'pathensemble_no_acc.txt',
            'test': {'fluxlength': (0, 1541.7987987987988)},
            }
        settings = {'analysis': SECTIONS['analysis']}
        settings['tis'] = {'detect': -0.9}
        filename = os.path.join(HERE, ens['file'])
        raw_data = PathEnsembleFile(filename, 'r', ensemble_settings=ens)
        with caplog.at_level(logging.WARNING):
            result = analyse_path_ensemble(raw_data, settings)
        assert 'pathlength' not in result
        assert 'No accepted paths' in caplog.text

    def test_path_analysis(self):
        """Test the path ensemble analysis."""
        ensembles = [
            {
                'ensemble_number': 0,
                'interfaces': [-float('inf'), -0.9, -0.9],
                'file': 'pathensemble000.txt',
                'test': {'fluxlength': (0, 1541.7987987987988)},
            },
            {
                'ensemble_number': 1,
                'interfaces': [-0.9, -0.9, 1],
                'file': 'pathensemble001.txt',
                'test': {'prun': (-1, 223/996)},
            },
            {
                'ensemble_number': 2,
                'interfaces': [-0.9, -0.8, 1],
                'file': 'pathensemble002.txt',
                'test': {'prun': (-1, 145/996)},
            },
            {
                'ensemble_number': 3,
                'interfaces': [-0.9, -0.7, 1],
                'file': 'pathensemble003.txt',
                'test': {'prun': (-1, 125/999)},
            },
            {
                'ensemble_number': 4,
                'interfaces': [-0.9, -0.71, 1],
                'file': 'pathensemble004.txt',
                'test': {'prun': (-1, 125/966)},
            },
        ]
        settings = {'analysis': SECTIONS['analysis']}
        results = []
        detects = [-0.9, -0.8, -0.7, -0.6, -0.6]
        raw_data = None
        val, key = None, None
        for ens, dec in zip(ensembles, detects):
            filename = os.path.join(HERE, ens['file'])
            settings['tis'] = {'detect': dec}
            raw_data = PathEnsembleFile(filename, 'r', ensemble_settings=ens)
            res = analyse_path_ensemble(raw_data, settings)
            for key, val in ens['test'].items():
                assert val[-1] == pytest.approx(res[key][val[0]])
            results.append(analyse_path_ensemble(raw_data, settings))
        match = match_probabilities(results[1:], detects[1:], settings)
        assert match['prob'] == pytest.approx(0.0005277540804156307)
        flux = retis_flux(results[0], results[1], 0.002)
        assert flux[0] == pytest.approx(0.26572079970779655)
        assert flux[1] == pytest.approx(0.02422284635774688)
        rate = retis_rate(match['prob'], match['relerror'], flux[0], flux[1])
        assert rate[0] == pytest.approx(0.00014023523629709417)
        assert rate[1] == pytest.approx(0.4734317724627789)

        # Re-check the last one:
        settings['analysis']['maxblock'] = 0
        res = analyse_path_ensemble(raw_data, settings)
        assert val[-1] == pytest.approx(res[key][val[0]])

    def test_skipping(self, caplog):
        """Test the skip_lines funcionality."""
        ens = {
            'ensemble_number': 0,
            'interfaces': [-float('inf'), -0.1, -0.1],
            'file': 'permeability/pathensemble.txt',
            }
        settings = {'analysis': SECTIONS['analysis']}
        # Skip more trajectories than there are
        # Make a copy to prevent overwrite
        settings = copy.deepcopy(settings)
        settings['tis'] = {'detect': -0.1}
        settings['analysis']['skip'] = 10000
        filename = os.path.join(HERE, ens['file'])
        raw_data = PathEnsembleFile(filename, 'r', ensemble_settings=ens)
        with caplog.at_level(logging.WARNING):
            result = analyse_path_ensemble(raw_data, settings)
        assert 'pathlength' not in result
        assert 'No accepted paths' in caplog.text
        # Make sure we can skip all trajectories
        settings['analysis']['skip'] = 49
        _ = analyse_path_ensemble(raw_data, settings)
        # Test more skipping than possible (Hard to hit normally)
        with pytest.raises(RuntimeError, match="Skipping more trajs than"):
            skip_paths([1, 1], 2, 3)
        # Test exact skipping of all
        weights, nacc = skip_paths([1, 2], 3, 3)
        assert weights == [0, 0]
        assert nacc == 1

        # See if this also works in ensemble 1
        ens['ensemble_number'] = 1
        raw_data = PathEnsembleFile(filename, 'r', ensemble_settings=ens)
        settings['analysis']['skip'] = 10000
        caplog.clear()
        with caplog.at_level(logging.WARNING):
            result = analyse_path_ensemble(raw_data, settings)
        assert 'pathlength' not in result
        assert 'No accepted paths' in caplog.text

    def test_permeability_path_analysis_interfaces(self, caplog):
        """Test the permeability calculations."""
        ens = {
            'ensemble_number': 0,
            'interfaces': [-float('inf'), -0.1, -0.1],
            'file': 'permeability/pathensemble.txt',
            }

        settings = dict(analysis=SECTIONS['analysis'],
                        simulation={'permeability': True})

        settings['analysis']['tau_ref_bin'] = [-0.175, -0.125]
        settings['tis'] = {'detect': -0.1}
        filename = os.path.join(HERE, ens['file'])
        raw_data = PathEnsembleFile(filename, 'r', ensemble_settings=ens)
        # Test no zero left
        result = analyse_path_ensemble(raw_data, settings)
        # Check if zero left is not overwritten
        assert result['interfaces'][0] == -float('inf')
        # Check if hr is properly calculated
        assert result['xi'][-1] == 30/49
        # The error is about 3%
        assert 0.03 < result['xierror'][4] < 0.04
        # Check tau/dz
        assert result['tau'][-1] == (8/49)/(-0.125+0.175)
        # Tau error is about 19%
        assert 0.19 < result['tauerror'][4] < 0.20

        # Test with zero left
        settings['simulation']['zero_left'] = -0.2
        result = analyse_path_ensemble(raw_data, settings)

        # Assert zero left overwritten and no data change
        assert result['interfaces'][0] == -0.2
        # Check if hr is properly calculated
        assert result['xi'][-1] == 30/49
        # The error is about 3%
        assert 0.03 < result['xierror'][4] < 0.04
        # Check tau/dz
        assert result['tau'][-1] == (8/49)/(-0.125+0.175)
        # Tau error is about 19%
        assert 0.19 < result['tauerror'][4] < 0.20

        # Now we are going to use the miss aligned order.txt (cycle 0 => 10)
        # It should not matter for the result as it is part of the init
        op_dir = filename.rsplit('/', 1)[0]
        correct = (op_dir+'/order.txt', op_dir+'/order_correct.txt')
        mis = (op_dir+'/order.txt', op_dir+'/order_misaligned.txt')
        os.rename(*correct)
        os.rename(*mis[::-1])
        raw_data2 = PathEnsembleFile(filename, 'r', ensemble_settings=ens)
        with turn_on_logging():
            with caplog.at_level(logging.WARNING,
                                 logger='pyretis.analysis.path_analysis'):
                result2 = analyse_path_ensemble(raw_data2, settings)
        assert any('align' in r.message.lower() for r in caplog.records)

        # Move files back
        os.rename(*mis)
        os.rename(*correct[::-1])

        assert (result['xi'] == result2['xi']).all()
        assert (result['tau'] == result2['tau']).all()

        # Now assert the other calculations
        result0 = {'out': result}
        # Make_fake_cros
        pcross = 0.5
        perr = 0.5
        xi, tau, perm, xi_err, tau_err, perm_err = perm_calculations(
            result0, pcross, perr)

        assert xi == result['xi'][-1]
        assert tau == result['tau'][-1]
        assert xi_err == result['xierror'][4]
        assert tau_err == result['tauerror'][4]
        # Check permeation formula and error propagation
        assert perm == result['xi'][-1]*(1/result['tau'][-1])*pcross
        assert perm_err == (xi_err**2+tau_err**2+perr**2)**(1/2)

    def test_path_analysis_pptis(self):
        """Test the path ensemble analysis for (RE)PPTIS simulations."""
        pptisdir = "pptis"
        ensembles = [
            {
                'ensemble_number': 0,
                'interfaces': [0.1, 0.2, 0.2],
                'file': os.path.join(pptisdir, 'pathensemble0.txt'),
                'test': {'fluxlength': (0, 91.98)},
            },
            {
                'ensemble_number': 1,
                'interfaces': [0.2, 0.2, 0.325],
                'file': os.path.join(pptisdir, 'pathensemble1.txt'),
                'test': {'prun_sl': (-1, 0.07526881720430108),
                         'prun_sr': (-1, 1.)},
            },
            {
                'ensemble_number': 2,
                'interfaces': [0.2, 0.325, 0.55],
                'file': os.path.join(pptisdir, 'pathensemble2.txt'),
                'test': {'prun_sl': (-1, 0.010638297872340425),
                         'prun_sr': (-1, 1.)},
            },
            {
                'ensemble_number': 3,
                'interfaces': [0.325, 0.55, 0.69],
                'file': os.path.join(pptisdir, 'pathensemble3.txt'),
                'test': {'prun_sl': (-1, 0.6451612903225806),
                         'prun_sr': (-1, 0.5507246376811594)},
            },
            {
                'ensemble_number': 4,
                'interfaces': [0.55, 0.69, 0.75],
                'file': os.path.join(pptisdir, 'pathensemble4.txt'),
                'test': {'prun_sl': (-1, 0.4807692307692308),
                         'prun_sr': (-1, 0.723404255319149)},
            },
            {
                'ensemble_number': 5,
                'interfaces': [0.69, 0.75, 0.9],
                'file': os.path.join(pptisdir, 'pathensemble5.txt'),
                'test': {'prun_sl': (-1, 0.19540229885057472),
                         'prun_sr': (-1, 0.8461538461538461)},
            }
        ]
        settings = {'analysis': SECTIONS['analysis']}
        detects = [0.2, 0.325, 0.55, 0.69, 0.75, 0.9]
        p_types = []
        results = []
        for ens, dec in zip(ensembles, detects):
            filename = os.path.join(HERE, ens['file'])
            settings['tis'] = {'detect': dec}
            raw_data = PathEnsembleFile(filename, 'r', ensemble_settings=ens)

            res = analyse_repptis_ensemble(raw_data, settings)
            for key, val in ens['test'].items():
                assert val[-1] == pytest.approx(res[key][val[0]])
            if ens['ensemble_number'] > 0:
                p_types.append(res['ptypes'])
            results.append(res)
        pcrossdict = repptis_running_pcross_analysis(p_types)
        assert pcrossdict['overall-prun'][-1] == pytest.approx(
            0.00010317722881789884)
        assert pcrossdict['overall-error'][-1] == pytest.approx(
            0.4873103759019344)
        flux = retis_flux(results[0], results[1], 0.01)
        assert flux[0] == pytest.approx(0.7346382685689098)
        assert flux[1] == pytest.approx(0.19580016239096693)
        rate = retis_rate(pcrossdict['overall-prun'][-1],
                          pcrossdict['overall-error'][-1],
                          flux[0], flux[1])
        assert rate[0] == pytest.approx(7.579794073451942e-05)
        assert rate[1] == pytest.approx(0.5251753098290262)

    def test_path_analysis_pptis_skiplines(self):
        """Test skiplines path ensemble analysis for REPPTIS simulations."""
        # we take the same pathensemble files as the previous test, but skip
        # 3 of the 100 lines.
        pptisdir = "pptis"
        ensembles = [
            {
                'ensemble_number': 0,
                'interfaces': [0.1, 0.2, 0.2],
                'file': os.path.join(pptisdir, 'pathensemble0.txt'),
                'test': {'fluxlength': (0, 91.21649484536083)},
            },
            {
                'ensemble_number': 1,
                'interfaces': [0.2, 0.2, 0.325],
                'file': os.path.join(pptisdir, 'pathensemble1.txt'),
                'test': {'prun_sl': (-1, 0.07692307692307693),
                         'prun_sr': (-1, 1.)},
            },
            {
                'ensemble_number': 2,
                'interfaces': [0.2, 0.325, 0.55],
                'file': os.path.join(pptisdir, 'pathensemble2.txt'),
                'test': {'prun_sl': (-1, 0.01098901098901099),
                         'prun_sr': (-1, 1.)},
            },
            {
                'ensemble_number': 3,
                'interfaces': [0.325, 0.55, 0.69],
                'file': os.path.join(pptisdir, 'pathensemble3.txt'),
                'test': {'prun_sl': (-1, 0.6333333333333333),
                         'prun_sr': (-1, 0.5757575757575758)},
            },
            {
                'ensemble_number': 4,
                'interfaces': [0.55, 0.69, 0.75],
                'file': os.path.join(pptisdir, 'pathensemble4.txt'),
                'test': {'prun_sl': (-1, 0.47058823529411764),
                         'prun_sr': (-1, 0.717391304347826)},
            },
            {
                'ensemble_number': 5,
                'interfaces': [0.69, 0.75, 0.9],
                'file': os.path.join(pptisdir, 'pathensemble5.txt'),
                'test': {'prun_sl': (-1, 0.21518987341772153),
                         'prun_sr': (-1, 0.8333333333333334)},
            }
        ]
        settings = {'analysis': SECTIONS['analysis'].copy()}
        settings['analysis']['skip'] = 3
        detects = [0.2, 0.325, 0.55, 0.69, 0.75, 0.9]
        p_types = []
        results = []
        for ens, dec in zip(ensembles, detects):
            filename = os.path.join(HERE, ens['file'])
            settings['tis'] = {'detect': dec}
            raw_data = PathEnsembleFile(filename, 'r', ensemble_settings=ens)
            res = analyse_repptis_ensemble(raw_data, settings)
            for key, val in ens['test'].items():
                assert val[-1] == pytest.approx(res[key][val[0]])
            if ens['ensemble_number'] > 0:
                p_types.append(res['ptypes'])
            results.append(res)
        pcrossdict = repptis_running_pcross_analysis(p_types)
        assert pcrossdict['overall-prun'][-1] == pytest.approx(
            0.00011041092698379775)
        assert pcrossdict['overall-error'][-1] == pytest.approx(
            0.5460014951716033)
        flux = retis_flux(results[0], results[1], 0.01)
        assert flux[0] == pytest.approx(0.7391039317281317)
        assert flux[1] == pytest.approx(0.21291472586483318)
        rate = retis_rate(pcrossdict['overall-prun'][-1],
                          pcrossdict['overall-error'][-1],
                          flux[0], flux[1])
        assert rate[0] == pytest.approx(8.160515023947258e-05)
        assert rate[1] == pytest.approx(0.5860463405053592)

    def test_path_analysis_pptis_local_nans(self):
        """Test REPPTIS analysis when some local probs are NaN.

        We still want a PDF report to be written about the other
        ensembles, so we make sure that NaNs are outputted without
        a crash occurring.

        """
        pptisdir = "pptis"
        ensembles = [
            {
                'ensemble_number': 0,
                'interfaces': [0.1, 0.2, 0.2],
                'file': os.path.join(pptisdir, 'pathensemble0.txt'),
                'test': {'fluxlength': (0, 91.21649484536083)},
            },
            {
                'ensemble_number': 1,
                'interfaces': [0.2, 0.2, 0.325],
                'file': os.path.join(pptisdir, 'pathensemble1.txt'),
                'test': {'prun_sl': (-1, 0.07692307692307693),
                         'prun_sr': (-1, 1.)},
            },
            {
                'ensemble_number': 2,
                'interfaces': [0.2, 0.325, 0.55],
                'file': os.path.join(pptisdir, 'pathensemble2_noRacc.txt'),
                'test': {'prun_sl': (-1, 0.010309278350515464),
                         'prun_sr': (-1, 0.0)},
            },
            {
                'ensemble_number': 3,
                'interfaces': [0.325, 0.55, 0.69],
                'file': os.path.join(pptisdir, 'pathensemble3_noLacc.txt'),
                'test': {'prun_sl': (-1, 0.0),
                         'prun_sr': (-1, 0.4742268041237113)},
            },
            {
                'ensemble_number': 4,
                'interfaces': [0.55, 0.69, 0.75],
                'file': os.path.join(pptisdir, 'pathensemble4.txt'),
                'test': {'prun_sl': (-1, 0.47058823529411764),
                         'prun_sr': (-1, 0.717391304347826)},
            },
            {
                'ensemble_number': 5,
                'interfaces': [0.69, 0.75, 0.9],
                'file': os.path.join(pptisdir, 'pathensemble5.txt'),
                'test': {'prun_sl': (-1, 0.21518987341772153),
                         'prun_sr': (-1, 0.8333333333333334)},
            }
        ]
        settings = {'analysis': SECTIONS['analysis'].copy()}
        settings['analysis']['skip'] = 3
        detects = [0.2, 0.325, 0.55, 0.69, 0.75, 0.9]
        p_types = []
        results = []
        for ens, dec in zip(ensembles, detects):
            filename = os.path.join(HERE, ens['file'])
            settings['tis'] = {'detect': dec}
            raw_data = PathEnsembleFile(filename, 'r', ensemble_settings=ens)
            res = analyse_repptis_ensemble(raw_data, settings)
            for key, val in ens['test'].items():
                assert val[-1] == pytest.approx(res[key][val[0]])
            if ens['ensemble_number'] > 0:
                p_types.append(res['ptypes'])
            results.append(res)
        pcrossdict = repptis_running_pcross_analysis(p_types)
        assert np.isnan(pcrossdict['overall-prun'][-1])
        assert np.isnan(pcrossdict['overall-error'][-1])
        flux = retis_flux(results[0], results[1], 0.01)
        assert flux[0] == pytest.approx(0.7391039317281317)
        assert flux[1] == pytest.approx(0.21291472586483318)
        rate = retis_rate(pcrossdict['overall-prun'][-1],
                          pcrossdict['overall-error'][-1],
                          flux[0], flux[1])
        assert np.isnan(rate[0])
        assert np.isnan(rate[1])
