# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Tests to raise coverage for report_md, report_path, and analysisio.

These tests exercise branches not covered by the main test suite:
- report_md: fmt='tex' and fmt='txt' branches in table helpers
- report_path: fmt='tex' branch in table helpers, generate_report functions
- analysisio: make-tis-files task, recursive_block_analysis short list,
              get_global_probz NaN/zero branch
"""
import os
import numpy as np
import pytest
from pyretis.inout.report.report_path import (
    _table_interface,
    _table_interface0,
    _table_probability,
    _table_probability_repptis,
    _table_efficiencies,
    _table_summary,
    generate_report_tis,
    generate_report_retis,
    generate_report_retis0,
)
from pyretis.inout.report.report_md import (
    _table_md_efficiency,
    _table_md_flux_cycles,
    _table_md_flux,
)
from pyretis.inout.analysisio.analysisio import (
    recursive_block_analysis,
    get_global_probz,
    get_path_simulation_files,
)

HERE = os.path.abspath(os.path.dirname(__file__))


# ---------------------------------------------------------------------------
# Minimal synthetic result dict for report_path table helpers
# ---------------------------------------------------------------------------

def _make_result(ensemble='001', interfaces=(-0.9, -0.9, 1.0)):
    """Return a minimal result dict for table helpers in report_path."""
    return {
        'out': {
            'ensemble': ensemble,
            'interfaces': list(interfaces),
            'detect': interfaces[1] if len(interfaces) > 1 else None,
            'prun': [0.1, 0.2, 0.3],
            'blockerror': [list(range(5)), 0, 0.01, [0.01] * 5, 0.05, 0, 0.1],
            'tis-cycles': 100,
            'efficiency': [0.5, 0.6],
            'pathlength': [
                ([0.1, 0.2], [1, 2], (5.0, 1.0))
            ],
            'prun_sl': [0.1, 0.2, 0.3],
            'prun_sr': [0.05, 0.1, 0.15],
            'blockerror_sl': [list(range(5)), 0, 0.01, [0.01] * 5, 0.05],
            'blockerror_sr': [list(range(5)), 0, 0.02, [0.02] * 5, 0.06],
            'fluxlength': (50, 5),
        },
        'figures': {},
    }


class TestReportPathTables:
    """Test fmt='tex' branches in report_path table helpers."""

    def test_table_interface_tex(self):
        """_table_interface with fmt='tex' generates latex output."""
        results = [_make_result('001', (-0.9, 0.0, 1.0))]
        table, txt = _table_interface(results, fmt='tex')
        assert len(table) == 1
        assert 'tabular' in txt or 'begin' in txt or 'hline' in txt

    def test_table_probability_tex(self):
        """_table_probability with fmt='tex' generates latex output."""
        results = [_make_result('001', (-0.9, 0.0, 1.0))]
        table, txt = _table_probability(results, fmt='tex')
        assert len(table) == 1
        assert 'tabular' in txt or 'begin' in txt or 'hline' in txt

    def test_table_efficiencies_tex(self):
        """_table_efficiencies with fmt='tex' generates latex output."""
        results = [_make_result('001', (-0.9, 0.0, 1.0))]
        table, txt = _table_efficiencies(results, fmt='tex')
        assert len(table) == 1
        assert 'tabular' in txt or 'begin' in txt or 'hline' in txt

    def test_table_probability_repptis_tex(self):
        """_table_probability_repptis with fmt='tex' generates latex output."""
        results = [_make_result()]
        table, txt = _table_probability_repptis(results, fmt='tex')
        assert len(table) == 1

    def test_table_summary_tex(self):
        """_table_summary with fmt='tex' generates latex output."""
        report = {
            'numbers': {
                'pcross': '0.001', 'perr': '5%',
                'flux': '1.0', 'flux-err': '2%',
                'rate': '1e-3', 'rate-err': '3%',
            },
            'text': {'flux-unit': 'ns'},
        }
        table, txt = _table_summary(report, fmt='tex')
        assert len(table) == 3

    def test_table_interface0_tex(self):
        """_table_interface0 with fmt='tex' generates latex output."""
        results = [_make_result('000', (-float('inf'), -0.9, -0.9))]
        table, txt = _table_interface0(results, fmt='tex')
        assert len(table) == 1

    def test_generate_report_tis(self):
        """generate_report_tis assembles report from analysis dict."""
        result = _make_result('001', (-0.9, 0.0, 1.0))
        analysis = {
            'pathensemble': [result],
            'matched': {
                'figures': {
                    'matched-probability': None,
                    'total-probability': None,
                    'overall-prun': None,
                    'overall-err': None,
                },
                'out': {'prob': 0.001, 'relerror': 0.1},
            },
        }
        report = generate_report_tis(analysis)
        assert 'figures' in report
        assert 'tables' in report

    def test_generate_report_retis_permeability(self):
        """generate_report_retis with permeability=True in result0 out dict."""
        result0 = _make_result('000', (-float('inf'), -0.9, -0.9))
        result0['out']['permeability'] = True
        result = _make_result('001', (-0.9, 0.0, 1.0))
        analysis = {
            'pathensemble0': result0,
            'pathensemble': [result],
            'matched': {
                'figures': {
                    'matched-probability': None,
                    'total-probability': None,
                    'overall-rrun': None,
                    'overall-err': None,
                },
                'out': {'prob': 0.001, 'relerror': 0.1},
            },
            'flux': {'value': 0.5, 'error': 0.02, 'unit': 'ns^-1'},
            'rate': {'value': 0.0005, 'error': 0.01},
            'xi': {
                'value': 1.0, 'error': 0.1,
                'fig': {'xirun': None, 'xierror': None},
            },
            'tau': {
                'value': 2.0, 'error': 0.2,
                'fig': {'taurun': None, 'tauerror': None, 'tauref': None},
            },
            'perm': {'value': 1e-6, 'error': 1e-7},
        }
        report = generate_report_retis(analysis)
        assert 'numbers' in report
        assert report['numbers']['permeability'] is True

    def test_generate_report_retis0_else_branch(self):
        """generate_report_retis0 else branch for pathensemble_repptis."""
        result = _make_result('000', (-float('inf'), -0.9, -0.9))
        analysis = {
            'pathensemble_repptis': result,
        }
        report = generate_report_retis0(analysis)
        assert 'ensemble' in report

    def test_generate_report_retis0_non_txt_output(self):
        """generate_report_retis0 warns and falls back to txt."""
        result = _make_result('000', (-float('inf'), -0.9, -0.9))
        analysis = {'pathensemble': result}
        report = generate_report_retis0(analysis, output='rst')
        assert 'ensemble' in report


class TestReportMdTex:
    """Test fmt='tex' and fmt='txt' branches in report_md table helpers."""

    def _make_flux_results(self, n=2):
        """Return minimal flux results dict."""
        return {
            'interfaces': [0.1 * i for i in range(n)],
            'runflux': [np.array([i * 0.01]) for i in range(n)],
            'errflux': [
                [0, 0, 0.001, [0.001], 0.01, 0, 0.1] for _ in range(n)
            ],
            'ncross': [10 * (i + 1) for i in range(n)],
            'neffcross': [5 * (i + 1) for i in range(n)],
            'neffc/nc': [0.5] * n,
            'cross_time': [2.0] * n,
            'pMD': [0.3] * n,
            'teffMD': [1.5] * n,
            'corrMD': [0.8] * n,
            '1-p': [0.7] * n,
            'times': {'A': 100, 'B': 200, 'OA': 300, 'OB': 400},
            'totalcycle': 1000,
        }

    def test_table_md_efficiency_tex(self):
        """_table_md_efficiency with fmt='tex' generates latex output."""
        results = self._make_flux_results()
        table, txt = _table_md_efficiency(results, fmt='tex')
        assert len(table) == 2
        assert 'tabular' in txt or 'begin' in txt or 'hline' in txt

    def test_table_md_efficiency_txt(self):
        """_table_md_efficiency with fmt='txt' generates rst output."""
        results = self._make_flux_results()
        table, txt = _table_md_efficiency(results, fmt='txt')
        assert len(table) == 2

    def test_table_md_efficiency_rst(self):
        """_table_md_efficiency with default fmt='rst' generates rst output."""
        results = self._make_flux_results()
        table, txt = _table_md_efficiency(results, fmt='rst')
        assert len(table) == 2

    def test_table_md_flux_cycles_tex(self):
        """_table_md_flux_cycles with fmt='tex' generates latex output."""
        results = self._make_flux_results()
        table, txt = _table_md_flux_cycles(results, fmt='tex')
        assert len(table) == 5
        assert 'tabular' in txt or 'begin' in txt or 'hline' in txt

    def test_table_md_flux_tex(self):
        """_table_md_flux with fmt='tex' generates latex output."""
        results = self._make_flux_results()
        table, txt = _table_md_flux(results, fmt='tex')
        assert len(table) == 2
        assert 'tabular' in txt or 'begin' in txt or 'hline' in txt


class TestAnalysisio:
    """Test uncovered branches in analysisio."""

    def test_recursive_block_analysis_short_list(self):
        """recursive_block_analysis returns NaN for short lists."""
        result = recursive_block_analysis([1.0, 2.0], minblocks=5)
        assert len(result) == 3
        assert np.isnan(result[0][0])
        assert np.isnan(result[1])
        assert np.isnan(result[2])

    def test_get_global_probz_zero_denominator(self):
        """get_global_probz returns NaN when denominator is zero."""
        pmps = np.array([1.0, 0.0])
        pmms = np.array([0.3, 0.0])
        ppps = np.array([0.2, 0.2])
        ppms = np.array([0.5, 0.5])
        result = get_global_probz(pmps, pmms, ppps, ppms)
        assert result[0] is np.nan or np.isnan(result[0])

    def test_get_path_simulation_files_make_tis(self, tmp_path):
        """get_path_simulation_files handles make-tis-files task."""
        settings = {
            'simulation': {
                'task': 'make-tis-files',
                'interfaces': [-1.0, 0.0, 1.0],
                'exe-path': str(tmp_path),
            },
            'tis': {'detect': 0.0},
            'output': {},
        }
        all_settings, all_files = get_path_simulation_files(settings)
        assert all_settings[0] is None
        assert all_files[0] is None
