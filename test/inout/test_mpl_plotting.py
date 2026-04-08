# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Tests to raise coverage for mpl_plotting.py.

These tests exercise branches not covered by the main test suite.
"""
import os
import logging
import numpy as np
import pytest

import matplotlib
matplotlib.use('Agg')  # must be called before any other matplotlib imports
from pyretis.inout.plotting.mpl_plotting import (  # noqa: E402
    MplPlotter,
    _mpl_read_style_file,
    mpl_set_style,
    mpl_savefig,
    mpl_plot_in_chunks,
    _mpl_plot_xy_chunk,
    mpl_simple_plot,
    mpl_linecollection_gradient,
    mpl_chunks_gradient,
    mpl_line_gradient,
    mpl_error_plot,
    _mpl_shoots_histogram,
    mpl_plot_orderp,
    mpl_plot_energy,
    mpl_plot_tau,
    mpl_plot_tau_ref,
    mpl_plot_matched,
    mpl_plot_pppath,
    get_color_map,
)
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.backends.backend_agg import (  # noqa: E402
    FigureCanvasAgg as FigureCanvas
)

logging.disable(logging.CRITICAL)


class TestMplSetStyle:
    """Test mpl_set_style branches."""

    def test_named_style(self):
        """mpl_set_style with a named matplotlib style (not 'pyretis')."""
        mpl_set_style('ggplot')

    def test_read_style_file(self, tmp_path):
        """_mpl_read_style_file reads a file and applies valid rcParams."""
        style_file = tmp_path / 'test.mplstyle'
        style_file.write_text(
            '# comment line\n'
            'lines.linewidth: 2\n'
            'bogus.key: value\n',
            encoding='utf-8'
        )
        _mpl_read_style_file(str(style_file))

    def test_read_style_file_color(self, tmp_path):
        """_mpl_read_style_file prepends '#' for color keys."""
        style_file = tmp_path / 'test.mplstyle'
        style_file.write_text('axes.facecolor: ffffff\n', encoding='utf-8')
        _mpl_read_style_file(str(style_file))


class TestMplSavefig:
    """Test mpl_savefig backup branch."""

    def test_savefig_with_backup(self, tmp_path):
        """mpl_savefig creates a backup when backup=True and file exists."""
        outfile = str(tmp_path / 'fig.png')
        # Write a dummy file first so backup is triggered
        with open(outfile, 'w') as f:
            f.write('placeholder')
        fig = Figure()
        canvas = FigureCanvas(fig)
        mpl_savefig(canvas, outfile, backup=True)
        # The original file was backed up and a new one written
        assert os.path.isfile(outfile)


class TestMplPlotInChunks:
    """Test mpl_plot_in_chunks chunk processing."""

    def test_large_series_chunked(self):
        """mpl_plot_in_chunks processes data larger than chunksize."""
        fig = Figure()
        canvas = FigureCanvas(fig)
        axs = fig.add_subplot(111)
        n = 25000
        series = {'type': 'xy', 'x': list(range(n)), 'y': list(range(n))}
        handle = mpl_plot_in_chunks(axs, series, chunksize=10000)
        assert handle is not None

    def test_chunk_with_color(self):
        """_mpl_plot_xy_chunk passes color when given."""
        fig = Figure()
        canvas = FigureCanvas(fig)
        axs = fig.add_subplot(111)
        series = {'x': list(range(10)), 'y': list(range(10))}
        handle = _mpl_plot_xy_chunk(axs, series, color='red')
        assert handle is not None

    def test_chunk_with_high(self):
        """_mpl_plot_xy_chunk clips at high when specified."""
        fig = Figure()
        canvas = FigureCanvas(fig)
        axs = fig.add_subplot(111)
        series = {'x': list(range(20)), 'y': list(range(20))}
        handle = _mpl_plot_xy_chunk(axs, series, low=0, high=5)
        assert handle is not None


class TestMplSimplePlot:
    """Test mpl_simple_plot with bar and line type series."""

    def test_bar_series(self):
        """mpl_simple_plot with bar type series."""
        series = [{'type': 'bar', 'x': [1, 2, 3], 'height': [0.3, 0.5, 0.2],
                   'width': 0.5, 'label': 'test bar', 'align': 'edge'}]
        canvas = mpl_simple_plot(series, fig_settings={})
        assert canvas is not None

    def test_hline_series(self):
        """mpl_simple_plot with hline type series."""
        series = [{'type': 'hline', 'y': 0.5, 'label': 'threshold'}]
        canvas = mpl_simple_plot(series, fig_settings={'ylabel': 'y'})
        assert canvas is not None


class TestMplGradient:
    """Test gradient plotting functions."""

    def test_linecollection_gradient(self):
        """mpl_linecollection_gradient plots a gradient line collection."""
        fig = Figure()
        canvas = FigureCanvas(fig)
        axs = fig.add_subplot(111)
        series = {'x': list(range(10)),
                  'y': [float(i) ** 0.5 for i in range(10)]}
        handle = mpl_linecollection_gradient(axs, series)
        assert handle is not None

    def test_chunks_gradient(self):
        """mpl_chunks_gradient plots large data in gradient chunks."""
        fig = Figure()
        canvas = FigureCanvas(fig)
        axs = fig.add_subplot(111)
        n = 25000
        series = {'x': list(range(n)), 'y': [float(i) for i in range(n)]}
        handle = mpl_chunks_gradient(axs, series, chunksize=10000)
        assert handle is not None

    def test_line_gradient_large(self):
        """mpl_line_gradient uses chunk path for n >= 10**6."""
        n = 10 ** 6 + 1
        series = [{'x': list(range(n)), 'y': list(range(n))}]
        canvas = mpl_line_gradient(series, fig_settings={})
        assert canvas is not None

    def test_line_gradient_small(self):
        """mpl_line_gradient uses linecollection path for small n."""
        series = [{'x': list(range(10)), 'y': list(range(10)),
                   'label': 'test'}]
        canvas = mpl_line_gradient(series,
                                   fig_settings={'xlabel': 'x', 'ylabel': 'y',
                                                 'title': 'T'})
        assert canvas is not None


class TestMplErrorPlot:
    """Test mpl_error_plot branches."""

    def test_error_plot_with_legend(self):
        """mpl_error_plot with 4-element tuple adds legend."""
        x = np.linspace(0, 1, 20)
        y = np.sin(x)
        err = 0.1 * np.ones_like(x)
        series = [(x, y, err, 'sin(x)')]
        canvas = mpl_error_plot(series,
                                fig_settings={'xlabel': 'x', 'ylabel': 'y',
                                              'title': 'T'})
        assert canvas is not None

    def test_error_plot_no_legend(self):
        """mpl_error_plot without legend (3-element tuple)."""
        x = np.linspace(0, 1, 10)
        y = np.cos(x)
        err = 0.05 * np.ones_like(x)
        series = [(x, y, err)]
        canvas = mpl_error_plot(series, fig_settings={})
        assert canvas is not None


class TestMplShootsHistogram:
    """Test _mpl_shoots_histogram."""

    def test_shoots_histogram_missing_key(self):
        """_mpl_shoots_histogram handles missing histogram keys gracefully."""
        # Only 'ACC' key; if empty, loop continues without breaking
        histograms = {}
        canvas = _mpl_shoots_histogram(histograms, '001')
        assert canvas is not None

    def test_shoots_histogram_acc(self):
        """_mpl_shoots_histogram plots ACC histogram."""
        histograms = {'ACC': (np.array([0.1, 0.2]), np.array([0, 1, 2]),
                              np.array([0.5, 1.5]))}
        canvas = _mpl_shoots_histogram(histograms, '001')
        assert canvas is not None


class TestMplPlotOrderp:
    """Test mpl_plot_orderp velocity branch."""

    def _make_results(self):
        """Return minimal order parameter analysis results."""
        n = 20
        block = [list(range(n)), 0, 0, [0.01] * n, 0.05, 0, 0.1]
        dist = (np.ones(n), np.linspace(0, 1, n), (0.5, 0.1))
        return [{
            'running': np.ones(n),
            'blockerror': block,
            'distribution': dist,
        }]

    def test_orderp_with_velocity_column(self):
        """mpl_plot_orderp generates velocity plot when col >= 3."""
        n = 20
        results = self._make_results()
        orderdata = np.column_stack([
            np.linspace(0, 1, n),
            np.random.rand(n),
            np.random.rand(n),
        ])
        canvas = mpl_plot_orderp(results, orderdata)
        assert canvas is not None

    def test_orderp_with_msd(self):
        """mpl_plot_orderp generates MSD plot when 'msd' key present."""
        n = 20
        results = self._make_results()
        results[0]['msd'] = np.column_stack([
            np.arange(n, dtype=float),
            0.01 * np.ones(n),
        ])
        orderdata = np.column_stack([
            np.linspace(0, 1, n),
            np.random.rand(n),
        ])
        canvas = mpl_plot_orderp(results, orderdata)
        assert canvas is not None


class TestMplPlotEnergy:
    """Test mpl_plot_energy Boltzmann branch."""

    def test_energy_with_boltzmann(self):
        """mpl_plot_energy generates Boltzmann distribution plot."""
        n = 20
        time = np.linspace(0, 1, n)
        energies = {
            'time': time,
            'vpot': np.random.rand(n),
            'ekin': np.random.rand(n),
        }
        block = [list(range(n)), 0, 0, [0.01] * n, 0.05, 0, 0.1]
        dist = (np.ones(n), np.linspace(0, 1, n), (0.5, 0.1))
        results = {
            'vpot': {
                'running': np.ones(n),
                'blockerror': block,
                'distribution': dist,
                'boltzmann-dist': (np.ones(n), np.linspace(0, 1, n)),
            },
            'ekin': {
                'running': np.ones(n),
                'blockerror': block,
                'distribution': dist,
            },
            'temp': {
                'running': np.ones(n),
                'blockerror': block,
                'distribution': dist,
            },
        }
        energies['temp'] = np.random.rand(n)
        canvas = mpl_plot_energy(results, energies)
        assert canvas is not None


class TestMplPlotTau:
    """Test mpl_plot_tau zero-tau branch and mpl_plot_tau_ref."""

    def test_tau_all_zeros(self):
        """mpl_plot_tau uses 'no' when tau is all zeros."""
        n = 20
        tau = np.zeros(n)
        tauerror = [list(range(n)), 0, 0, [0.01] * n, 0.05, 0, 0.1]
        canvas_run, canvas_err = mpl_plot_tau(tau, tauerror)
        assert len(canvas_run) == 1
        assert len(canvas_err) == 1

    def test_tau_nonzero(self):
        """mpl_plot_tau uses 'in' when tau has non-zero values."""
        n = 20
        tau = np.ones(n) * 0.5
        tauerror = [list(range(n)), 0, 0, [0.01] * n, 0.05, 0, 0.1]
        canvas_run, canvas_err = mpl_plot_tau(tau, tauerror)
        assert len(canvas_run) == 1

    def test_tau_ref_with_two_bins(self):
        """mpl_plot_tau_ref adds vertical lines when len(bins)==2."""
        tau_ref = np.array([0.1, 0.2, 0.3])
        tau_ref_bins = [(0.0, 0.33), (0.33, 0.66), (0.66, 1.0)]
        bins = (0.1, 0.9)
        result = mpl_plot_tau_ref(tau_ref, tau_ref_bins, bins)
        assert 'canvas' in result
        assert result['name'] == 'tauref'

    def test_tau_ref_no_bins(self):
        """mpl_plot_tau_ref without vlines when len(bins) != 2."""
        tau_ref = np.array([0.1, 0.2])
        tau_ref_bins = [(0.0, 0.5), (0.5, 1.0)]
        result = mpl_plot_tau_ref(tau_ref, tau_ref_bins, bins=())
        assert 'canvas' in result


class TestMplPlotMatched:
    """Test mpl_plot_matched with reptis=True and many ensembles."""

    def _make_matched(self, with_rrun=False):
        """Return minimal matched dict."""
        n = 20
        overall_prob = np.column_stack([
            np.linspace(0, 1, n),
            np.exp(-np.linspace(0, 5, n)),
        ])
        matched = {
            'overall-prob': overall_prob,
            'overall-cycle': list(range(n)),
            'overall-error': [list(range(n)), 0, 0, [0.01] * n, 0.05, 0, 0.1],
            'matched-prob': [overall_prob],
            'pcross': np.exp(-np.linspace(0, 5, 5)),
        }
        if with_rrun:
            matched['overall-rrun'] = np.ones(n)
        else:
            matched['overall-prun'] = np.ones(n)
        return matched

    def test_matched_reptis_true(self):
        """mpl_plot_matched with reptis=True constructs pcross array."""
        n = 20
        detect = np.linspace(0, 1, 5)
        matched = self._make_matched()
        canvas = mpl_plot_matched(['001'], detect, matched, reptis=True)
        assert canvas is not None

    def test_matched_with_rrun(self):
        """mpl_plot_matched uses overall-rrun when present."""
        detect = [0.3, 0.6]
        matched = self._make_matched(with_rrun=True)
        canvas = mpl_plot_matched(['001'], detect, matched, reptis=False)
        assert canvas is not None

    def test_matched_many_ensembles(self):
        """mpl_plot_matched overrides color cycle for many ensembles."""
        n = 20
        overall_prob = np.column_stack([
            np.linspace(0, 1, n),
            np.exp(-np.linspace(0, 5, n)),
        ])
        detect = list(np.linspace(0, 1, 15))
        matched = {
            'overall-prob': overall_prob,
            'overall-prun': np.ones(n),
            'overall-cycle': list(range(n)),
            'overall-error': [list(range(n)), 0, 0, [0.01] * n, 0.05, 0, 0.1],
            'matched-prob': [overall_prob] * 15,
        }
        path_ensembles = [f'{i:03d}' for i in range(15)]
        canvas = mpl_plot_matched(path_ensembles, detect, matched)
        assert canvas is not None

    def test_matched_no_rrun_no_prun(self):
        """mpl_plot_matched returns early when neither rrun nor prun."""
        n = 20
        overall_prob = np.column_stack([
            np.linspace(0, 1, n),
            np.exp(-np.linspace(0, 5, n)),
        ])
        matched = {'overall-prob': overall_prob}
        detect = [0.5]
        canvas = mpl_plot_matched([], detect, matched)
        assert canvas is not None


class TestMplPlotterOutputPpGlobalCross:
    """Test MplPlotter.output_pp_global_cross."""

    def test_output_pp_global_cross(self, tmp_path):
        """output_pp_global_cross generates a figure."""
        plotter = MplPlotter('png', style='pyretis', out_dir=str(tmp_path))
        pp = np.exp(-np.linspace(0, 3, 10))
        files = plotter.output_pp_global_cross(pp)
        assert 'pp_Pc' in files


class TestMplPlotPppath:
    """Test mpl_plot_pppath with pp_Pc results."""

    def test_pppath_with_pp_Pc(self):
        """mpl_plot_pppath generates pp_Pc canvas when key is present."""
        from pyretis.core.pathensemble import PathEnsemble
        ens = PathEnsemble(1, [-0.9, -0.9, 1.0])
        n = 10
        results = {
            'pp_Pc': np.exp(-np.linspace(0, 3, n)),
            'pathlength': [([0.1] * n, list(range(n)), (5.0, 1.0))],
            'shoots': [{}],
        }
        canvas = mpl_plot_pppath(results, ens)
        assert canvas is not None


class TestGetColorMap:
    """Test get_color_map for more than 20 colors."""

    def test_get_color_map_large(self):
        """get_color_map with ncolors > 20 uses tab20."""
        colors = get_color_map(25)
        assert len(colors) == 25
