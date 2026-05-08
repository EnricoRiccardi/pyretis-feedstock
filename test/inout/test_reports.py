# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the FileIO class."""
import logging
import numpy as np
import os
import pytest
from pyretis.analysis import match_all_histograms
from pyretis.inout.plotting import create_plotter
from pyretis.inout.report.report import generate_report, get_template
from .help import turn_on_logging

logging.disable(logging.CRITICAL)


HERE = os.path.abspath(os.path.dirname(__file__))


class TestReporting:
    """Test the generate_report."""
    def test_generate_report(self, caplog):
        """Test analysisio function"""
        with pytest.raises(ValueError) as ext:
            generate_report(report_type='fake', analysis_results={},
                            output='ext')
        assert 'nkown report type fake' in str(ext.value)
        with turn_on_logging():
            with caplog.at_level(logging.WARNING,
                                 logger='pyretis.inout.report.report'):
                with pytest.raises(ValueError) as ext:
                    generate_report(report_type='retis0',
                                    analysis_results={},
                                    output='ext')
                assert 'Could not locate template' in str(ext.value)

    def test_create_plotter(self):
        """Test create_plotter."""
        with pytest.raises(ValueError) as ext:
            create_plotter({'plotter': 'Fake'})
        assert 'Unknown plotter: Fake' in str(ext.value)

        with pytest.raises(ValueError) as ext:
            create_plotter({'plotter': 'mpl', 'output': 'kik'})
        assert 'Output format "kik" is not support' in str(ext.value)

    def test_get_template(self):
        """Test get_template."""
        # Any file would do
        ifile = os.path.join(HERE, 'initial-gro.rst')
        nav, path = get_template(output='', report_type='', template=ifile)

        assert (nav, path) == ('initial-gro.rst', HERE)

    def test_match_all_histograms(self):
        """Test match histograms for umbrella sampling."""
        histograms = [[np.array([1, 5, 1, 0, 1, 1]),
                       np.array([+0.76, -0.74, -0.72, -0.71, +0.68, -0.66]),
                       np.array([-0.26, -0.34, +0.42, -0.52, -0.18, -0.16])],
                      [np.array([1, 3, 1, 0, 3, 1]),
                       np.array([-0.76, -0.74, +0.72, -0.71, -0.68, -0.66]),
                       np.array([-0.26, +0.34, -0.42, +0.32, +0.18, -0.16])],
                      [np.array([1, 2, 1, 0, 7, 1]),
                       np.array([-0.76, -0.74, +0.72, -0.71, -0.68, +0.66]),
                       np.array([-0.26, +0.34, -0.42, -0.22, -0.18, -0.16])]]
        windows = [[-1.0, -0.4], [-0.5, -0.2], [-0.3, 0.0]]

        histograms_s, mix, hist_avg = match_all_histograms(histograms, windows)
        assert all(histograms_s[0] == [1, 5, 1, 0, 1, 1])
        assert all(histograms_s[1] == [1., 3., 1., 0., 3., 1.])
        assert all(histograms_s[2] == [1., 2., 1., 0., 7., 1.])
        assert np.all(np.array(mix) == [1.0, 1.0, 1.0])
        assert np.all(hist_avg == np.array([1., 3., 0., 0., 7., 1.]))
