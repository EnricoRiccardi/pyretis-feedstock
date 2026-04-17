# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the common methods in pyretis.pyvisa.statistical_methods."""
import unittest
import pytest
from unittest import mock
import numpy as np
import pandas
from pyretis.pyvisa import HAS_PYVISA_REQ
if HAS_PYVISA_REQ:
    from pyretis.pyvisa import statistical_methods

Dataframe = pandas.DataFrame(np.random.rand(40, 2))
TRUE_FALSE = pandas.DataFrame(np.random.randint(0, 2, 40))
cluster_data = np.column_stack([Dataframe[0], Dataframe[1]])
settings = {'op1': 'op1', 'op2': 'op2', 'fol': '000'}
colormap = 'viridis'


@unittest.skipIf(HAS_PYVISA_REQ is False, "PyVisA reqs not installed")
class TestMethods:
    """Testing class of pyretis.pyvisa.statistical_methods."""

    @mock.patch(f"{__name__}.statistical_methods.plt")
    def test_correlation_matrix(self, mock_plt):
        """Test if correlation matrix is generated."""
        statistical_methods.correlation_matrix(Dataframe)
        # Assert plt.show has been called
        mock_plt.show.assert_called_once_with()
        # Assert plt.figure got called
        assert mock_plt.figure.called

    @mock.patch(f"{__name__}.statistical_methods.plt")
    def test_pyvisa_pca(self, mock_plt):
        """Test if PCA is performed."""
        statistical_methods.pyvisa_pca(2, settings, Dataframe, colormap)
        # Assert plt.show has been called
        mock_plt.show.assert_called_once_with()
        # Assert plt.figure got called
        assert mock_plt.figure.called

    @mock.patch(f"{__name__}.statistical_methods.plt")
    def test_k_means(self, mock_plt):
        """Test if k-means cluster plot is generated."""
        statistical_methods.k_means(2, cluster_data, settings, colormap)
        # Assert plt.show has been called
        mock_plt.show.assert_called_once_with()
        # Assert plt.figure got called
        assert mock_plt.figure.called

    @mock.patch(f"{__name__}.statistical_methods.plt")
    @mock.patch(f"{__name__}.statistical_methods.MiniBatchKMeans")
    def test_k_means_uses_keyword_api(self, mock_kmeans, mock_plt):
        """Test that K-means uses the explicit scikit-learn keyword API."""
        mock_kmeans.return_value.fit_predict.return_value = \
            np.zeros(len(cluster_data), dtype=int)

        statistical_methods.k_means(2, cluster_data, settings, colormap)

        mock_kmeans.assert_called_once_with(n_clusters=2)
        mock_plt.show.assert_called_once_with()

    @mock.patch(f"{__name__}.statistical_methods.plt")
    def test_hierarchical(self, mock_plt):
        """Test if hierarchical cluster plot is generated."""
        statistical_methods.hierarchical(2, cluster_data, settings, colormap)
        # Assert plt.show has been called
        mock_plt.show.assert_called_once_with()
        # Assert plt.figure got called
        assert mock_plt.figure.called

    @mock.patch(f"{__name__}.statistical_methods.plt")
    def test_gaussian_mixture(self, mock_plt):
        """Test if Gaussian mixture cluster plot is generated."""
        statistical_methods.gaussian_mixture(2,
                                             cluster_data,
                                             settings,
                                             colormap)
        # Assert plt.show has been called
        mock_plt.show.assert_called_once_with()
        # Assert plt.figure got called
        assert mock_plt.figure.called

    @mock.patch(f"{__name__}.statistical_methods.plt")
    @mock.patch(f"{__name__}.statistical_methods.GaussianMixture")
    def test_gaussian_mixture_uses_keyword_api(self, mock_gaussian, mock_plt):
        """Test that GaussianMixture uses the current component keyword."""
        mock_gaussian.return_value.fit_predict.return_value = \
            np.zeros(len(cluster_data), dtype=int)

        statistical_methods.gaussian_mixture(2, cluster_data, settings,
                                             colormap)

        mock_gaussian.assert_called_once_with(n_components=2)
        mock_plt.show.assert_called_once_with()

    @mock.patch(f"{__name__}.statistical_methods.plt")
    def test_spectral(self, mock_plt):
        """Test if spectral cluster plot is generated."""
        statistical_methods.spectral(2, cluster_data, settings, colormap)
        # Assert plt.show has been called
        mock_plt.show.assert_called_once_with()
        # Assert plt.figure got called
        assert mock_plt.figure.called

    @mock.patch(f"{__name__}.statistical_methods.plt")
    def test_random_forest(self, mock_plt):
        """Test if random forest plot is generated."""
        statistical_methods.random_forest(Dataframe, TRUE_FALSE, 3)
        # Assert plt.show has been called
        mock_plt.show.assert_called_once_with()
        # Assert plt.subplots got called
        assert mock_plt.figure.called

    @mock.patch(f"{__name__}.statistical_methods.graphviz")
    def test_decision_tree(self, mock_graphviz):
        """Test if decision tree plot is generated."""
        statistical_methods.decision_tree(Dataframe, TRUE_FALSE, 3)
        # Assert graphviz.Source has been called
        assert mock_graphviz.Source.called


if __name__ == '__main__':
    unittest.main()
