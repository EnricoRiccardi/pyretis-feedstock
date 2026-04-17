#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Collection of statistical methods for PyVisA.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
correlation_matrix (:py:func:`.correlation_matrix`)
    Method for displaying the correlations from the simulation variables.

gaussian_mixture (:py:func:`.gaussian_mixture`)
    Generic method for performing gaussian mixture clustering
    on simulation data.

hierarchical (:py:func:`.hierarchical`)
    Generic method for performing hierarchical clustering on simulation data.

k_means (:py:func:`.k_means`)
    Generic method for performing K-means clustering on simulation data.

pyvis_pca (:py:func:`.pyvisa_pca`)
    Generic method for performing pca on simulation data.

spectral (:py:func:`.spectral`)
    Generic method for performing spectral clustering on simulation data.

"""
import graphviz
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import AgglomerativeClustering, SpectralClustering, \
    MiniBatchKMeans
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.tree import DecisionTreeClassifier, export_graphviz


def _plot_pca_results(n_pca, principal_df, loadings, pca_info, cmap):
    """Plot the results of a PCA analysis.

    Parameters
    ----------
    n_pca : integer
        Number of principal components.
    principal_df : object like pandas.DataFrame
        Principal components data.
    loadings : object like pandas.DataFrame
        PCA loadings.
    pca_info : dict
        Dict with keys 'model' (fitted PCA), 'cols' (column labels),
        'features' (feature labels).
    cmap : string
        Matplotlib colormap.

    """
    pca_model = pca_info['model']
    cols = pca_info['cols']
    features = pca_info['features']
    # Plot PC1 vs PC2
    if n_pca > 1:
        fig_1 = plt.figure()
        axis_1 = fig_1.add_subplot(111)
        axis_1.scatter(principal_df['PC1'], principal_df['PC2'],
                       edgecolors=None, alpha=0.5,
                       cmap=cmap)
        axis_1.set_xlabel('PC1')
        axis_1.set_ylabel('PC2')

    # Plot the loadings
    fig_2 = plt.figure()
    axis_2 = fig_2.add_subplot(111)
    fig_2.colorbar(axis_2.matshow(loadings))
    axis_2.set_xticks(np.arange(len(cols)))
    axis_2.set_yticks(np.arange(len(features)))
    axis_2.set_xticklabels(cols)
    axis_2.set_yticklabels(features)

    # Plot Cumulative explained variance
    fig_3 = plt.figure()
    axis_3 = fig_3.add_subplot(111)
    axis_3.plot(np.cumsum(pca_model.explained_variance_ratio_))
    axis_3.set_xlabel('Principal components')
    axis_3.set_ylabel('Cumulative explained variance')
    plt.show()


def pyvisa_pca(n_pca, settings, data, cmap):
    """Perform PCA on ensemble data.

    This function performs PCA on data from the selected ensembles
    and stores the data to a hdf-file for further analysis, plots the
    first two principal components against each other,
    the cumulative explained variance and the loadings.
    The hdf5 file contains the keywords:

    * `sim_data`: pandas Dataframe containing the simulation data.
    * `PC`: Dataframe containing the data from the principal components.
    * `loadings`: pandas Dataframe containing the loadings.
    * `explained`: pandas Series containing the simulation data.
    * `correlation_matrix`: Dataframe containing the correlation matrix.


    Parameters
    ----------
    n_pca : integer
        Number of clusters.
    settings : dict
        Settings from GUI.
    data : object like pandas.DataFrame
        Simulation data from chosen ensembles.
    cmap : string
        Matplotlib colormap.

    """
    basename = f'PCA_data_{settings["fol"]}_{n_pca}.hdf'
    data.to_hdf(basename, key='sim_data', mode='w')

    features = data.columns
    data = StandardScaler().fit_transform(data)

    cols = [f'PC{i}' for i in range(1, n_pca + 1)]

    pca_model = PCA(n_components=n_pca)
    principal_comps = pca_model.fit_transform(data)
    principal_df = pd.DataFrame(principal_comps, columns=cols)
    explained = pd.Series(pca_model.explained_variance_ratio_)
    explained.index += 1
    loadings = pd.DataFrame(pca_model.components_.T, columns=cols,
                            index=features)
    load_corr = pca_model.components_.T * np.sqrt(
        pca_model.explained_variance_)
    load_corr_mat = pd.DataFrame(load_corr, columns=cols,
                                 index=features)

    principal_df.to_hdf(basename, key='PC', mode='a')
    loadings.to_hdf(basename, key='loadings', mode='a')
    load_corr_mat.to_hdf(basename, key='correlation_matrix', mode='a')
    explained.to_hdf(basename, key='explained', mode='a')

    pca_info = {'model': pca_model, 'cols': cols, 'features': features}
    _plot_pca_results(n_pca, principal_df, loadings, pca_info, cmap)


def k_means(n_clusters, data, settings, cmap):
    """Perform K-means clustering.

    This function performs K-means clustering on simulation data
    with a chosen amount of clusters, and plots the results.

    Parameters
    ----------
    n_clusters : integer
        Number of clusters.
    data : object like numpy.ndarray
        Simulation data from chosen ensembles.
    settings : dict
        Settings from GUI.
    cmap : string
        Matplotlib colormap.

    """
    model = MiniBatchKMeans(n_clusters=n_clusters)
    labels = model.fit_predict(data)
    fig = plt.figure()
    axis = fig.add_subplot(111)
    axis.scatter(data[:, 0], data[:, 1],
                 c=labels,
                 s=5,
                 cmap=cmap)
    axis.set_title('Kmeans clustering')
    axis.set_xlabel(settings['op1'])
    axis.set_ylabel(settings['op2'])
    plt.show()


def hierarchical(n_clusters, data, settings, cmap):
    """Perform hierarchical clustering.

    This function performs hierarchical clustering on simulation data
    with a chosen amount of clusters, and plots the results.

    Parameters
    ----------
    n_clusters : integer
        Number of clusters.
    data : object like numpy.ndarray
        Simulation data from chosen ensembles.
    settings : dict
        Settings from GUI.
    cmap : string
        Matplotlib colormap.

    """
    model = AgglomerativeClustering(n_clusters=n_clusters,
                                    linkage='average')
    labels = model.fit_predict(data)
    fig = plt.figure()
    axis = fig.add_subplot(111)
    axis.scatter(data[:, 0], data[:, 1],
                 c=labels,
                 s=5,
                 cmap=cmap)
    axis.set_title('Hierarchical clustering')
    axis.set_xlabel(settings['op1'])
    axis.set_ylabel(settings['op2'])
    plt.show()


def spectral(n_clusters, data, settings, cmap):
    """Perform spectral clustering.

    This function performs spectral clustering on simulation data
    with a chosen amount of clusters, and plots the results.

    Parameters
    ----------
    n_clusters : integer
        Number of clusters.
    data : object like numpy.ndarray
        Simulation data from chosen ensembles.
    settings : dict
        Settings from GUI.
    cmap : string
        Matplotlib colormap.

    """
    model = SpectralClustering(
        n_clusters=n_clusters,
        affinity='nearest_neighbors',
        n_neighbors=30,
        assign_labels='kmeans')
    labels = model.fit_predict(data)
    fig = plt.figure()
    axis = fig.add_subplot(111)
    axis.scatter(data[:, 0], data[:, 1],
                 c=labels,
                 s=5,
                 cmap=cmap)
    axis.set_title('Spectral clustering')
    axis.set_xlabel(settings['op1'])
    axis.set_ylabel(settings['op2'])
    plt.show()


def gaussian_mixture(n_clusters, data, settings, cmap):
    """Perform gaussian mixture clustering.

    This function performs gaussian mixture clustering on simulation data
    with a chosen amount of clusters, and plots the results.

    Parameters
    ----------
    n_clusters : integer
        Number of clusters.
    data : object like numpy.ndarray
        Simulation data from chosen ensembles.
    settings : dict
        Settings from GUI.
    cmap : string
        Matplotlib colormap.

    """
    model = GaussianMixture(n_components=n_clusters)
    labels = model.fit_predict(data)
    fig = plt.figure()
    axis = fig.add_subplot(111)
    axis.scatter(data[:, 0], data[:, 1],
                 c=labels,
                 s=5,
                 cmap=cmap)
    axis.set_title('Gaussian mixture clustering')
    axis.set_xlabel(settings['op1'])
    axis.set_ylabel(settings['op2'])
    plt.show()


def correlation_matrix(dataframe):
    """Show correlation matrix of simulation variables.

    Function which computes the Pearson correlations between parameters
    and shows it in a correlation plot.

    Parameters
    ----------
    dataframe : object like pandas.DataFrame
        Dataframe containing simulation data.

    """
    titles = list(dataframe.columns)
    fig = plt.figure()
    axis = fig.add_subplot(111)
    corr = axis.matshow(dataframe.corr())
    fig.colorbar(corr)
    axis.set_xticks(np.arange(len(titles)))
    axis.set_yticks(np.arange(len(titles)))
    axis.set_xticklabels(titles)
    axis.set_yticklabels(titles)
    plt.show()


def random_forest(xdata, ydata, depth):
    """Create random forest classification.

    Parameters
    ----------
    xdata : object like pandas.DataFrame
        Pandas dataframe containing the op and cv data from selected
        frames in the simulation.
    ydata : object like pandas.DataFrame
        Pandas dataframe containing True/False for each frame in the
        selected ensembles. True if the frame is reactive, else False.
    depth : integer
        Depth of random forest model.

    """
    rf_model = RandomForestClassifier(n_estimators=500,
                                      max_depth=depth,
                                      random_state=16)
    ydata = ydata.values.ravel()
    rf_model.fit(xdata, ydata)
    important_features = np.argsort(rf_model.feature_importances_)
    feature_imp = pd.Series(rf_model.feature_importances_,
                            index=xdata.columns).sort_values(ascending=False)
    figure = plt.figure()
    axis = figure.add_subplot(111)
    x_ticks = [i + 0.5 for i in range(len(important_features))]
    axis.bar(x_ticks, feature_imp)
    axis.set_xticks(x_ticks)
    axis.set_xticklabels(feature_imp.keys(), rotation=90, ha='center')
    axis.set_ylabel('Gini feature importance')
    plt.show()


def decision_tree(xdata, ydata, depth):
    """Create a decision tree.

    Parameters
    ----------
    xdata : object like pandas.DataFrame
        Pandas dataframe containing the op and cv data from selected
        frames i the simulation.
    ydata : object like pandas.DataFrame
        Pandas dataframe containing True/False for each frame in the
        selected ensembles. True if the frame is reactive, else False.
    depth : integer
        Depth of the decision tree.

    """
    names = xdata.columns.astype(str).tolist()
    x_train, _, y_train, _ = train_test_split(xdata, ydata,
                                              test_size=0.2,
                                              stratify=ydata,
                                              random_state=16,
                                              shuffle=True)

    clf = DecisionTreeClassifier(criterion='entropy', max_depth=depth)
    clf.fit(x_train, y_train)

    dot_data = export_graphviz(clf, out_file=None,
                               feature_names=names,
                               class_names=['reactive', 'unreactive'],
                               filled=True, rounded=True,
                               special_characters=True)
    graph = graphviz.Source(dot_data)
    graph.render('decision_tree', view=True, format='png', cleanup=True)
