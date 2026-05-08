# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This package defines analysis tools for the PyRETIS program.

The analysis tools are intended to be used for analysis of the
simulation output from the PyRETIS program. The typical use of this
package is in post-processing of the results from a simulation (or
several simulations).

Package structure
-----------------

Modules
~~~~~~~

__init__.py
    This file, imports from the other modules. The method to analyse
    results from MD flux simulations is defined here since it will
    make use of analysis tools from `energy_analysis.py` and
    `order_analysis.py`.

analysis.py (:py:mod:`pyretis.analysis.analysis`)
    General methods for numerical analysis.

energy_analysis.py (:py:mod:`pyretis.analysis.energy_analysis`)
    Defines methods useful for analysing the energy output.

histogram.py (:py:mod:`pyretis.analysis.histogram`)
    Defines methods useful for generating histograms.

order_analysis.py (:py:mod:`pyretis.analysis.order_analysis`)
    Defines methods useful for analysis of order parameters.

path_analysis.py (:py:mod:`pyretis.analysis.path_analysis`)
    Defines methods for analysis of path ensembles.

Important methods defined in this package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

analyse_energies (:py:func:`.analyse_energies`)
    Analyse energy data from a simulation. It will calculate a running
    average, a distribution and do a block error analysis.

analyse_flux (:py:func:`.analyse_flux`)
    Analyse flux data from a MD flux simulation. It will calculate a
    running average, a distribution and do a block error analysis.

analyse_orderp (:py:func:`.analyse_orderp`)
    Analyse order parameter data. It will calculate a running average,
    a distribution and do a block error analysis. In addition, it will
    analyse the mean square displacement (if requested).

analyse_path_ensemble (:py:func:`.analyse_path_ensemble`)
    Analyse the results from a single path ensemble. It will calculate
    a running average of the probabilities, a crossing probability,
    perform a block error analysis, analyse lengths of paths,
    type of Monte Carlo moves and calculate an efficiency.

match_probabilities (:py:func:`.match_probabilities`)
    Method to match probabilities from several path simulations.
    Useful for obtaining the overall crossing probability.

histogram (:py:func:`.histogram`)
    Generates histogram, basically a wrapper around
    numpy's ``histogram``.

match_all_histograms (:py:func:`.match_all_histograms`)
    Method to match histograms from umbrella simulations.

retis_flux (:py:func:`.retis_flux`)
    Method for calculating the initial flux for RETIS simulations.

retis_rate (:py:func:`.retis_rate`)
    Method for calculating the rate constant for RETIS simulations.
"""
from .analysis import running_average, block_error, block_error_corr
from .energy_analysis import analyse_energies
from .flux_analysis import analyse_flux
from .histogram import histogram, match_all_histograms
from .order_analysis import analyse_orderp
from .path_analysis import (analyse_path_ensemble, analyse_repptis_ensemble,
                            match_probabilities, retis_flux, retis_rate,
                            perm_calculations)


def analyse_md_flux(crossdata, energydata, orderdata, settings):
    """Analyse the output from a MD-flux simulation.

    The obtained results will be returned as a convenient structure for
    plotting or reporting.

    Parameters
    ----------
    crossdata : numpy.array
        This is the data containing information about crossings.
    energydata : numpy.array
        This is the raw data for the energies.
    orderdata : numpy.array
        This is the raw data for the order parameter.
    settings : dict
        The settings for the analysis (e.g block length for error
        analysis) and some settings from the simulation (interfaces,
        time step etc.).

    Returns
    -------
    results : dict
        This dict contains the results from the different analysis as a
        dictionary. This dict can be used further for plotting or for
        generating reports.

    """
    results = {'flux': analyse_flux(crossdata, settings),
               'energy': analyse_energies(energydata, settings),
               'order': analyse_orderp(orderdata, settings)}
    return results
