# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Methods for analysing energy data from simulations.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

analyse_energies (py:func:`.analyse_energies`)
    Run the analysis for energies (kinetic, potential etc.).
"""
import numpy as np
from scipy.stats import gamma
from pyretis.analysis.analysis import analyse_data
from pyretis.core.units import CONSTANTS

np.set_printoptions(legacy='1.21')

__all__ = ['analyse_energies']


def analyse_energies(energies, settings):
    """Run the energy analysis on several energy types.

    The function will run the energy analysis on several energy types
    and collect the energies into a structure which is convenient for
    plotting the results.

    Parameters
    ----------
    energies : dict
        This dict contains the energies to analyse.
    settings : dict
        This dictionary contains settings for the analysis.

    Returns
    -------
    results : dict
        For each energy key `results[key]` contains the result from the
        energy analysis.

    See Also
    --------
    :py:func:`.analyse_data` in :py:mod:`pyretis.analysis.analysis.py`.

    """
    results = {}
    for key in energies:
        results[key] = analyse_data(energies[key], settings)
    # For the energy analysis it is also useful to add some
    # theoretical distributions:
    # Guess n particles, but overwrite it if explicitly stated.
    npart = len(settings['particles'].get('mass', []))
    npart = settings['particles'].get('npart', npart)
    alp = 0.5 * npart * settings['system']['dimensions']
    beta = settings['system'].get('beta')
    if beta is None:
        beta = (1.0 / (settings['system']['temperature'] *
                CONSTANTS['kB'][settings['system'].get('units', 'lj')]))
    scale = {'ekin': 1.0 / beta,
             'temp': settings['system']['temperature'] / alp}
    for key, value in scale.items():
        if key in results:
            dist = results[key]['distribution']
            pos = np.linspace(min(0.0, dist[1].min()), dist[1].max(), 1000)
            tdist = gamma.pdf(pos, alp, loc=0, scale=value)
            results[key]['boltzmann-dist'] = [tdist, pos, (0.0, value)]
    return results
