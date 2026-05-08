# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Module defining functions for analysis of order parameters.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

analyse_orderp (:py:func:`.analyse_orderp`)
    Run a simple order parameter analysis.
"""
from pyretis.analysis.analysis import analyse_data, mean_square_displacement


__all__ = ['analyse_orderp']


def analyse_orderp(orderdata, settings):
    """Run the analysis on several order parameters.

    The results are collected into a structure which is convenient for
    plotting.

    Parameters
    ----------
    orderdata : numpy.arrays
        The data read from the order parameter file.
    settings : dict
        This dictionary contains settings for the analysis.

    Returns
    -------
    results : numpy.array
        For each order parameter `key`, `results[key]` contains the
        result of the analysis.

    See Also
    --------
    `analyse_data` in `pyretis.analysis.analysis.py`

    Note
    ----
    We here (and in the plotting routines) make certain assumptions
    about the structure, i.e. the positions are assumed to have a
    specific meaning: column zero is the time, column one the order
    parameter and so on.

    """
    results = []
    _, col = orderdata.shape
    for i in range(col):
        if i == 0:  # first column is just the time, skip it
            pass
        else:
            result = analyse_data(orderdata[:, i], settings)
            if i == 1:  # assume that we want the MSD analysis here:
                ndt = settings['analysis']['maxordermsd']
                result['msd'] = mean_square_displacement(
                    orderdata[:, i],
                    ndt=ndt
                )
            results.append(result)
    return results
