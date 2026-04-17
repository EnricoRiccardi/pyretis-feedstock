# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Methods for analysis of crossings for flux data.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

analyse_flux (:py:func:`.analyse_flux`)
    Run analysis for simulation flux data. This will calculate
    the initial flux for a simulation.
"""
import numpy as np
from numpy import divide  # pylint: disable=no-name-in-module
from pyretis.analysis.analysis import running_average, block_error_corr
from pyretis.core.path import check_crossing

np.set_printoptions(legacy='1.21')

__all__ = ['analyse_flux', 'find_crossings']


def analyse_flux(fluxdata, settings):
    """Run the analysis on the given flux data.

    This will run the flux analysis and collect the results into a
    structure which is convenient for plotting and reporting the
    results.

    Parameters
    ----------
    fluxdata : list of tuples of integers
        This array contains the data obtained from a MD simulation
        for the fluxes.
    settings : dict
        This dict contains the settings for the analysis. Note that
        this dictionary also needs some settings from the simulation,
        in particular the number of cycles, the interfaces and
        information about the time step.

    Returns
    -------
    results : dict
        This dict contains the results from the flux analysis.
        The keys are defined in the `results` variable.

    """
    end_step = (settings['simulation'].get('endcycle') or
                settings['simulation'].get('steps'))
    time_subcycles = settings['engine'].get('subcycles', 1)
    time_step = settings['engine']['timestep']*time_subcycles
    interfaces = list(settings['simulation']['interfaces'])
    results = {'eff_cross': [],  # effective crossings times
               'ncross': None,  # number of crossings
               'neffcross': [],  # number of effective crossings
               'times': {},  # time spent in the different states
               'flux': [],  # store raw flux data
               'runflux': [],  # running average of flux
               'errflux': [],  # block error analysis
               'interfaces': interfaces,
               'totalcycle': end_step,  # store total number of cycles
               'cross_time': [],  # steps per crossing
               'neffc/nc': [],  # Effective crossings per crossing
               'pMD': [],  #
               '1-p': [],  #
               'teffMD': [],  #
               'corrMD': []}  #
    if not fluxdata:
        nint = len(results['interfaces'])
        results['eff_cross'] = [[] for _ in range(nint)]
        results['ncross'] = [0] * nint
        results['neffcross'] = [0] * nint
        results['times'] = {'A': end_step, 'B': 0, 'OA': end_step, 'OB': 0}
        results['cross_time'] = [np.inf] * nint
        results['neffc/nc'] = [np.nan] * nint
        for _ in results['interfaces']:
            results['flux'].append(np.column_stack(
                (np.array([0]), np.array([0]), np.array([0.0]))
            ))
            results['runflux'].append(np.array([0.0]))
            results['errflux'].append([[], np.nan, np.nan, [], np.nan,
                                       np.nan, np.nan])
            results['pMD'].append(0.0)
            results['1-p'].append(np.inf)
            results['teffMD'].append(np.nan)
            results['corrMD'].append(np.nan)
        return results
    ret = _effective_crossings(fluxdata, len(results['interfaces']), end_step)
    results['eff_cross'] = ret[0]
    results['ncross'] = ret[1]
    results['neffcross'] = ret[2]
    results['times'] = ret[3]
    analysis = settings['analysis']
    for i in range(len(results['interfaces'])):
        time, ncross, flux = _calculate_flux(results['eff_cross'][i],
                                             results['times']['OA'],
                                             analysis['skipcross'],
                                             time_step)
        results['flux'].append(np.column_stack((time, ncross, flux)))
        # now it's also a good time to obtain running averages etc.:
        results['runflux'].append(running_average(flux))
        block_error = block_error_corr(flux,
                                       maxblock=analysis['maxblock'],
                                       blockskip=analysis['blockskip'])
        results['errflux'].append(block_error)

    # do some additional statistics:
    results['cross_time'] = [divide(float(end_step), float(neff))
                             for neff in results['neffcross']]

    results['neffc/nc'] = [divide(float(neff), float(ncr)) for neff, ncr
                           in zip(results['neffcross'], results['ncross'])]
    for flux, error in zip(results['runflux'], results['errflux']):
        results['pMD'].append(flux[-1] * time_step)
        results['1-p'].append(divide(float(1.0 - results['pMD'][-1]),
                                     results['pMD'][-1]))
        results['teffMD'].append(end_step * error[4]**2)
        results['corrMD'].append(
            divide(results['teffMD'][-1], results['1-p'][-1])
        )
    return results


def _effective_crossings(fluxdata, nint, end_step):
    """Analyse flux data and obtain effective crossings.

    Parameters
    ----------
    fluxdata : list of tuples of ints
        The list contains the data obtained from a ``md-flux``
        simulation.
    nint : int
        The number of interfaces used.
    end_step : int
        This is the last step done in the simulation.

    Returns
    -------
    eff_cross : list of lists
        `eff_cross[i]` is the effective crossings times for
        interface `i`.
    ncross : list of ints
        `ncross[i]` is the number of crossings for interface `i`.
    neffcross : list of ints
        `neffcross[i]` is the number of effective crossings for
        interface `i`.
    time_in_state : dict
        The time spent in the different states which are labeled with
        the keys 'A', 'B', 'OA', 'OB'. 'O' is taken to mean the
        'overall' state.

    Note
    ----
    We do here `intf - 1`. This is just to be compatible with the old
    FORTRAN code where the interfaces are numbered 1, 2, 3 rather than
    0, 1, 2. If this is to be changed in the future the `-1` can just
    be removed. Such a change will also require changes to the
    formatter for flux data.

    """
    # First line is used to determine if we start in B or A
    overallstate_a = not (fluxdata[0][1] == 2 and fluxdata[0][2] < 0)
    firstcross = [True] * nint
    ncross = [0] * nint
    neffcross = [0] * nint
    eff_cross = [[] for _ in range(nint)]
    end = {'A': 0, 'B': 0, 'OA': 0, 'OB': 0}
    start = {'A': 0, 'B': 0, 'OA': 0, 'OB': 0}
    time_in_state = {'A': 0, 'B': 0, 'OA': 0, 'OB': 0}
    time, intf, sign = None, None, None
    for (time, intf, sign) in fluxdata:
        if sign > 0:  # positive direction
            if intf - 1 == 0:  # moving out of A
                end['A'] = time
                time_in_state['A'] += (end['A'] - start['A'])
            elif intf - 1 == 2:  # moving into B
                start['B'] = time
                if overallstate_a:  # if we came from A
                    end['OA'] = time
                    start['OB'] = time
                    time_in_state['OA'] += (end['OA'] - start['OA'])
                    overallstate_a = False
            ncross[intf - 1] += 1
            if firstcross[intf - 1]:
                firstcross[intf - 1] = False
                neffcross[intf - 1] += 1
                eff_cross[intf - 1].append((time - time_in_state['OB'], time))
        elif sign < 0:
            if intf - 1 == 0:  # moving into A
                firstcross = [True] * nint
                start['A'] = time
                if not overallstate_a:  # if we came from B
                    end['OB'] = time
                    start['OA'] = time
                    time_in_state['OB'] += (end['OB'] - start['OB'])
                    overallstate_a = not overallstate_a
            elif intf - 1 == 2:  # moving out of B
                end['B'] = time
                time_in_state['B'] += (end['B'] - start['B'])
    # Now, just add up the remaining:
    state = 'OA' if overallstate_a else 'OB'
    time_in_state[state] += (end_step - start[state])
    if intf - 1 == 0 and sign < 0:
        # Note that the sign < 0 works for sign=None
        time_in_state['A'] += (end_step - start['A'])
    elif intf - 1 == 2 and sign > 0:
        # Note that the sign > 0 works for sign=None
        time_in_state['B'] += (end_step - start['B'])
    return eff_cross, ncross, neffcross, time_in_state


def _calculate_flux(effective_cross, time_in_state, time_window, time_step):
    """Calculate the flux in different time windows.

    Parameters
    ----------
    effective_cross : list
        The number of effective crossings, obtained from
        ``_effective_crossings``.
    time_in_state : int
        Time spent in over-all state ``A``.
    time_window : int
        This is the time window we consider for calculating the flux.
    time_step : float
        This is the time-step for the simulation.

    Returns
    -------
    time : np.array
        The times for which we have calculated the flux.
    ncross : np.array
        The number of crossings within a time window.
    flux : np.array
        The flux within a time window.

    """
    max_windows = int(1.0 * time_in_state / time_window)
    ncross = np.zeros(max_windows, dtype=int)
    for crossing in effective_cross:
        idx = int(np.floor((crossing[0] - 0.0) / time_window))
        idx = min(idx, max_windows - 1)
        ncross[idx] += 1
    flux = (1.0 * ncross) / (time_step * time_window)
    time = np.arange(1, max_windows+1) * time_window
    return time, ncross, flux


def find_crossings(order, interfaces):
    """Find crossings with interfaces for given order parameter data.

    Parameters
    ----------
    order : numpy.array (1D)
        Order parameters, as a function of time.
    interfaces : list of floats
        The interfaces for which we will investigate crossings.

    Returns
    -------
    out : list of tuple
        Each tuple contains the crossings on the following form:
        (step, interface-number, direction), where direction = '+' if
        the interface was crossed while moving to the right and '-' if
        the movement was towards the left.

    """
    leftside_prev = None  # previous
    cross = []
    for step, orderi in enumerate(order):
        leftside, ncross = check_crossing(step, orderi,
                                          interfaces, leftside_prev)
        leftside_prev = list(leftside)
        cross.extend(ncross)
    return cross
