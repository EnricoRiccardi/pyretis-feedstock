# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This is a simple RETIS example.

This example will just consider a simple 1D potential where we
aim to calculate the crossing probability and the rate constant.

Have fun!
"""
import colorama
import numpy as np
from tqdm import tqdm
from pyretis.initiation import initiate_path_simulation
from pyretis.core.properties import Property
from pyretis.inout.settings import (fill_up_tis_and_retis_settings,
                                    add_default_settings,
                                    add_specific_default_settings)
from pyretis.setup import create_simulation
from pyretis.analysis.path_analysis import _pcross_lambda_cumulative

np.set_printoptions(legacy='1.21')

INTERFACES = [-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, 1.0]
# Let us define the simulation:
SETTINGS = {}
# Basic settings for the simulation:
SETTINGS['simulation'] = {
    'task': 'retis',
    'steps': 150,
    'interfaces': INTERFACES
}
# Basic settings for the simulation box:
SETTINGS['box'] = {'cell': [3]}

# Basic settings for the system:
SETTINGS['system'] = {'units': 'lj',
                      'temperature': 0.09,
                      'dimensions': 1}

# Basic settings for the Langevin integrator:
SETTINGS['engine'] = {
    'class': 'Langevin',
    'gamma': 0.3,
    'high_friction': False,
    'seed': 0,
    'timestep': 0.002
}
# Particles:
SETTINGS['particles'] = {'mass': {'Ar': 1.0},
                         'name': ['Ar'],
                         'ptype': [0],
                         'position': {'generate': '1d'}}

# Potential parameters:
# The potential is: `V_\text{pot} = a x^4 - b (x - c)^2`
SETTINGS['potential'] = [
    {
        'class': 'DoubleWell',
        'a': 1.0, 'b': 2.0, 'c': 0.0,
    }
]
# Settings for the order parameter:
SETTINGS['orderparameter'] = {
    'class': 'PositionVelocity',
    'dim': 'x', 'index': 0,
    'periodic': False
}
# TIS specific settings:
SETTINGS['tis'] = {
    'freq': 0.5,
    'maxlength': 20000,
    'aimless': True,
    'allowmaxlength': False,
    'sigma_v': -1,
    'seed': 0,
    'zero_momentum': False,
    'rescale_energy': False
}
SETTINGS['initial-path'] = {'method': 'kick'}
# RETIS specific settings:
SETTINGS['retis'] = {
    'swapfreq': 0.5,
    'relative_shoots': None,
    'nullmoves': True,
    'swapsimul': True
}

# For convenience:
TIMESTEP = SETTINGS['engine']['timestep']
ANALYSIS = {'ngrid': 100, 'nblock': 5}
add_default_settings(SETTINGS)
add_specific_default_settings(SETTINGS)
fill_up_tis_and_retis_settings(SETTINGS)


def step_txt(ensembles, retis_result, prun):
    """Return some text information about the RETIS step.

    Parameters
    ----------
    ensembles : list of dicts, it contains:
       *  path_ensemble: objects like :py:class:`.PathEnsemble`
               The different path ensembles we are simulating.
    retis_result : list of lists
        The results from a RETIS simulation step.
    prun : list of floats
        The running average for the ensemble crossing probabilities.

    Returns
    -------
    out : list of strings
        Each string contains information about the move performed
        in an ensemble.

    """
    txt = []
    force = 0  # counter for force evaluations
    # The force counter will just count the number of
    # force evaluations for generating a path, it will not include
    # the force evaluation that might be performed prior to propagation,
    # i.e. the initial force evaluation since this can be stored in memory,
    # along the path -> i.e. we can avoid it.
    for ensemble in ensembles:
        path_ensemble = ensemble['path_ensemble']
        name = path_ensemble.ensemble_name
        idx = path_ensemble.ensemble_number
        name_of_move = retis_result[f'move-{idx}']
        accepted = retis_result[f'accept-{idx}']
        line = []
        if name_of_move == 'swap':
            idx2 = retis_result[f'all-{idx}']['swap-with']
            name2 = ensembles[idx2]['path_ensemble'].ensemble_name
            move = f'{name_of_move} {name2},'
            if idx == 0 or (idx == 1 and idx2 == 0):
                # Evaluate forces when swapping [0^-] <-> [0^+]:
                force += path_ensemble.paths[-1]['length'] - 2
        elif name_of_move == 'tis':
            trial_path = retis_result[f'path-{idx}']
            tis_move = trial_path.generated[0]
            move = f'{name_of_move} ({tis_move}),'
            if tis_move == 'sh':
                force += path_ensemble.paths[-1]['length'] - 1
        else:
            move = f'{name_of_move},'
        line.append(f'{name}: {move:11s}')
        line.append(f'{accepted},')
        if idx > 0:
            line.append(f'p = {prun[idx]:<8.6g}')
        txt.append(' '.join(line))
    return txt, force


def probability_path_ensemble(path_ensemble, step, prun, orderp):
    """Update running estimates of probabilities for an ensemble.

    Parameters
    ----------
    path_ensemble : objects like :py:class:`.PathEnsemble`
        The ensemble to analyse.
    step : int
        The current simulation step.
    prun : float
        Current running estimate for the crossing probability
        for this ensemble.
    orderp : numpy.array
        The maximum order parameters for all accepted paths.

    Returns
    -------
    prun : float
        Updated running average for the crossing probability.
    ordernew : numpy.array
        Updated list of all order parameters seen.
    lamb : numpy.array
        Coordinates used for obtaining the crossing probability
        as a function of the order parameter.
    pcross : numpy.array
        The crossing probability as a function of the order
        parameter.

    """
    ordermax = path_ensemble.last_path.ordermax[0]
    success = 1 if ordermax > path_ensemble.detect else 0
    if prun is None:
        prun = success
    else:
        npath = step + 1
        prun = float(success + prun * (npath - 1)) / float(npath)
    if orderp is None:
        ordernew = np.array([ordermax])
    else:
        ordernew = np.append(orderp, ordermax)
    ordermax = min(max(ordernew), max(path_ensemble.interfaces))
    ordermin = path_ensemble.interfaces[1]
    pcross, lamb = _pcross_lambda_cumulative(ordernew, ordermin, ordermax,
                                             ANALYSIS['ngrid'])
    return prun, ordernew, lamb, pcross


def simple_block(data, nblock):
    """Block error analysis for running average."""
    length = len(data)
    if length < nblock:
        return float('inf')
    skip, _ = divmod(len(data), nblock)
    run_avg = []
    block_avg = []
    for i in range(nblock):
        idx = (i + 1) * skip - 1
        run_avg.append(data[idx])
        if i == 0:
            block_avg.append(run_avg[i])
        else:
            block_avg.append((i + 1)*run_avg[i] - i*run_avg[i-1])
    var = sum([(x - data[-1])**2 for x in block_avg]) / (nblock - 1)
    err = np.sqrt(var / nblock)
    return err


def analyse_path_ensembles(ensembles, step, variables):
    """Calculate some properties for the path ensembles.

    Here we obtain the crossing probabilities and the initial flux.

    Parameters
    ----------
    ensembles : a list of dicts. They contain:
      *  objects like :py:class:`.PathEnsemble`
             These objects contain information about accepted paths for the
             different path ensembles.
    step : integer
        The current step number.
    variables : dict of objects
        This dict contains variables we can use for updating derived data
        for the different path ensembles.

    Returns
    -------
    out : dict of floats
        A dictionary with computed initial flux, crossing probability and
        errors in these.

    """
    calc = {'flux': 0, 'fluxe': 0,
            'pcross': 0, 'pcrosse': float('inf'),
            'kab': 0, 'kabe': float('inf')}
    length0 = variables['length0']
    length1 = variables['length1']
    length0.add(ensembles[0]['path_ensemble'].last_path.length)
    length1.add(ensembles[1]['path_ensemble'].last_path.length)
    flux = 1.0 / ((length0.mean + length1.mean - 4.0) * TIMESTEP)
    calc['flux'] = flux
    variables['flux'].append(flux)
    fluxe = simple_block(variables['flux'], ANALYSIS['nblock'])
    calc['fluxe'] = fluxe
    accprob = 1.0
    matched_prob = []
    matched_lamb = []
    for i, ensemble in enumerate(ensembles):
        if i == 0:
            continue
        prun = variables['prun'][i]
        orderp = variables['orderp'][i]
        prun, ordernew, lamb, pcross = probability_path_ensemble(
            ensemble['path_ensemble'], step, prun, orderp)
        variables['orderp'][i] = ordernew
        variables['prun'][i] = prun
        variables['pcross'][i] = pcross
        variables['lamb'][i] = lamb
        idx = np.where(lamb <= ensemble['path_ensemble'].detect)[0]
        matched_lamb.extend(lamb[idx])
        matched_prob.extend(pcross[idx] * accprob)
        accprob *= prun
    if matched_lamb[-1] < INTERFACES[-1]:
        pcross = 0
        pcrosse = float('inf')
    else:
        pcross = matched_prob[-1]
        variables['match'].append(pcross)
        pcrosse = simple_block(variables['match'], ANALYSIS['nblock'])
    calc['pcross'] = pcross
    calc['pcrosse'] = pcrosse
    calc['kab'] = flux * pcross
    if pcross == 0 or flux == 0:
        calc['kabe'] = float('inf')
    else:
        calc['kabe'] = calc['kab'] * np.sqrt(pcrosse**2 / pcross**2 +
                                             fluxe**2 / flux**2)
    return calc


def main():
    """Just run the simulation."""
    colorama.init(autoreset=True)
    simulation = create_simulation(SETTINGS)
    print(simulation)
    print('# INITIATING TRAJECTORIES...')
    print('# GENERATING INITIAL PATHS')
    for i, _ in enumerate(initiate_path_simulation(simulation, SETTINGS)):
        path_ensemble = simulation.ensembles[i]['path_ensemble']
        name = path_ensemble.ensemble_name
        print(f'Info about ensemble {name}:')
        print(path_ensemble)
        print('Info about the initial path:')
        print(path_ensemble.last_path)
        print('')

    # We make a dictionary of these variable for easier access:
    # Set up some variables for storing results:
    variables = {'length0': Property('Path length in [0^-]'),
                 'length1': Property('Path length in [0^+]'),
                 'flux': [],
                 'match': [],
                 'orderp': [None for _ in INTERFACES],
                 'prun': [None for _ in INTERFACES],
                 'lamb': [None for _ in INTERFACES],
                 'pcross': [None for _ in INTERFACES]}
    # Let us look at the resulting path ensembles and paths:
    ensembles = simulation.ensembles
    analyse_path_ensembles(ensembles, 0, variables)
    for ensemble in ensembles:
        path_ensemble = ensemble['path_ensemble']
        name = path_ensemble.ensemble_name
        print(f'Info about ensemble {name}:')
        print(path_ensemble)
        print('Info about the initial path:')
        print(path_ensemble.last_path)
        print('')
    # Run the rest of the simulation.
    ftot = 0
    for result in tqdm(simulation.run(), initial=1,
                       total=SETTINGS['simulation']['steps']):
        step = result['cycle']['step']
        print(f'# Current cycle: {step}')
        anr = analyse_path_ensembles(ensembles, step, variables)
        retis_txt, force = step_txt(ensembles, result,
                                    variables['prun'])
        ftot += force
        for line in retis_txt:
            print(f'# {line}')
        print('# Flux: {flux:<8.6g} +- {fluxe:<8.6g}'.format(**anr))
        print(('# Crossing probability: {pcross:<8.6g} +-'
               '{pcrosse:<8.6g}').format(**anr))
        print('# K_AB: {kab:<8.6g} +- {kabe:<8.6g}'.format(**anr))
        print(f'# No. of force evaluations: {force:g}')
        print('')
    print(f'# Total number of force evaluations: {ftot}')


if __name__ == '__main__':
    main()
