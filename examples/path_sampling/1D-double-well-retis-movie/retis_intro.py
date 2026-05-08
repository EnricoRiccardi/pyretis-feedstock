# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This is a simple RETIS example.

Here we do nothing fancy, we will just get to know
some of the objects in PyRETIS.

Have fun!
"""
import colorama
import numpy as np
from pyretis.inout.settings import (fill_up_tis_and_retis_settings,
                                    _add_default_settings,
                                    _add_specific_default_settings)
from pyretis.core import System, create_box, Particles
from pyretis.initiation import initiate_path_simulation
from pyretis.setup import (
    create_force_field,
    create_orderparameter,
    create_simulation
)

# Let us define the simulation:
SETTINGS = {}
# Basic settings for the simulation:
SETTINGS['simulation'] = {
    'task': 'retis',
    'steps': 10,
    'interfaces': [-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, 1.0],
}
# Basic settings for the system:
SETTINGS['system'] = {'units': 'reduced', 'temperature': 0.07, 'dimensions': 1}
# Basic settings for the Langevin integrator:
SETTINGS['engine'] = {
    'class': 'Langevin',
    'gamma': 0.3,
    'high_friction': False,
    'seed': 0,
    'timestep': 0.002
}
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
    'class': 'Position',
    'dim': 'x',
    'index': 0,
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


def set_up_system(settings):
    """Set up the system.

    Parameters
    ----------
    settings : dict
        The settings required to set up the system.

    Returns
    -------
    sys : object like :py:class:`.System`
        A system object we can use in a simulation.

    """
    box = create_box(periodic=[False])
    sys = System(temperature=settings['system']['temperature'],
                 units=settings['system']['units'], box=box)
    sys.forcefield = create_force_field(settings)
    sys.order_function = create_orderparameter(settings)
    sys.particles = Particles(dim=1)
    sys.add_particle(np.array([-1.0]), mass=1, name='Ar', ptype=0)
    return sys


def print_step_results(ensembles, retis_result):
    """Print out RETIS results.

    Parameters
    ----------
    ensembles : list of dicts. They contain:
        *  PathEnsemble: objects like :py:class:`.PathEnsemble`
           The different path ensembles we are simulating.
    retis_result : list of lists
        The results of a RETIS simulation step.

    """
    for ensemble in ensembles:
        path_ensemble = ensemble['path_ensemble']
        name = path_ensemble.ensemble_name
        idx = path_ensemble.ensemble_number
        print(f'Move in {name}')
        accepted = retis_result[f'accept-{idx}']
        name_of_move = retis_result[f'move-{idx}']
        print(f'\tType: {name_of_move}')
        if name_of_move == 'swap':
            idx2 = retis_result[f'all-{idx}']['swap-with']
            name2 = ensembles[idx2]['path_ensemble'].ensemble_name
            print(f'\tSwapping: {name2} -> {name}')
        elif name_of_move == 'tis':
            trial_path = retis_result[f'path-{idx}']
            if trial_path.generated[0] == 'sh':
                tis_move = 'shooting'
            elif trial_path.generated[0] == 'tr':
                tis_move = 'time-reversal'
            else:
                tis_move = 'unknown'
            print(f'\tTIS move: {tis_move}')
        print(f'\tResult: {accepted}')


def main():
    """Just run the simulation."""
    colorama.init(autoreset=True)
    _add_default_settings(SETTINGS)
    _add_specific_default_settings(SETTINGS)
    SETTINGS['system']['obj'] = set_up_system(SETTINGS)
    fill_up_tis_and_retis_settings(SETTINGS)
    simulation = create_simulation(SETTINGS)
    print(simulation)
    print('# INITIATING TRAJECTORIES...')

    ensembles = simulation.ensembles
    for i, _ in enumerate(initiate_path_simulation(simulation, SETTINGS)):
        path_ensemble = ensembles[i]['path_ensemble']
        name = path_ensemble.ensemble_name
        print(f'Info about ensemble {name}:')
        print(path_ensemble)
        print('Info about the initial path:')
        print(path_ensemble.last_path)
        print('')
    # We can interact directly with points in trajectories,
    # here is a simple example
    path = ensembles[2]['path_ensemble'].last_path
    first = True
    for i, point in enumerate(path.phasepoints):
        order = point.order[0]
        pos = point.particles.pos
        vel = point.particles.vel
        if order > -0.8 and first:
            print('First crossing of -0.8 for [1^+]:')
            print(f'\tStep: {i}')
            print(f'\tlambda: {order}')
            print(f'\tPosition and velocity: {pos}  {vel}')
            first = False

    # Let us do one more step:
    print('Running a single RETIS step...')
    result = simulation.step()
    for ensemble in ensembles:
        path_ensemble = ensemble['path_ensemble']
        name = path_ensemble.ensemble_name
        idx = path_ensemble.ensemble_number
        print(f'Move in {name}')
        status = result[f'status-{idx}']
        accepted = result[f'accept-{idx}']
        name_of_move = result[f'move-{idx}']
        # `status` is equal to "ACC" if the move is accepted
        # otherwise it will be one of:
        # 'MCR': 'Momenta change rejection',
        # 'BWI': 'Backward trajectory end at wrong interface',
        # 'BTL': 'Backward trajectory too long (detailed balance condition)',
        # 'BTX': 'Backward trajectory too long (max-path exceeded)',
        # 'KOB': 'Kicked outside of boundaries',
        # 'FTL': 'Forward trajectory too long (detailed balance condition)',
        # 'FTX': 'Forward trajectory too long (max-path exceeded)',
        # 'NCR': 'No crossing with middle interface'
        print(f'\tType: {name_of_move}')
        if name_of_move == 'swap':
            print(status, accepted)
            # If this is the case, the result is on the form
            # [move, accepted, .., swap-with] where swap-with is the
            # ensemble we are trying to swap with.
            idx2 = result[f'all-{idx}']['swap-with']
            name2 = ensembles[idx2]['path_ensemble'].ensemble_name
            print(f'\tSwapping: {name2} -> {name}')
        elif name_of_move == 'tis':
            trial_path = result[f'path-{idx}']
            if trial_path.generated[0] == 'sh':
                tis_move = 'shooting'
            elif trial_path.generated[0] == 'tr':
                tis_move = 'time-reversal'
            else:
                tis_move = 'unknown'
            print(f'\tTIS move: {tis_move}')
        print(f'\tResult: {status}')

    # Run the rest of the simulation:
    while not simulation.is_finished():
        result = simulation.step()
        print(f"Simulation step: {result['cycle']['step']}")
        print_step_results(ensembles, result)
        print('')


if __name__ == '__main__':
    main()
