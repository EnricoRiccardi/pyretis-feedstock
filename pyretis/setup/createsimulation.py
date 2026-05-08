# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This module handles the creation of simulations from settings.

The different simulations are defined as objects which inherit
from the base :py:class:`.Simulation` class defined in
:py:mod:`pyretis.simulation.simulation`. Here, we are treating
each simulation with a special case.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

create_simulation (:py:func:`.create_simulation`)
    Method for creating a simulation object from settings.

create_ensemble (:py:func:`.create_ensemble`)
    Method for creating an ensemble dictionary from settings.

create_ensembles (:py:func:`.create_ensembles`)
    Method for creating a list of ensemble from settings.

create_nve_simulation (:py:func:`.create_mve_simulation`)
    Method for creating a nve simulation object from settings.

create_mdflux_simulation (:py:func:`.create_mdflux_simulation`)
    Method for creating a mdflux simulation object from settings.

create_umbrellaw_simulation (:py:func:`.create_umbrellaw_simulation`)
    Method for creating a umbrellaw simulation object from settings.

create_tis_simulation (:py:func:`.create_tis_simulation`)
    Method for creating a tis simulation object from settings.

create_retis_simulation (:py:func:`.create_retis_simulation`)
    Method for creating a retis simulation object from settings.

prepare_system (:py:func:`.prepare_system`)
    Method for creating a system object from settings.

prepare_engine (:py:func:`.prepare_engine`)
    Method for creating an engine object from settings.

"""
import logging
import os
import colorama
from pyretis.core.pathensemble import get_path_ensemble_class
from pyretis.setup.createforcefield import create_force_field
from pyretis.inout.settings import (add_default_settings,
                                    add_specific_default_settings,
                                    settings_from_restart)
from pyretis.setup.common import create_orderparameter
from pyretis.core.units import units_from_settings
from pyretis.core.random_gen import create_random_generator
from pyretis.simulation.md_simulation import (
    SimulationMD,
    SimulationNVE,
    SimulationMDFlux,
)
from pyretis.simulation.mc_simulation import UmbrellaWindowSimulation
from pyretis.simulation.path_simulation import (
    SimulationTIS,
    SimulationRETIS
)
from pyretis.setup import (
    create_system,
    create_engine
)
from pyretis.inout.checker import (
    check_engine,
    check_ensemble,
)
from pyretis.inout.screen import REFERENCE

logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = ['create_simulation', 'create_ensemble', 'create_ensembles',
           'create_nve_simulation', 'create_mdflux_simulation',
           'create_umbrellaw_simulation', 'create_tis_simulation',
           'create_retis_simulation', 'prepare_system', 'prepare_engine']


def create_ensemble(settings):
    """Create the path ensemble from (ensemble) simulation settings.

    Parameters
    ----------
    settings : dict
        This dict contains the settings needed to create the path
        ensemble.

    Returns
    -------
    ensemble : dict of
        objects that contains all the information needed in the ensemble.

    """
    i_ens = settings['tis']['ensemble_number']

    rgen_ens = create_random_generator(settings['tis'])
    rgen_path = create_random_generator(settings['system'])

    system = prepare_system(settings)
    engine = prepare_engine(settings)
    klass = get_path_ensemble_class(settings['particles']['type'])
    interfaces = settings['simulation']['interfaces']
    exe_dir = settings['simulation'].get('exe_path', os.path.abspath('.'))
    path_ensemble = klass(i_ens, interfaces, rgen=rgen_path, exe_dir=exe_dir)
    logger.log(REFERENCE, '\n--- Creating ensemble: %s ---',
               path_ensemble.ensemble_name)

    # for PPRETIS / PPTIS: correct the starting condition of the paths
    if settings['simulation']['task'] in ['pptis', 'repptis']:
        change_start_cond = True
        # but skip this change if we are dealing with the [0-] ensemble
        if i_ens == 0:  # the first ensemble
            if settings['simulation'].get('flux', False) or \
                    settings['simulation'].get('zero_left', False):
                change_start_cond = False

        if change_start_cond:
            # adapt the starting condition of the paths
            logger.info("adapted start_condition to R/L in ensemble %s", i_ens)
            path_ensemble.start_condition = ['R', 'L']
            path_ensemble.must_cross_M = True  # used in the shooting move

    order_function = create_orderparameter(settings)
    # code for computing the permeability
    if i_ens == 0 and settings['simulation'].get("permeability", False):
        # the starting condition in [0-'] can be R or L
        path_ensemble.start_condition = ['R', 'L']
        # check mirror move settings
        if hasattr(order_function, 'mirror_pos'):
            offset = getattr(order_function, 'offset', 0)
            moved_interfaces = [i - offset for i in interfaces]
            left, right = moved_interfaces[0], moved_interfaces[2]
            correct_mirror = (left+right)/2.
            if abs(order_function.mirror_pos-correct_mirror) > 1E-5:
                msg = "Order function should have a mirror at "
                msg += f"{correct_mirror}, found one at "
                msg += f"{order_function.mirror_pos} instead."
                raise ValueError(msg)
    # Add check to see if mirror makes sense
    engine.can_use_order_function(order_function)
    ensemble = {'engine': engine,
                'system': system,
                'order_function': order_function,
                'interfaces': interfaces,
                'exe_path': exe_dir,
                'path_ensemble': path_ensemble,
                'rgen': rgen_ens}
    ensemble['system'].order = engine.calculate_order(ensemble)

    return ensemble


def create_ensembles(settings):
    """Create a list of dictionary from ensembles simulation settings.

    Parameters
    ----------
    settings : dict
        This dict contains the settings needed to create the path
        ensemble.

    Returns
    -------
    ensembles : list of dicts
        List of ensembles.

    """
    # Example:
    # ensembles, assuming len(interfaces) = 3
    # (RE)TIS:   flux=False, zero_ensemble=False  [1+]
    # (RE)TIS:   flux=False, zero_ensemble=True   [0+], [1+]
    # (RE)TIS:   flux=True,  zero_ensemble=True   [0-], [0+], [1+]
    # /      :   flux=True,  zero_ensemble=False  doesn't make sense
    # so number of ensembles can be 1, 2, or 3
    ensembles = []
    intf = settings['simulation']['interfaces']
    # for PPTIS or REPPTIS: memory
    mem = settings.get(settings['simulation']['task'], {}).get('memory', 1)

    j_ens = 0
    # add [0-] if it exists
    if settings['simulation'].get('flux', False) or \
            settings['simulation'].get('zero_left', False):
        reactant, middle, product = float('-inf'), intf[0], intf[0]
        if settings['simulation'].get('zero_left', False):
            reactant = settings['simulation']['zero_left']
            if settings['simulation'].get('permeability', False):
                middle = (settings['simulation']['zero_left'] + intf[0]) / 2
        settings['ensemble'][j_ens]['simulation']['interfaces'] =\
            [reactant, middle, product]
        j_ens += 1

    # add [0+] if it exists
    if settings['simulation'].get('zero_ensemble', True):
        reactant, middle, product = intf[0], intf[0], intf[-1]

        # for PPTIS or REPPTIS: overwrite 'product' interface
        if settings['simulation']['task'] in ['pptis', 'repptis']:
            # mem is an integer, DEFAULT is 1
            # intf cannot be more to the right than intf[-1]
            product = intf[min(mem, len(intf)-1)]

        settings['ensemble'][j_ens]['simulation']['interfaces'] = \
            [reactant, middle, product]
        j_ens += 1

    # j_ens is the number of ensembles that is skipped

    # set interfaces and set detect for [1+], [2+], ...
    reactant, product = intf[0], intf[-1]
    for i, i_ens in enumerate(range(j_ens, len(settings['ensemble']))):
        middle = intf[i + 1]    # the lambda_i interface

        # for PPTIS / REPPTIS: overwrite 'reactant' and 'product' intf
        if settings['simulation']['task'] in ['pptis', 'repptis']:
            # mem is an integer, DEFAULT is 1
            # intf cannot be more to the left than intf[0]
            reactant = intf[max(0, i + 1 - mem)]
            # intf cannot be more to the right than intf[-1]
            product = intf[min(i + 1 + mem, len(intf)-1)]
            # middle is not changed, actually

        settings['ensemble'][i_ens]['simulation']['interfaces'] = \
            [reactant, middle, product]
        settings['ensemble'][i_ens]['tis']['detect'] = intf[i + 2]  # next intf

    # create all ensembles
    for i in range(len(settings['ensemble'])):
        ensembles.append(create_ensemble(settings['ensemble'][i]))

    return ensembles


def create_nve_simulation(settings):
    """Set up and create a NVE simulation.

    Parameters
    ----------
    settings : dict
        The settings needed to set up the simulation.

    Returns
    -------
    SimulationNVE : object like :py:class:`.SimulationNVE`
        The object representing the simulation to run.

    """
    ensemble = {'engine': prepare_engine(settings),
                'system': prepare_system(settings),
                'rgen': create_random_generator(settings['simulation'])}

    key_check('steps', settings)

    controls = {'steps': settings['simulation']['steps'],
                'startcycle': settings['simulation'].get('startcycle', 0)}

    return SimulationNVE(ensemble, settings, controls)


def create_md_simulation(settings):
    """Set up and create a generic MD simulation.

    Parameters
    ----------
    settings : dict
        The settings needed to set up the simulation.

    Returns
    -------
    out : object like :py:class:`.SimulationMD`
        The object representing the simulation to run.

    """
    ensemble = {'engine': prepare_engine(settings),
                'system': prepare_system(settings),
                'order_function': create_orderparameter(settings)}

    key_check('steps', settings)

    controls = {'steps': settings['simulation']['steps'],
                'startcycle': settings['simulation'].get('startcycle', 0)}

    return SimulationMD(ensemble, settings, controls)


def create_mdflux_simulation(settings):
    """Set up and create a MD FLUX simulation.

    Parameters
    ----------
    settings : dict
        The settings needed to set up the simulation.

    Returns
    -------
    SimulationMDFlux : object like :py:class:`.SimulationMDFlux`
        The object representing the simulation to run.

    """
    engine = prepare_engine(settings)
    order_function = create_orderparameter(settings)
    ensemble = {'system': prepare_system(settings),
                'engine': engine,
                'order_function': order_function,
                'rgen': create_random_generator(settings['simulation'])}
    engine.can_use_order_function(order_function)
    for key in ('steps', 'interfaces'):
        key_check(key, settings)
    controls = {'steps': settings['simulation']['steps'],
                'startcycle': settings['simulation'].get('startcycle', 0)}
    return SimulationMDFlux(ensemble, settings, controls)


def create_umbrellaw_simulation(settings):
    """Set up and create a Umbrella Window simulation.

    Parameters
    ----------
    settings : dict
        The settings needed to set up the simulation.
        Note that mincycle is trasmitted as steps to the object,
        but is has a different meaning than the other
        simulations.

    Returns
    -------
    UmbrellaWindowSimulation : obj like :py:class:`.UmbrellaWindowSimulation`
        The object representing the simulation to run.

    """
    ensemble = {'system': prepare_system(settings),
                'rgen': create_random_generator(settings['simulation'])}

    for key in ('umbrella', 'overlap', 'maxdx', 'mincycle'):
        key_check(key, settings)

    controls = {'steps': settings['simulation']['mincycle'],
                'startcycle': settings['simulation'].get('startcycle', 0)}

    return UmbrellaWindowSimulation(ensemble, settings, controls)


def create_tis_simulation(settings):
    """Set up and create a single TIS simulation.

    Parameters
    ----------
    settings : dict
        The settings needed to set up the simulation.

    Returns
    -------
    SimulationTIS : object like :py:class:`.SimulationSingleTIS`
        The object representing the simulation to run.

    """
    check_ensemble(settings)
    ensembles = create_ensembles(settings)

    key_check('steps', settings)

    controls = {'rgen': create_random_generator(settings['simulation']),
                'steps': settings['simulation']['steps'],
                'startcycle': settings['simulation'].get('startcycle', 0)}

    return SimulationTIS(ensembles, settings, controls)


def create_retis_simulation(settings):
    """Set up and create a RETIS simulation.

    Parameters
    ----------
    settings : dict
        This dictionary contains the settings for the simulation.

    Returns
    -------
    SimulationRETIS : object like :py:class:`.SimulationRETIS`
        The object representing the simulation to run.

    """
    check_ensemble(settings)
    ensembles = create_ensembles(settings)

    key_check('steps', settings)

    controls = {'rgen': create_random_generator(settings['simulation']),
                'steps': settings['simulation']['steps'],
                'startcycle': settings['simulation'].get('startcycle', 0)}

    return SimulationRETIS(ensembles, settings, controls)


def prepare_system(settings):
    """Create a system from given settings.

    Parameters
    ----------
    settings : dict
        This dictionary contains the settings for the simulation.

    Returns
    -------
    syst : object like :py:class:`.syst`
        This object will correspond to the selected simulation type.

    """
    if settings.get('system', {}).get('obj', False):
        return settings['system']['obj']

    logtxt = 'Initializing unit system.'
    logger.progress(logtxt)
    units_from_settings(settings)

    logtxt = 'Creating system from settings.'
    logger.progress(logtxt)
    system = create_system(settings)

    logtxt = 'Creating force field.'
    logger.progress(logtxt)
    system.forcefield = create_force_field(settings)
    system.particles.vpot = system.evaluate_potential()

    system.extra_setup()
    return system


def prepare_engine(settings):
    """Create an engine from given settings.

    Parameters
    ----------
    settings : dict
        This dictionary contains the settings for the simulation.

    Returns
    -------
    engine : object like :py:class:`.engine`
        This object will correspond to the selected simulation type.

    """
    if settings.get('engine', {}).get('obj', False):
        return settings['engine']['obj']

    logtxt = units_from_settings(settings)
    logger.progress(logtxt)

    check_engine(settings)
    engine = create_engine(settings)
    logtxt = f'Created engine "{engine}" from settings.'
    logger.progress(logtxt)
    return engine


def key_check(key, settings):
    """Check for the presence of a key in settings."""
    # todo These checks shall be done earlier, when cleaning the input.
    if key not in settings['simulation']:
        msgtxt = f'Simulation setting "{key}" is missing!'
        logger.critical(msgtxt)
        raise ValueError(msgtxt)


def create_simulation(settings):
    """Create simulation(s) from given settings.

    This function will set up some common simulation types.
    It is meant as a helper function to automate some very common set-up
    task. It will here check what kind of simulation we are to perform
    and then call the appropriate function for setting that type of
    simulation up.

    Parameters
    ----------
    settings : dict
        This dictionary contains the settings for the simulation.

    Returns
    -------
    simulation : object like :py:class:`.Simulation`
        This object will correspond to the selected simulation type.

    """
    sim_type = settings['simulation']['task'].lower()

    sim_map = {
        'md': create_md_simulation,
        'md-nve': create_nve_simulation,
        'md-flux': create_mdflux_simulation,
        'umbrellawindow': create_umbrellaw_simulation,
        'make-tis-files': create_tis_simulation,
        'tis': create_tis_simulation,
        'explore': create_tis_simulation,
        'retis': create_retis_simulation,
        'pptis': create_tis_simulation,
        'repptis': create_retis_simulation,
    }

    # Improve setting quality
    add_default_settings(settings)
    add_specific_default_settings(settings)

    if settings['simulation'].get('restart', False):
        settings, info_restart = settings_from_restart(settings)

    if sim_type not in sim_map:  # TODO put in check_sim_type
        msgtxt = f'Unknown simulation task {sim_type}'
        logger.error(msgtxt)
        raise ValueError(msgtxt)

    simulation = sim_map[sim_type](settings)
    msgtxt = f'{simulation}'
    lines = msgtxt.split('\n')
    sim_name = (colorama.Fore.BLUE + lines[0] +
                colorama.Fore.WHITE)
    rest = '\n'.join(lines[1:])
    logger.reference('\nCreated simulation:\n%s\n%s', sim_name, rest)

    if settings['simulation'].get('restart', False):
        simulation.load_restart_info(info_restart)

    return simulation
