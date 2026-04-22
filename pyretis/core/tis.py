# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This file contains functions used in TIS.

This module defines the functions needed to perform TIS simulations
(algorithms implemented as described by van Erp et al. [TIS]_ ) with
advanced shooting moves, such as Stone Skipping and Web Throwing
(algorithms implemented as described by Riccardi et al. [SS+WT]_ ).

Methods defined here
~~~~~~~~~~~~~~~~~~~~

make_tis_step (:py:func:`.make_tis_step`)
    A method that will perform a single TIS step.

make_tis_step_ensemble (:py:func:`.make_tis_step_ensemble`)
    A method to perform a TIS step for a path ensemble. It will handle
    adding of the path to a path ensemble object.

shoot (:py:func:`.shoot`)
    A method that will perform a shooting move.

select_shoot (:py:func:`.select_shoot`)
    A method that randomly selects the shooting point.

wire_fencing (:py:func:`.wire_fencing`)
    A method that will do a wire fencing shooting move.

stone_skipping (:py:func:`.stone_skipping`)
    A method that will do a stone skipping shooting move.

web_throwing (:py:func:`.web_throwing`)
    A method that will do a web throwing shooting move.

extender (:py:func:`.extender`)
    A method that will extend a path to target interfaces.

time_reversal (:py:func:`.time_reversal`)
    A method for performing the time reversal move.

ss_wt_wf_acceptance (:py:func:`.ss_wt_wf_acceptance`)
    A method to check the acceptance of a newly generated path
    according to super detailed balance rules for Stone Skipping,
    Wire Fencing (basic and High Acceptance version) and Web Throwing.


References
~~~~~~~~~~
.. [SS+WT] E. Riccardi, O. Dahlen, T. S. van Erp,
   J. Phys. Chem. letters, 8, 18, 4456, (2017),
   https://doi.org/10.1021/acs.jpclett.7b01617

.. [TIS] T. S. van Erp, D. Moroni, P. G. Bolhuis,
   J. Chem. Phys. 118, 7762 (2003),
   https://dx.doi.org/10.1063%2F1.1562614

"""
import logging
from pyretis.core.common import (null_move,
                                 counter,
                                 relative_shoots_select,
                                 segments_counter,
                                 compute_weight,
                                 select_and_trim_a_segment,
                                 crossing_counter,
                                 priority_checker,
                                 crossing_finder,
                                 wirefence_weight_and_pick)

from pyretis.core.montecarlo import metropolis_accept_reject
from pyretis.core.path import Path, paste_paths
from pyretis.inout.restart import write_ensemble_restart

logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())

__all__ = ['extender',
           'make_tis',
           'make_tis_step',
           'make_tis_step_ensemble',
           'paste_paths',
           'select_shoot',
           'shoot',
           'ss_wt_wf_acceptance',
           'wire_fencing',
           'stone_skipping',
           'time_reversal',
           'web_throwing']


def _order_value(order):
    """Return the first order parameter component (scalar or vector)."""
    try:
        return order[0]
    except (TypeError, IndexError):
        return order


def make_tis(ensembles, rgen, settings, cycle):
    """Perform a TIS step in a given path ensemble.

    This function will execute the TIS steps. When used in the RETIS method,
    three different options can be given:

    1) If `relative_shoots` is given in the input settings as a list of
       probability for each ensemble to be picked, then an ensemble
       will be ranomly picked and TIS performed on. For all the
       other ensembles we again have two options based on the given
       `settings['nullmoves']`:

       a) Do a 'null move' in all other ensembles.
       b) Do nothing for all other ensembles.

       Performing the null move in an ensemble will simply just accept
       the previously accepted path in that ensemble again.

    2) If `relative_shoots` is not given in the input settings, TIS moves
       will be executed for all path ensembles.

    Parameters
    ----------
    ensembles : list of dictionaries of objects
        This is a list of the ensembles we are using in the TIS method
        with their properties.
        They contain:

        * `path_ensemble`: object like :py:class:`.PathEnsemble`
          This is used for storing results for the simulation. It
          is also used for defining the interfaces for this simulation.
        * `system`: object like :py:class:`.System`
          System is used here since we need access to the temperature
          and to the particle list
        * `order_function`: object like :py:class:`.OrderParameter`
          The class used for calculating the order parameter(s).
        * `engine`: object like :py:class:`.EngineBase`
          The engine to use for propagating a path.

    rgen : object like :py:class:`.RandomGenerator`
        This is the random generator that will be used.
    settings : dict
        This dictionary contains the TIS settings.
    cycle : int
        The current cycle number.

    Yields
    ------
    out : dict
        This dictionary contains the results of the TIS move(s).

    """
    relative = settings.get('retis', {}).get('relative_shoots', None)
    nullmoves = settings.get('retis', {}).get('nullmoves', False)

    if relative is not None:
        idx, ensemble = relative_shoots_select(ensembles, rgen, relative)
        accept, trial, status = make_tis_step_ensemble(
            ensemble, settings['ensemble'][idx], cycle
        )
        result = {
            'ensemble_number': idx,
            'mc-move': trial.get_move(),
            'status': status,
            'trial': trial,
            'accept': accept
        }
        yield result
        # Do null moves in the other ensembles, if requested:
        if nullmoves:
            for ensemble in ensembles:
                other = ensemble['path_ensemble'].ensemble_number
                if other != idx:
                    accept, trial, status = null_move(ensemble, cycle)
                    result = {
                        'ensemble_number': other,
                        'mc-move': 'nullmove',
                        'status': status,
                        'trial': trial,
                        'accept': accept,
                    }
                    if 0 == cycle % settings['ensemble'][idx].get(
                            'output', {}).get('restart-file', 1):
                        write_ensemble_restart(ensemble, settings)
                    yield result
    else:  # Do TIS in all ensembles:
        prio_skip = priority_checker(ensembles, settings)
        for i, (ens, e_set) in enumerate(zip(ensembles, settings['ensemble'])):
            if prio_skip[i]:
                # Here we skip the ensembles with higher cycle numbers.
                yield
            else:
                accept, trial, status = make_tis_step_ensemble(ens,
                                                               e_set,
                                                               cycle)
                result = {
                    'ensemble_number': ens['path_ensemble'].ensemble_number,
                    'mc-move': trial.get_move(),
                    'status': status,
                    'trial': trial,
                    'accept': accept
                }
                yield result


def make_tis_step_ensemble(ensemble, settings, cycle):
    """Perform a TIS step for a given path ensemble.

    This function will run `make_tis_step` for the given path_ensemble.
    It will handle adding of the path. This function is intended for
    convenience when working with path ensembles. If we are using the
    path ensemble ``[0^-]`` then the start condition should be 'R' for
    right.

    Parameters
    ----------
    ensemble : dictionary of objects
        It contains:

        * `path_ensemble`: object like :py:class:`.PathEnsemble`
          This is used for storing results for the simulation. It
          is also used for defining the interfaces for this simulation.
        * `system`: object like :py:class:`.System`
          System is used here since we need access to the temperature
          and to the particle list.
        * `order_function`: object like :py:class:`.OrderParameter`
          The class used for obtaining the order parameter(s).
        * `engine`: object like :py:class:`.EngineBase`
          The engine to use for propagating a path.
        * `rgen`: object like :py:class:`.RandomGenerator`
          This is the random generator that will be used.

    settings : dict
        This dictionary contains the ensemble settings.
    cycle : int
        The current cycle number.

    Returns
    -------
    out[0] : boolean
        True if the new path can be accepted, False otherwise.
    out[1] : object like :py:class:`.PathBase`
        The generated path.
    out[2] : string
        The status of the path.

    """
    path_ensemble = ensemble['path_ensemble']
    engine = ensemble['engine']
    start_cond = path_ensemble.start_condition
    logger.info('TIS move in: %s', path_ensemble.ensemble_name)
    engine.exe_dir = path_ensemble.directory['generate']
    accept, trial, status = make_tis_step(ensemble,
                                          settings['tis'],
                                          start_cond)

    # If we are exploring, we accept everything.
    if settings['simulation']['task'] == 'explore':
        if trial.length > 2:
            accept, status = True, 'EXP'
            trial.status = 'EXP'

    if accept:
        logger.info('The move was accepted!')
    else:
        logger.info('The move was rejected! (%s)', status)

    path_ensemble.add_path_data(trial, status, cycle=cycle)
    if cycle % settings.get('output', {}).get('restart-file', 1) == 0:
        write_ensemble_restart(ensemble, settings)
    return accept, trial, status


def make_tis_step(ensemble, tis_settings, start_cond):
    """Perform a TIS step and generate a new path.

    The new path will be generated from the input path, either by
    performing a time-reversal move or by a shooting move. This is
    determined pseudo-randomly by drawing a random number from a
    uniform distribution using the given random generator.

    Parameters
    ----------
    ensemble : dictionary of objects
        It contains:

        * `path_ensemble`: object like :py:class:`.PathEnsemble`
          This is the path ensemble to perform the TIS step for.
        * `path`: object like :py:class:`.PathBase`
          This is the input path which will be used for generating a
          new path.
        * `system`: object like :py:class:`.System`
          System is used here since we need access to the temperature
          and to the particle list.
        * `order_function`: object like :py:class:`.OrderParameter`
          The class used for obtaining the order parameter(s).
        * `engine`: object like :py:class:`.EngineBase`
          The engine to use for propagating a path.
        * `rgen`: object like :py:class:`.RandomGenerator`
          This is the random generator that will be used.
        * `interfaces`: list of floats
          These are the interface positions on form
          [left, middle, right].

    tis_settings : dict
        This dictionary contains the settings for the TIS method. Here we
        explicitly use:

        * `freq`: float, the frequency of how often we should do time
          reversal moves.
        * `shooting_move`: string, the label of the shooting move to perform.

    start_cond : string
        The starting condition for the path. This is determined by the
        ensemble we are generating for - it is 'R'ight or 'L'eft.

    Returns
    -------
    out[0] : boolean
        True if the new path can be accepted.
    out[1] : object like :py:class:`.PathBase`
        The generated path.
    out[2] : string
        The status of the path.

    """
    interfaces = ensemble['interfaces']
    order_function = ensemble['order_function']
    path_ensemble = ensemble['path_ensemble']
    rgen = ensemble['rgen']
    last_path = path_ensemble.last_path
    # Get weights
    moves = ['tr', 'mirror', 'targetswap', 'sh']
    weights = [tis_settings.get('freq', 0),
               tis_settings.get('mirror_freq', 0),
               tis_settings.get('target_freq', 0),
               0]

    # Mirror and target swap only makes sense for ensemble 000
    if path_ensemble.ensemble_number != 0:
        weights[moves.index('mirror')] = 0
        weights[moves.index('targetswap')] = 0

    # Add sh weight (make sure it is non-negative)
    weights[-1] = max(1-sum(weights), 0)

    if sum(weights) > 1:
        msg = "move probabilities sum to more than 1. Got the following "
        msg += f"probabilities: {list(zip(moves, weights))}"
        raise ValueError(msg)

    choice = rgen.choice(moves, 1, p=weights)[0]
    if choice == 'tr':
        logger.info('Performing a time reversal move')
        accept, new_path, status = time_reversal(last_path, order_function,
                                                 interfaces, start_cond)
    elif choice == 'sh':
        logger.info('Selecting the shooting move.')
        accept, new_path, status = select_shoot(ensemble, tis_settings,
                                                start_cond)
    elif choice == 'mirror':
        logger.info("Selecting the mirror move")
        accept, new_path, status = mirror(ensemble)

    elif choice == 'targetswap':
        logger.info("Selecting the target swap move")
        accept, new_path, status = target_swap(ensemble, tis_settings)

    else:
        raise RuntimeError(f"Illegal move choice: {choice}")

    return accept, new_path, status


def select_shoot(ensemble, tis_settings, start_cond):
    """Select the shooting move to generate a new path.

    The new path will be generated from the input path, either by
    performing a normal shooting or web-throwing. This is
    determined pseudo-randomly by drawing a random number from a
    uniform distribution using the given random generator.

    Parameters
    ----------
    ensemble : dictionary of objects
        It contains:

        * `path_ensemble`: object like :py:class:`.PathEnsemble`
          This is the path ensemble to perform the TIS step for.
        * `path`: object like :py:class:`.PathBase`
          This is the input path which will be used for generating a
          new path.
        * `system`: object like :py:class:`.System`
          System is used here since we need access to the temperature
          and to the particle list.
        * `order_function`: object like :py:class:`.OrderParameter`
          The class used for obtaining the order parameter(s).
        * `engine`: object like :py:class:`.EngineBase`
          The engine to use for propagating a path.
        * `rgen`: object like :py:class:`.RandomGenerator`
          This is the random generator that will be used.
        * `interfaces`: list of floats
          These are the interface positions on form
          [left, middle, right].

    tis_settings : dict
        This dictionary contains the settings for the TIS method. Here we
        explicitly use:

        * `freq`: float, the frequency of how often we should do time
          reversal moves.
        * `shooting_move`: string, the label of the shooting move to perform.

    start_cond : string
        The starting condition for the path. This is determined by the
        ensemble we are generating for - it is 'R'ight or 'L'eft.

    Returns
    -------
    out[0] : boolean
        True if the new path can be accepted.
    out[1] : object like :py:class:`.PathBase`
        The generated path.
    out[2] : string
        The status of the path.

    """
    shooting_move = tis_settings.get('shooting_move', 'sh')

    if shooting_move == 'wt':
        logger.info('Performing a Web Throwing move')
        accept, new_path, status = web_throwing(ensemble, tis_settings)

    elif shooting_move == 'ss':
        logger.info('Performing a Stone Skipping move')
        accept, new_path, status = stone_skipping(ensemble, tis_settings,
                                                  start_cond)
    elif shooting_move == 'wf':
        logger.info('Performing a Wire Fencing move')
        accept, new_path, status = wire_fencing(ensemble, tis_settings,
                                                start_cond)
    elif shooting_move == '00':
        logger.info('Performing a null move')
        accept, new_path, status = null_move(ensemble, 0)
    else:
        logger.info('Performing a shooting move')
        accept, new_path, status = shoot(ensemble, tis_settings, start_cond)

    return accept, new_path, status


def time_reversal(path, order_function, interfaces, start_condition):
    """Perform a time-reversal move.

    Parameters
    ----------
    path : object like :py:class:`.PathBase`
        This is the input path which will be used for generating a
        new path.
    order_function : object like :py:class:`.OrderParameter`
        The class used for obtaining the order parameter(s).
    interfaces : list/tuple of floats
        These are the interface positions of the form
        ``[left, middle, right]``.
    start_condition : string
        The starting condition, 'L'eft or 'R'ight.

    Returns
    -------
    out[0] : boolean
        True if the path can be accepted.
    out[1] : object like :py:class:`.PathBase`
        The time reversed path.
    out[2] : string
        Status of the path, this is one of the strings defined in
        `.path._STATUS`.

    """
    new_path = path.reverse(order_function)
    start, _, _, _ = new_path.check_interfaces(interfaces)
    # Explicitly set how this was generated. Keep the old if loaded.
    if path.get_move() == 'ld':
        new_path.set_move('ld')
    else:
        new_path.generated = ('tr', 0, 0, 0)

    if start in set(start_condition):
        accept = True
        status = 'ACC'
    else:
        accept = False
        status = 'BWI'  # Backward trajectory end at wrong interface.
    new_path.status = status
    return accept, new_path, status


def mirror(ensemble):
    """Perform a mirror move.

    This function performs a mirror move by mirroring the order function in
    place. It then calculates the new order parameter for each phase point in
    the path. The generated path is then returned. If the path is loaded, then
    the newly generated path will also carry the generate-attribute 'ld'.
    This move is always acceptable, so the status is set to 'ACC', and the
    bool accept is set to True.

    Parameters
    ----------
    ensemble : dictionary of objects

    Returns
    -------
    out[0] : boolean
        True if the path can be accepted.
    out[1] : object like :py:class:`.PathBase`
        The generated path.
    out[2] : string
        Status of the path,

    """
    order_function = ensemble['order_function']
    # This also mirrors the ensemble['of'] inplace
    order_function.mirror()
    path = ensemble['path_ensemble'].last_path
    new_path = path.copy()
    # Calculate the new orderparameter
    for phasepoint in new_path.phasepoints:
        phasepoint.order = order_function.calculate(phasepoint)

    # This does not fix loaded paths
    if path.get_move() == 'ld':
        new_path.set_move('ld')
    else:
        new_path.generated = ('mr', 0, 0, 0)

    # This is an always accept move
    accept = True
    status = 'ACC'
    new_path.status = status
    return accept, new_path, status


def get_target_swap_choices(allowed_indices, path, ensemble):
    """Get the valid choices for the target swap move.

    Parameters
    ----------
    allowed_indices : list of integers
        The allowed indices for the target swap.
    path : object like :py:class:`.PathBase`
        The path to use for the target swap.
    ensemble : dictionary of objects
        The ensemble to use for the target swap.

    Returns
    -------
    choices : list of tuples
        The choices for the target swap.
        tuple = (index, frame)

    """
    order_function = ensemble['order_function']
    left, _, right = tuple(ensemble['interfaces'])
    choices = []
    old_idx = order_function.index
    for idx in allowed_indices:
        # swap out index on the order function
        order_function.index = idx
        for frame, phasepoint in enumerate(path.phasepoints):
            if left <= order_function.calculate(phasepoint)[0] < right:
                choices.append((idx, frame))
    order_function.index = old_idx
    return choices


def target_swap_backward(choice, choices, old_path, new_path, ensemble):
    """Perform the backward part of the target swap move.

    This function loops backward from a selected frame (the "trial point")
    until it detects a frame where the new target has a position outside a
    certain interval. The number of steps taken in this process is recorded as
    `n_n`. If there are not enough frames, it propagates backwards to
    generate more.

    Parameters
    ----------
    choice : tuple
        The choice to use for the target swap.
        tuple = (index, frame)
    choices : list of tuples
        The choices for the target swap.
        tuple = (index, frame)
    old_path : object like :py:class:`.PathBase`
        The path to use for the target swap.
    new_path : object like :py:class:`.PathBase`
        The new path to use for the target swap.
    ensemble : dictionary of objects
        The ensemble to use for the target swap.

    Returns
    -------
    n_n : integer
        The number of valid frames in the new path where a swap could
        potentially occur.
    n_o : integer
        The number of valid frames in the old path where a swap could
        potentially occur.

    """
    order_function = ensemble['order_function']
    engine = ensemble['engine']
    old_idx = order_function.index
    order_function.index = choice[0]
    n_n = 0
    first_frame = choice[1]
    n_o = 0

    # First loop back:
    for i in range(choice[1], -1, -1):
        if (choice[0], i) not in choices:
            break
        n_n += 1
        first_frame = i

        # Prevent double frames if we need to extend
        if i != 0:
            phasepoint = old_path.phasepoints[i].copy()
            phasepoint.order = order_function.calculate(phasepoint)
            new_path.append(phasepoint)

    if first_frame != 0:
        # We know we can append at least 1 more frame that would end up in a
        # state
        extra_frame = first_frame - 1
        logger.debug("Found path starting at frame %s, skipping backward MD",
                     extra_frame)
        phasepoint = old_path.phasepoints[extra_frame].copy()
        phasepoint.order = order_function.calculate(phasepoint)
        new_path.append(phasepoint)
        # Not all old_points are valid
        if extra_frame != 0:
            n_o -= extra_frame-1
            logger.debug("n_old updated from %s to %s",
                         n_o+(extra_frame-1), n_o)

    else:
        logger.debug("Need more frames, propagating backwards")
        ensemble['system'] = old_path.phasepoints[0].copy()
        engine.propagate(new_path, ensemble, reverse=True)
    order_function.idx = old_idx
    return n_n, n_o


def target_swap_forward(choice, choices, old_path, new_path, ensemble):
    """Perform the forward part of the target swap move.

    This function loops forward from a selected frame (the "trial point") until
    it detects a frame where the new target has a position outside a certain
    interval. The number of steps taken in this process is recorded as `n_n`.
    If there are not enough frames, it propagates forwards to generate more.

    Parameters
    ----------
    choice : tuple
        The choice to use for the target swap.
        tuple = (index, frame)
    choices : list of tuples
        The choices for the target swap.
        tuple = (index, frame)
    old_path : object like :py:class:`.PathBase`
        The path to use for the target swap.
    new_path : object like :py:class:`.PathBase`
        The new path to use for the target swap.
    ensemble : dictionary of objects
        The ensemble to use for the target swap.

    Returns
    -------
    n_n : integer
        The number of valid frames in the new path where a swap could
        potentially occur.
    n_o : integer
        The number of valid frames in the old path where a swap could
        potentially occur.

    """
    order_function = ensemble['order_function']
    engine = ensemble['engine']
    old_idx = order_function.index
    order_function.index = choice[0]
    n_n = 0
    last_frame = choice[1]
    n_o = 0

    for i in range(choice[1]+1, old_path.length, 1):
        if (choice[0], i) not in choices:
            break
        n_n += 1
        last_frame = i
        # Prevent last point if we need to extend
        if i != old_path.length-1:
            phasepoint = old_path.phasepoints[i].copy()
            phasepoint.order = order_function.calculate(phasepoint)
            new_path.append(phasepoint)

    if last_frame < old_path.length-1:
        # We know we can append at least one frame that will result in a valid
        # path
        extra_frame = last_frame + 1
        logger.debug("Found path ending at frame %s, skipping forward MD",
                     extra_frame)
        phasepoint = old_path.phasepoints[extra_frame].copy()
        phasepoint.order = order_function.calculate(phasepoint)
        new_path.append(phasepoint)
        # Not all old_points are valid
        # -2 because we already dismissed the last frame
        # Need to catch when we append the last frame here as we don't want to
        # add that one
        if extra_frame != old_path.length-1:
            n_o -= (old_path.length-2)-extra_frame
            logger.debug("n_old updated from %s to %s",
                         n_o+((old_path.length-2)-extra_frame), n_o)
    else:
        logger.debug("Need more frames, propagating forwards")
        ensemble['system'] = old_path.phasepoints[-1].copy()
        engine.propagate(new_path, ensemble, reverse=False)

    order_function.index = old_idx
    return n_n, n_o


def target_swap(ensemble, tis_settings):
    """Perform a target swap move.

    This function attempts a target swap move starting from a given index to
    allowed indices. It finds a trial point from which a new path is generated.
    It records Z(o->n) and Z(n->o), which are counts of valid frames where a
    target swap could potentially occur (from old to new, and from new to old,
    respectively). These are used to calculate the acceptance probability.

    Once a trial point is chosen, it loops forwards and backwards in the
    trajectory of the old path (along the new choice index) to construct a new
    path. This may include MD integration if necessary. If there is no valid
    swap, it returns with status 'TSS'.

    Parameters
    ----------
    ensemble : dictionary of objects
        It contains:

        * `path_ensemble`: object like :py:class:`.PathEnsemble`
            This is the path ensemble to perform the target swap move for.
        * `order_function`: object like :py:class:`.OrderParameter`
            The class used for obtaining the order parameter(s).

    tis_settings : dict
        This dictionary contains the settings for the TIS method. Here we
        explicitly use:

        * `freq`: float, the frequency of how often we should do time
          reversal moves.
        * `shooting_move`: string, the label of the shooting move to perform.

    Returns
    -------
    tuple : A tuple containing:
        - accept (bool): Whether the move was accepted or not.
        - new_path (object like :py:class:`.PathBase`): The new path.
        - status (str): The status of the move.
          Can be 'ACC', 'TSS', 'BTX', or 'FTX'

    Raises
    ------
    ValueError : If the current order parameter index is not part of the target
        indices.

    """
    # Grab required variables
    target_indices = tis_settings['target_indices']
    order_f = ensemble['order_function']
    old_index = order_f.index
    if old_index not in target_indices:
        raise ValueError(f"Current orderparameter index: {old_index}, is not "
                         f"part of the target indices: {target_indices}.")
    old_path = ensemble['path_ensemble'].last_path
    maxlen = tis_settings['maxlength']
    new_path = old_path.empty_path(maxlen=maxlen)
    rgen = ensemble['rgen']

    allowed_indices = [i for i in target_indices if i != old_index]

    logger.debug("Attempt a target swap starting from index %s "
                 "to the allowed indices: %s",
                 old_index, allowed_indices)

    # 2 Find trial point
    choices = get_target_swap_choices(allowed_indices=allowed_indices,
                                      path=old_path,
                                      ensemble=ensemble)
    # If there is no valid swap
    if len(choices) == 0:
        logger.info("Found no valid swapping frames")
        accept = False
        status = 'TSS'
        new_path = old_path.copy()
        new_path.generated = ('ts', 0, 0, 0)
        new_path.status = status
        return accept, new_path, status
    choice_idx = rgen.choice(len(choices), 1)[0]

    choice = choices[choice_idx]
    logger.debug("Choose %s out of %s.", choice, choices)
    logger.info("Attempting target swap from index %s to index %s",
                old_index, choice[0])

    # Record Z(o->n)
    z_o_n = len(choices)
    # Figure out n_o
    n_n = 0
    # default n_o is the old path min the 2 frames in the states
    n_o = old_path.length-2

    # First loop backwards
    diff_n_n, diff_n_o = target_swap_backward(choice=choice,
                                              choices=choices,
                                              old_path=old_path,
                                              new_path=new_path,
                                              ensemble=ensemble)
    n_n += diff_n_n
    n_o += diff_n_o

    new_path = new_path.reverse()
    # Needed because path.reverse does not clone attributes
    new_path.maxlen = maxlen
    if new_path.length >= maxlen:
        order_f.index = old_index
        logger.debug("Backward trajectory to long")
        accept = False
        status = 'BTX'
        new_path.generated = ('ts', 0, 0, 0)
        new_path.status = status
        return accept, new_path, status

    # Loop forward:
    diff_n_n, diff_n_o = target_swap_forward(choice=choice,
                                             choices=choices,
                                             old_path=old_path,
                                             new_path=new_path,
                                             ensemble=ensemble)
    n_n += diff_n_n
    n_o += diff_n_o
    if new_path.length >= maxlen:
        order_f.index = old_index
        logger.debug("Forward trajectory to long")
        accept = False
        status = 'FTX'
        new_path.generated = ('ts', 0, 0, 0)
        new_path.status = status
        return accept, new_path, status

    # Need to calculate Z(n->o)
    allowed_indices = [i for i in target_indices if i != choice[0]]
    new_choices = get_target_swap_choices(allowed_indices=allowed_indices,
                                          path=new_path,
                                          ensemble=ensemble)

    z_n_o = len(new_choices)

    logger.debug("Going to calculate acceptance\n"
                 "Old choices: %s\n"
                 "n_n: %s, Z_o_n: %s\n"
                 "New choices: %s\n"
                 "n_o: %s, Z_n_o: %s\n",
                 choices,
                 n_n, z_o_n,
                 new_choices,
                 n_o, z_n_o)

    acc = min(1, n_o/n_n*z_o_n/z_n_o)
    ran_num = rgen.rand()[0]
    logger.info("prob = %s, random = %s", acc, ran_num)
    if ran_num >= acc:
        order_f.index = old_index
        logger.debug("Failed acceptance")
        accept = False
        status = "TSA"
        new_path.generated = ('ts', 0, 0, 0)
        new_path.status = status
        return accept, new_path, status

    order_f.index = choice[0]
    logger.debug("Accepted target swap")
    accept = True
    status = "ACC"
    new_path.generated = ('ts', 0, 0, 0)
    new_path.status = status
    return accept, new_path, status


def prepare_shooting_point(path, ensemble, tis_settings):
    """Select and modify velocities for a shooting move.

    This method will randomly select a shooting point from a given
    path and modify its velocities.

    Parameters
    ----------
    path : object like :py:class:`.PathBase`
        This is the input path which will be used for generating a
        new path.
    ensemble : dictionary of objects
        It contains:

        * `order_function`: object like :py:class:`.OrderParameter`
          The class used for obtaining the order parameter(s).
        * `engine`: object like :py:class:`.EngineBase`
          The engine to use for propagating a path.
        * `rgen`: object like :py:class:`.RandomGenerator`
          This is the random generator that will be used.
        * `interfaces`: list/tuple of floats
          These are the interface positions of the form
          ``[left, middle, right]``.

    tis_settings : dict
        This contains the settings for TIS. Here, we use the
        settings which dictates how we modify the velocities.
        * `aimless`: boolean, is the shooting aimless or not?

    Returns
    -------
    out[0] : object like :py:class:`.System`
        The shooting point with modified velocities.
    out[1] : integer
        The index of the shooting point in the original path.
    out[2] : float
        The change in kinetic energy when modifying the velocities.

    """
    if path.length < 3:
        ens_name = ensemble['path_ensemble'].ensemble_name
        if path.get_move() == 'ld':
            msg = (
                f'Path of length {path.length} was loaded for ensemble '
                f'{ens_name} but contains no interior frames to shoot from. '
                'The loaded file does not provide a suitable initial path for '
                'this ensemble.\n'
                'Suggestion: add "load_and_kick = True" in the [initial-path] '
                'section of the input file so PyRETIS will kick a frame from '
                'the loaded path to generate a valid initial path.'
            )
        else:
            msg = (
                f'Path of length {path.length} was generated for ensemble '
                f'{ens_name} but contains no interior frames to shoot from. '
                'A shootable path requires at least 3 frames.\n'
                'Suggestions to resolve this:\n'
                '  - Reduce the timestep\n'
                '  - Decrease the number of subcycles\n'
                '  - Move the interfaces further apart'
            )
        raise ValueError(msg)

    shooting_point, idx = path.get_shooting_point(
        criteria=tis_settings.get('shooting_move', 'rnd'),
        interfaces=ensemble.get('interfaces'))
    engine = ensemble['engine']
    orderp = shooting_point.order
    logger.info('Shooting from order parameter/index: %f, %d', orderp[0], idx)
    # Copy the shooting point, so that we can modify velocities without
    # altering the original path:
    shooting_copy = shooting_point.copy()
    ensemble['system'] = shooting_copy
    # Modify the velocities:
    dek, _, = engine.modify_velocities(
        ensemble,
        {'sigma_v': tis_settings['sigma_v'],
         'aimless': tis_settings['aimless'],
         'zero_momentum': tis_settings['zero_momentum'],
         'rescale': tis_settings['rescale_energy']})
    orderp = engine.calculate_order(ensemble)
    shooting_copy.order = orderp
    return shooting_copy, idx, dek


def one_step_crossing(ensemble, interface):
    """Create a path of one step and check the crossing with the interface.

    This function will do a single step to try to cross the interface.
    Note that a step might involve several substeps, depending on the
    input file selection/sampling strategy.
    This task is the Achilles` Heel of Stone Skipping. If not wisely
    done, a large number of attempts to cross the interface
    in one step will be done, destroying the sampling efficiency.


    Parameters
    ----------
    ensemble : dictionary of objects
        It contains:

        * `system`: object like :py:class:`.System`
          System is used here since we need access to the temperature
          and to the particle list.
        * `order_function`: object like :py:class:`.OrderParameter`
          The class used for obtaining the order parameter(s).
        * `engine`: object like :py:class:`.EngineBase`
          The engine to use for propagating a path.
        * `rgen`: object like :py:class:`.RandomGenerator`
          This is the random generator that will be used.
    interface : float
        This is the interface position to be crossed.

    Returns
    -------
    out[0] : boolean
        True if the path can be accepted.
    out[1] : object like :py:class:`.PathBase`
        Returns the generated path.

    """
    # The trial path we need to generate. Note, 1 step = 2 points
    trial_path = Path(rgen=ensemble['rgen'], maxlen=2)
    interfaces = [-float("inf"), float("inf"), float("inf")]
    sub_ensemble = {'interfaces': interfaces,
                    'system': ensemble['system'],
                    'order_function': ensemble['order_function']}
    ensemble['engine'].propagate(trial_path, sub_ensemble)
    if crossing_counter(trial_path, interface) == 0:
        return False, trial_path

    return True, trial_path


def check_kick(shooting_point, interfaces, trial_path, rgen, dek,
               tis_settings):
    """Check the modification of the shooting point.

    After generating velocities for a shooting point, we
    do some additional checking to see if the shooting point is
    acceptable.

    Parameters
    ----------
    shooting_point : object like :py:class:`.System`
        The shooting point with modified velocities.
    interfaces : list of floats
        The interfaces used for TIS, in the format
        ``[left, middle, right]``.
    trial_path : object like :py:class:`.PathBase`
        The path we are currently generating.
    rgen : object like :py:class:`.RandomGenerator`
        This is the random generator that will be used to check if
        we accept the shooting point based on the change in kinetic
        energy.
    dek : float
        The change in kinetic energy when modifying the velocities.
    tis_settings : dict
        This contains the settings for TIS.

    Returns
    -------
    out : boolean
        True if the kick was OK, False otherwise.

    """
    # 1) Check if the kick was too violent:
    left, _, right = interfaces
    if 'exp' in tis_settings.get('shooting_move', {}):
        return True
    if not left <= _order_value(shooting_point.order) < right:
        # Shooting point was velocity dependent and was kicked outside
        # of boundaries when modifying velocities.
        trial_path.append(shooting_point)
        trial_path.status = 'KOB'
        return False
    # 2) If the kick is not aimless, we check if we reject it or not:
    if not tis_settings['aimless']:
        accept_kick = metropolis_accept_reject(rgen, shooting_point, dek)
        # If one wish to implement a bias call, this can be done here.
        if not accept_kick:
            trial_path.append(shooting_point)
            trial_path.status = 'MCR'  # Momenta Change Rejection.
            return False
    return True


def shoot_backwards(path_back, trial_path, ensemble,
                    tis_settings, start_cond):
    """Shoot in the backward time direction.

    Parameters
    ----------
    path_back : object like :py:class:`.PathBase`
        The path we will fill with phase points from the propagation.
    trial_path : object like :py:class:`.PathBase`
        The current trial path generated by the shooting.
    ensemble : dictionary of objects
        It contains:

        * `order_function`: object like :py:class:`.OrderParameter`
          The class used for obtaining the order parameter(s).
        * `engine`: object like :py:class:`.EngineBase`
          The engine to use for propagating a path.
        * `interfaces`: list/tuple of floats
          These are the interface positions of the form
          ``[left, middle, right]``.
        * `system`: object like :py:class:`.System`
           The system that originates the path.

    tis_settings : dict
        This contains the settings for TIS.
    start_cond : string
        The starting condition for the current ensemble, 'L'eft or
        'R'ight.

    Returns
    -------
    out : boolean
        True if the backward path was generated successfully, False
        otherwise.

    """
    logger.debug('Propagating backwards for the shooting move.')
    path_back.time_origin = trial_path.time_origin
    engine = ensemble['engine']
    success_back, _ = engine.propagate(path_back, ensemble, reverse=True)
    if not success_back:
        # Something went wrong, most probably the path length was exceeded.
        trial_path.status = 'BTL'  # BTL = backward trajectory too long.
        # Add the failed path to trial path for analysis:
        trial_path += path_back
        if path_back.length >= tis_settings['maxlength'] - 1:
            # BTX is backward trajectory longer than maximum memory.
            trial_path.status = 'BTX'
        return False
    # Backward seems OK so far, check if the ending point is correct:
    left, _, right = ensemble['interfaces']
    if path_back.get_end_point(left, right) not in set(start_cond):
        # Nope, backward trajectory end at wrong interface.
        trial_path += path_back  # Store path for analysis.
        trial_path.status = 'BWI'
        return False
    return True


def shoot(ensemble, tis_settings, start_cond, shooting_point=None):
    """Perform a shooting-move.

    This function will perform the shooting move from a randomly
    selected time-slice.

    Parameters
    ----------
    ensemble: dict
        It contains:

        * `path_ensemble`: object like :py:class:`.PathEnsemble`
          This is the path ensemble to perform the TIS step for.
        * `system`: object like :py:class:`.System`
          System is used here since we need access to the temperature
          and to the particle list.
        * `order_function`: object like :py:class:`.OrderParameter`
          The class used for obtaining the order parameter(s).
        * `engine`: object like :py:class:`.EngineBase`
          The engine to use for propagating a path.
        * `rgen`: object like :py:class:`.RandomGenerator`
          This is the random generator that will be used.
        * `interfaces`: list of floats
          These are the interface positions on form
          [left, middle, right].

    tis_settings : dict
        This contains the settings for TIS:

        * `aimless`: boolean, is the shooting aimless or not?
        * `allowmaxlength`: boolean, should paths be allowed to reach
          maximum length?
        * `maxlength`: integer, maximum allowed length of paths.

    start_cond : string or tuple of strings
        The starting condition for the current ensemble, 'L'eft or
        'R'ight. ('L', 'R'), the default option, implies no directional
        difference.
    shooting_point: object like :py:class:`.System`, optional
        If given, it is the shooting point from which the path is generated.

    Returns
    -------
    out[0] : boolean
        True if the path can be accepted.
    out[1] : object like :py:class:`.PathBase`
        Returns the generated path.
    out[2] : string
        Status of the path, this is one of the strings defined in
        :py:const:`.path._STATUS`.

    """
    path_ensemble = ensemble['path_ensemble']
    path = path_ensemble.last_path
    interfaces = ensemble['interfaces']
    trial_path = path.empty_path()  # The trial path we will generate.
    if shooting_point is None:
        shooting_point, idx, dek = prepare_shooting_point(
            path, ensemble, tis_settings
        )
        kick = check_kick(shooting_point, interfaces, trial_path, path.rgen,
                          dek, tis_settings)
    else:
        kick = True
        idx = getattr(shooting_point, 'idx', 0)

    # Store info about this point, just in case we have to return
    # before completing a full new path:
    trial_path.generated = ('sh', _order_value(shooting_point.order), idx, 0)
    trial_path.time_origin = path.time_origin + idx
    # We now check if the kick was OK or not:
    if not kick:
        return False, trial_path, trial_path.status
    # OK: kick was either aimless or it was accepted by Metropolis
    # we should now generate trajectories, but first check how long
    # it should be (if the path comes from a load, it is assumed to not
    # respect the detail balance anyway):
    if path.get_move() == 'ld' or tis_settings['allowmaxlength']:
        maxlen = tis_settings['maxlength']
    else:
        maxlen = min(int((path.length - 2) / path.rgen.rand()[0]) + 2,
                     tis_settings['maxlength'])
    # Since the forward path must be at least one step, the maximum
    # length for the backward path is maxlen-1.
    # Generate the backward path:
    path_back = path.empty_path(maxlen=maxlen - 1)
    # Set ensemble state to the selected shooting point:
    ensemble['system'] = shooting_point.copy()
    if not shoot_backwards(path_back, trial_path, ensemble,
                           tis_settings, start_cond):
        return False, trial_path, trial_path.status

    # Everything seems fine, now propagate forward.
    # Note that the length of the forward path is adjusted to
    # account for the fact that it shares a point with the backward
    # path (i.e. the shooting point). The duplicate point is just
    # counted once when the paths are merged by the method
    # `paste_paths` by setting `overlap=True` (which indicates that
    # the forward and backward paths share a point).
    path_forw = path.empty_path(maxlen=maxlen - path_back.length + 1)
    logger.debug('Propagating forwards for shooting move...')
    # Set ensemble state to the selected shooting point:
    # change the system state.
    ensemble['system'] = shooting_point.copy()
    success_forw, _ = ensemble['engine'].propagate(path_forw, ensemble,
                                                   reverse=False)
    path_forw.time_origin = trial_path.time_origin
    # Now, the forward propagation could have failed by exceeding the
    # maximum length for the forward path. However, it could also fail
    # when we paste together so that the length is larger than the
    # allowed maximum. We paste first and ask later:
    trial_path = paste_paths(path_back, path_forw, overlap=True,
                             maxlen=tis_settings['maxlength'])

    # Also update information about the shooting:
    trial_path.generated = ('sh', _order_value(shooting_point.order), idx,
                            path_back.length - 1)
    trial_path.weight = 1.
    if not success_forw:
        trial_path.status = 'FTL'
        # If we reached this point, the backward path was successful,
        # but the forward was not. For the case where the forward was
        # also successful, the length of the trial path cannot exceed
        # the maximum length given in the TIS settings. Thus we only
        # need to check this here, i.e. when given that the backward
        # was successful and the forward not:
        if trial_path.length == tis_settings['maxlength']:
            trial_path.status = 'FTX'  # exceeds "memory".

    # Deal with the rejections for path properties.
    # Make sure we did not hit the left interface on {0-} if this is disallowed
    # Cases:
    # Permeability=true: starting/ending at L is allowed in [0-]
    # RETIS: [0-] is the only ensemble that allows paths starting R
    # (RE)PPTIS: starting/ending at L and R is allowed for body ensembles.
    #            exception: [0^{+-}'] is not allowed to start AND end at R.
    #                       As this ensemble has must_cross_M=True, and L=M,
    #                       This is automatically taken care of.
    elif ('L' not in set(path_ensemble.start_condition) and
            'L' in trial_path.check_interfaces(interfaces)[:2]):
        trial_path.status = '0-L'

    # Last check - Did we cross the middle interface?
    # if we are using PPTIS or REPPTIS
    elif getattr(path_ensemble, 'must_cross_M', False) and \
            (not trial_path.check_interfaces(interfaces)[-1][1] or
             not (trial_path.check_interfaces(interfaces)[-1][0] or
                  trial_path.check_interfaces(interfaces)[-1][2])):
        # not for [0-]
        # detect if: minorder < middle interf <= maxorder
        # or
        # if we do cross middle interface, the
        # path still has to come from left or right
        # so we need cross[0] or cross[2] to be true
        trial_path.status = 'NCR'

    # if we are not using PPTIS or REPPTIS
    # Special treatment if paths can start everywhere
    elif set(('R', 'L')) != set(path_ensemble.start_condition) and \
            not trial_path.check_interfaces(interfaces)[-1][1]:
        # No, we did not cross the middle interface:
        trial_path.status = 'NCR'
    else:  # If nothing went wrong, then...
        trial_path.status = 'ACC'
        return True, trial_path, trial_path.status

    return False, trial_path, trial_path.status


def wire_fencing(ensemble, tis_settings, start_cond):
    """Perform a wire_fencing move.

    This function will perform the non-famous wire fencing move
    from an initial path.

    Parameters
    ----------
    ensemble: dict
        It contains:

        * `path_ensemble`: object like :py:class:`.PathEnsemble`
          This is the path ensemble to perform the TIS step for.
        * `system`: object like :py:class:`.System`
          System is used here since we need access to the temperature
          and to the particle list.
        * `order_function`: object like :py:class:`.OrderParameter`
          The class used for obtaining the order parameter(s).
        * `engine`: object like :py:class:`.EngineBase`
          The engine to use for propagating a path.
        * `rgen`: object like :py:class:`.RandomGenerator`
          This is the random generator that will be used.
        * `interfaces`: list of floats
          These are the interface positions on form
          [left, middle, right].

    tis_settings : dict
        This contains the settings for TIS. Keys used here:

        * `aimless`: boolean, is the shooting aimless or not?
        * `allowmaxlength`: boolean, should paths be allowed to reach
          maximum length?
        * `maxlength`: integer, maximum allowed length of paths.
        * `high_accept`: boolean, the option for High Acceptance WF.

    start_cond : string
        The starting condition for the current ensemble, 'L'eft or
        'R'ight.

    Returns
    -------
    out[0] : boolean
        True if the path can be accepted.
    out[1] : object like :py:class:`.PathBase`
        Returns the generated path.
    out[2] : string
        Status of the path, this is one of the strings defined in
        :py:const:`.path._STATUS`.

    """
    trial_path = ensemble['path_ensemble'].last_path
    engine = ensemble['engine']
    intf_cap = tis_settings.get('interface_cap', ensemble['interfaces'][2])
    # Wire fencing needs middle != cap interface. For the [0^-] ensemble
    # middle == right == lambda_0, making WF inapplicable.
    if ensemble['interfaces'][1] == intf_cap:
        logger.info('Wire fencing not applicable (middle == cap '
                    'interface = %s). Falling back to standard shooting.',
                    ensemble['interfaces'][1])
        return shoot(ensemble, tis_settings, start_cond)

    wf_int = [ensemble['interfaces'][1], ensemble['interfaces'][1], intf_cap]
    n_frames, new_segment = wirefence_weight_and_pick(trial_path, wf_int[0],
                                                      wf_int[2],
                                                      return_seg=True)

    if n_frames == 0:
        logger.warning('Wire fencing move not usable. N frames of Path = 0 '
                       'between interfaces %s and %s.', wf_int[0], wf_int[-1])
        return False, trial_path, 'NSG'
    sub_ens = {'interfaces': wf_int, 'engine': engine,
               'rgen': ensemble['rgen'],
               'order_function': ensemble['order_function'],
               'path_ensemble': ensemble['path_ensemble']}
    sub_settings = tis_settings.copy()
    sub_settings['allowmaxlength'] = True
    succ_seg = 0

    for i in range(tis_settings['n_jumps']):
        logger.debug('Trying a new web with Wire Fencing, jump %i', i)
        # Select the shooting point:
        sh_pt, _, _ = prepare_shooting_point(new_segment, sub_ens,
                                             sub_settings)
        engine.dump_phasepoint(sh_pt, str(counter()) + '_wf_shoot')

        success, trial_seg, status = shoot(sub_ens,
                                           sub_settings,
                                           ('L', 'R'),
                                           sh_pt)
        start, end, _, _ = trial_seg.check_interfaces(wf_int)
        logger.info('Jump %s, len %s, status %s, intf: %s %s',
                    i, trial_seg.length, status, start, end)
        if not success:
            # This handles R to R (start_cond = L) paths. Counter + 1, no ups.
            logger.debug('Wire Fencing Fail.')
        else:
            logger.debug('Acceptable Wire Fence link.')
            succ_seg += 1
            new_segment = trial_seg.copy()

    if succ_seg == 0:
        # No usable segments were generated.
        trial_path.status = 'NSG'
        success = False
    else:
        success, trial_path, _ = extender(new_segment, ensemble,
                                          tis_settings, start_cond)
    if success:
        success, trial_path = ss_wt_wf_acceptance(trial_path, ensemble,
                                                  tis_settings,
                                                  start_cond)

    trial_path.generated = ('wf', _order_value(sh_pt.order), succ_seg,
                            trial_path.length)

    logger.debug('WF move %s', trial_path.status)
    if not success:
        return False, trial_path, trial_path.status

    # This might get triggered when accepting 0-L paths.
    left, _, right = ensemble['interfaces']
    if start_cond != trial_path.get_start_point(left, right):
        raise RuntimeError('WF: Path has an implausible start.')

    trial_path.status = 'ACC'
    return True, trial_path, trial_path.status


def ss_wt_wf_acceptance(trial_path, ensemble, tis_settings,
                        start_cond='L'):
    """Weights, possibly reverses and accept/rejects generated SS/WT/WF paths.

    Parameters
    ----------
    trial_path : object like :py:class:`.PathBase`
        This is the new path that will obtain weights, and might be reversed
        and accepted.
    ensemble : dict
        It contains:

        * `path_ensemble`: object like :py:class:`.PathEnsemble`
          This is the path ensemble to perform the TIS step for.
        * `order_function`: object like :py:class:`.OrderParameter`
          The class used for obtaining the order parameter(s).
        * `rgen`: object like :py:class:`.RandomGenerator`
          This is the random generator that will be used.
        * `interfaces` : list/tuple of floats
          These are the interface positions of the form
          ``[left, middle, right]``.

    tis_settings : dict
        This contains the settings for TIS. KEys used here;

        * `high_accept` : boolean, the option for High Acceptance SS/WF.
        * `shooting_move` : string, the label of the shooting move to perform.

    start_cond : string, optional
        The starting condition for the current ensemble, 'L'eft or 'R'ight.

    Returns
    -------
    out[0] : boolean
        True if the path can be accepted.
    out[1] : object like :py:class:`.PathBase`
        Returns the weighed and possibly reversed path.

    """
    intf = list(ensemble['interfaces'])
    move = tis_settings['shooting_move']

    if move == 'wt' or not tis_settings.get('high_accept', False):
        trial_path.weight = 1.
    else:
        if move == 'wf':
            intf[2] = tis_settings.get('interface_cap', intf[2])
        trial_path.weight = compute_weight(trial_path, intf, move)
        if start_cond != trial_path.get_start_point(intf[0], intf[2]):
            trial_path = trial_path.reverse(ensemble['order_function'])

    success = ss_wt_wf_metropolis_acc(trial_path, ensemble, tis_settings,
                                      start_cond)

    return success, trial_path


def stone_skipping(ensemble, tis_settings, start_cond):
    """Perform a stone_skipping move.

    This function will perform the famous stone skipping move
    from an initial path.

    Parameters
    ----------
    ensemble: dict
        It contains:

        * `path_ensemble`: object like :py:class:`.PathEnsemble`
          This is the path ensemble to perform the TIS step for.
        * `path`: object like :py:class:`.PathBase`
          This is the input path which will be used for generating a
          new path.
        * `system`: object like :py:class:`.System`
          System is used here since we need access to the temperature
          and to the particle list.
        * `order_function`: object like :py:class:`.OrderParameter`
          The class used for obtaining the order parameter(s).
        * `engine`: object like :py:class:`.EngineBase`
          The engine to use for propagating a path.
        * `rgen`: object like :py:class:`.RandomGenerator`
          This is the random generator that will be used.
        * `interfaces`: list of floats
          These are the interface positions on form
          [left, middle, right].

    tis_settings : dict
        This contains the settings for TIS. Keys used here:

        * `aimless`: boolean, is the shooting aimless or not?
        * `allowmaxlength`: boolean, should paths be allowed to reach
          maximum length?
        * `maxlength`: integer, maximum allowed length of paths.
        * `high_accept`: boolean, the option for High Acceptance SS.

    start_cond : string
        The starting condition for the current ensemble, 'L'eft or
        'R'ight.

    Returns
    -------
    out[0] : boolean
        True if the path can be accepted.
    out[1] : object like :py:class:`.PathBase`
        Returns the generated path.
    out[2] : string
        Status of the path, this is one of the strings defined in
        :py:const:`.path._STATUS`.

    """
    path_old = ensemble['path_ensemble'].last_path
    intf = ensemble['interfaces']
    ph_pt1, ph_pt2 = crossing_finder(path_old, intf[1])
    if ph_pt1 is None and ph_pt2 is None:
        return False, path_old, 'NCR'
    sub_ens = {'interfaces': [intf[1], intf[1], intf[2]],
               'order_function': ensemble['order_function']}
    osc_try = 0  # One step crossing attempt counter
    success = False
    for i in range(tis_settings['n_jumps']):
        logger.debug('Trying a new stone skipping move, jump %i', i)
        # Here we choose between the two
        # possible shooting points that describe a crossing.
        sh_pt = ph_pt1 if ensemble['rgen'].rand() >= 0.5 else ph_pt2
        ensemble['engine'].dump_phasepoint(sh_pt, str(counter()) + '_ss_shoot')
        # To continue, we must be sure that the new path
        # CROSSES the interface in ONLY ONE step.
        # Generate paths until it succeed. That is
        # what makes this version of the SS move useless for large systems.
        for j in range(tis_settings['maxlength']):
            # This function can become actually fun to work on.
            # e.g. have a 50% chance to give random v for each particle
            # Modify the velocities:
            # todo modify_v could just use system directly
            logger.debug("jump %s, try %s, start: %s", i, j,
                         _order_value(sh_pt.order))
            ensemble['system'] = sh_pt.copy()
            ensemble['engine'].modify_velocities(ensemble, tis_settings)
            # A path of two frames is going to be generated.
            success, path = one_step_crossing(ensemble, intf[1])
            osc_try += 1
            if osc_try > 5 * tis_settings['n_jumps'] and \
                    path_old.get_move() in {'ld', 'ki', 'is'}:
                logger.info('Performing a shooting move before the use of ss')
                success, trial_path, status = shoot(ensemble,
                                                    tis_settings,
                                                    start_cond, sh_pt)
                trial_path.set_move('is')
                return success, trial_path, status

            if success:
                break
        else:  # In case we reached maxlength in jumps attempts.
            success = False
            path.status = 'NSS'
            trial_path = path
            break

        # Depending on the shooting point (before or after the interface),
        # a backward path or a continuation has to be generated.
        new_segment = path.empty_path(maxlen=tis_settings['maxlength'] - 1)
        if path.get_end_point(intf[1], intf[2]) == start_cond:
            path = path.reverse(ensemble['order_function'])
        sub_ens['system'] = path.phasepoints[1].copy()
        success, _ = ensemble['engine'].propagate(new_segment, sub_ens)
        new_segment.phasepoints.insert(0, path.phasepoints[0].copy())

        if not success:
            new_segment.status = 'XSS'
            trial_path = new_segment
            break

        ph_pt1, ph_pt2 = crossing_finder(new_segment, intf[1], last_frame=True)

    logger.debug('SS web: %s, one step crossing tries: %s', success, osc_try)

    if success:
        if ensemble['rgen'].rand() < 0.5:
            new_segment = new_segment.reverse(ensemble['order_function'])
        success, trial_path, _ = extender(new_segment, ensemble,
                                          tis_settings, start_cond)

    if success:
        success, trial_path = ss_wt_wf_acceptance(trial_path, ensemble,
                                                  tis_settings,
                                                  start_cond)

    trial_path.generated = ('ss', _order_value(sh_pt.order), osc_try,
                            trial_path.length)

    logger.debug('SS move: %s', trial_path.status)
    if not success:
        return False, trial_path, trial_path.status

    # This might get triggered when accepting 0-L paths.
    if start_cond != trial_path.get_start_point(intf[0], intf[2]):
        raise RuntimeError('SS: Path has an implausible start.')

    trial_path.status = 'ACC'

    return True, trial_path, trial_path.status


def ss_wt_wf_metropolis_acc(path_new, ensemble,
                            tis_settings, start_cond='L'):
    """Accept or reject the path_new.

    Super detailed balance rule is used in the original version
    and in the High Acceptance one for SS and WF.

    In the regular version, P acc = min (1, Cold/Cnew), where
    for Stone Skipping C is crossing, for Web Throwing C is segment,
    for Wire Fencing C is number of phasepoint between ensemble and
    right interface.

    In the High Acceptance version, P acc = 1 for SS and Wf. It
    also allows to accept paths that go from B to A, by reversing them.
    NB. To respect super detailed balance, the weights have to be changed
    accordingly. This is done elsewhere.

    Parameters
    ----------
    path_new : object like :py:class:`.PathBase`
        This is the new path that might get accepted.
    ensemble : dict
        It contains:

        * `path_ensemble`: object like :py:class:`.PathEnsemble`
          This is the path ensemble to perform the TIS step for.
        * `order_function`: object like :py:class:`.OrderParameter`
          The class used for obtaining the order parameter(s).
        * `rgen`: object like :py:class:`.RandomGenerator`
          This is the random generator that will be used.
        * `interfaces` : list/tuple of floats
          These are the interface positions of the form
          ``[left, middle, right]``.

    tis_settings : dict
        This contains the settings for TIS. Keys used here:

        * `high_accept` : boolean, the option for High Acceptance SS/WF.

    start_cond : string, optional
        The starting condition for the current ensemble, 'L'eft or
        'R'ight.

    Returns
    -------
    out[0] : boolean
        True if the path can be accepted.

    """
    interfaces = ensemble['interfaces']
    path_old = ensemble['path_ensemble'].last_path
    if tis_settings.get('shooting_move') == 'wt':
        sour_int = tis_settings['interface_sour']
        cr_old = segments_counter(path_old, sour_int, interfaces[1])
        cr_new = segments_counter(path_new, sour_int, interfaces[1])
        if ensemble['rgen'].rand() >= min(1.0, cr_old / cr_new):
            path_new.status = 'WTA'
            return False

    else:
        if not tis_settings.get('high_accept', False):
            if tis_settings.get('shooting_move') == 'ss':
                cr_old = crossing_counter(path_old, interfaces[1])
                cr_new = crossing_counter(path_new, interfaces[1])
                if ensemble['rgen'].rand() >= min(1.0, cr_old / cr_new):
                    path_new.status = 'SSA'
                    return False
            elif tis_settings.get('shooting_move') == 'wf':
                wf_cap = tis_settings.get('interface_cap', interfaces[2])
                cr_old, _ = wirefence_weight_and_pick(path_old, interfaces[1],
                                                      wf_cap)
                cr_new, _ = wirefence_weight_and_pick(path_new, interfaces[1],
                                                      wf_cap)
                if ensemble['rgen'].rand() >= min(1.0, cr_old / cr_new):
                    path_new.status = 'WFA'
                    return False

    if start_cond != path_new.get_start_point(interfaces[0], interfaces[2]):
        path_new.status = 'BWI'
        return False
    path_new.status = 'ACC'
    return True


def web_throwing(ensemble, tis_set, start_cond='L'):
    """Perform a web_throwing move.

    This function performs the great web throwing move from an initial path.

    Parameters
    ----------
    ensemble : dict
        It contains:

        * `path_ensemble`: object like :py:class:`.PathEnsemble`
          This is the path ensemble to perform the TIS step for.
        * `system`: object like :py:class:`.System`
          System is used here since we need access to the temperature
          and to the particle list.
        * `order_function`: object like :py:class:`.OrderParameter`
          The class used for obtaining the order parameter(s).
        * `engine`: object like :py:class:`.EngineBase`
          The engine to use for propagating a path.
        * `rgen`: object like :py:class:`.RandomGenerator`
          This is the random generator that will be used.
        * `interfaces`: list of floats
          These are the interface positions on form
          [left, middle, right].

    tis_set : dict
        This contains the settings for TIS. Keys used here:

        * `aimless`: boolean, is the shooting aimless or not?
        * `allowmaxlength`: boolean, should paths be allowed to reach
          maximum length?
        * `maxlength`: integer, maximum allowed length of paths.

    start_cond : string, optional
        The starting condition for the current ensemble, 'L'eft or
        'R'ight.


    Returns
    -------
    out[0] : boolean
        True if the path can be accepted.
    out[1] : object like :py:class:`.PathBase`
        Returns the generated path.
    out[2] : string
        Status of the path, this is one of the strings defined in
        :py:const:`.path._STATUS`.

    """
    path_old = ensemble['path_ensemble'].last_path
    interfaces = ensemble['interfaces']
    sour = tis_set['interface_sour']
    if not interfaces[0] < sour <= interfaces[1]:
        raise ValueError('SOUR interface is not correctly positioned')

    ccnt = segments_counter(path_old, sour, interfaces[1])
    if ccnt == 0:
        return False, path_old, 'NSG'

    seg_i = int(ensemble['rgen'].rand()[0] * ccnt)
    wt_int = [sour, sour, ensemble['interfaces'][1]]
    source_seg = select_and_trim_a_segment(path_old, sour, wt_int[2], seg_i)
    sub_ens = {'interfaces': wt_int,
               'order_function': ensemble['order_function']}

    shoots, save_acc = [0], 0

    key = ensemble['rgen'].rand() >= 0.5  # Start from a random side
    for _ in range(tis_set['n_jumps']):
        if ensemble['rgen'].rand() >= 0.5:
            shoots[-1] += 1  # One more on the Same side
        else:
            shoots.append(1)  # A move in the other side

    for n_virtual in shoots:
        key = not key  # Change side, key controls also path reverse
        for _ in range(n_virtual):
            if key:
                pre_shooting_point = source_seg.phasepoints[-1]
                shooting_point = source_seg.phasepoints[-2]
            else:
                pre_shooting_point = source_seg.phasepoints[0]
                shooting_point = source_seg.phasepoints[1]

            prefix = str(counter())
            ensemble['engine'].dump_phasepoint(pre_shooting_point,
                                               prefix + '_wt_pre_shoot')
            ensemble['engine'].dump_phasepoint(shooting_point,
                                               prefix + '_wt_shoot')

            new_seg = path_old.empty_path(maxlen=tis_set['maxlength'])
            new_seg.append(pre_shooting_point)
            logger.debug('Trying a new web')
            sub_ens['system'] = shooting_point.copy()
            ensemble['engine'].propagate(new_seg, sub_ens, reverse=key)
            start = new_seg.get_start_point(wt_int[0], wt_int[-1])
            end = new_seg.get_end_point(wt_int[0], wt_int[-1])
            logger.debug('WT web starts %s, ends %s, reverse %s',
                         start, end, key)
            if segments_counter(new_seg, sour, wt_int[2], reverse=key) == 1:
                logger.debug('Web successful')
                source_seg = new_seg.reverse(ensemble['order_function'],
                                             rev_v=False) if key else new_seg
                source_seg.status = 'ACC'
                save_acc += 1
                break

    logger.debug('WT segments accepted: %s', save_acc)

    accept, trial_path, _ = extender(source_seg, ensemble, tis_set, start_cond)

    trial_path.generated = ('wt',
                            _order_value(source_seg.phasepoints[1].order),
                            save_acc, trial_path.length)
    # Also Check that we did not get a B to A or a B to B path.
    if accept:
        accept, trial_path = ss_wt_wf_acceptance(trial_path, ensemble, tis_set)
    logger.debug('WT move: %s', trial_path.status)

    # Set the path flags
    if accept and path_old.get_move() == 'ld' and save_acc == 0:
        trial_path.set_move('ld')

    return accept, trial_path, trial_path.status


def extender(source_seg, ensemble, tis_set, start_cond=('R', 'L')):
    """Extend a path to the given interfaces.

    This function will perform the web throwing move from an initial path.

    Parameters
    ----------
    source_seg : object like :py:class:`.PathBase`
        This is the input path which will be prolonged.
    ensemble : dictionary of objects
        It contains:

        * `order_function`: object like :py:class:`.OrderParameter`
          The class used for obtaining the order parameter(s).
        * `engine`: object like :py:class:`.EngineBase`
          The engine to use for propagating a path.
        * `system`: object like :py:class:`.System`
          System is used here since we need access to the temperature
          and to the particle list.
        * `interfaces`: list of floats
          These are the interface positions on form
          [left, middle, right].

    tis_set : dict
        This contains the settings for TIS. Keys used here:

        * `aimless`: boolean, is the shooting aimless or not?
        * `allowmaxlength`: boolean, should paths be allowed to reach
          maximum length?
        * `maxlength`: integer, maximum allowed length of paths.
        * `high_accept`: boolean, the option for High Acceptance WF and SS.

    start_cond : string or tuple of strings, optional
        The starting condition for the current ensemble, 'L'eft or
        'R'ight. ('L', 'R'), the default option, implies no directional
        difference.


    Returns
    -------
    out[0] : boolean
        True if the path can be accepted.
    out[1] : object like :py:class:`.PathBase`
        Returns the generated path.
    out[2] : string
        Status of the path, this is one of the strings defined in
        :py:const:`.path._STATUS`.

    """
    interfaces = ensemble['interfaces']
    ensemble['system'] = source_seg.phasepoints[0].copy()

    # Extender
    if interfaces[0] <= _order_value(ensemble['system'].order) < \
            interfaces[-1]:
        back_segment = source_seg.empty_path(maxlen=tis_set['maxlength'])
        logger.debug('Trying to extend backwards')
        source_seg_copy = source_seg.copy()

        if not shoot_backwards(back_segment, source_seg_copy, ensemble,
                               tis_set, start_cond):
            if not tis_set.get('high_accept', False):
                return False, source_seg_copy, source_seg_copy.status

        trial_path = paste_paths(back_segment, source_seg, overlap=True,
                                 maxlen=tis_set['maxlength'])
    else:
        trial_path = source_seg.copy()

    ensemble['system'] = trial_path.phasepoints[-1].copy()
    if interfaces[0] <= _order_value(ensemble['system'].order) < \
            interfaces[-1]:
        forth_segment = source_seg.empty_path(maxlen=tis_set['maxlength'])
        ensemble['engine'].propagate(forth_segment, ensemble)

        trial_path.phasepoints = trial_path.phasepoints[:-1] + \
            forth_segment.phasepoints

    if trial_path.length >= tis_set['maxlength']:
        trial_path.status = 'FTX'  # exceeds "memory".
        return False, trial_path, trial_path.status
    trial_path.status = 'ACC'
    return True, trial_path, trial_path.status
