# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This file contains functions used for initiation of paths.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

generate_initial_path_kick (:py:func:`.generate_initial_path_kick`)
    Function for generating an initial path by repeatedly kicking a
    phase point.

initiate_kick (:py:func:`.initiate_kick`)
    Helper method, selects either :py:func:`.initiate_kicki` or
    :py:func:`.initiate_kick_max`.

initiate_kicki (:py:func:`.initiate_kicki`)
    A method for initiating a path ensemble by repeatedly modifying
    velocities to find the crossing with the interfaces.

initiate_kick_max (:py:func:`.initiate_kick_max`)
    A method similar to py:meth:`.initiate_kick`. Here, if possible,
    we will use points from the previous paths, closest to the target
    interface.

initiate_path_ensemble_kick (:py:func:`.initiate_path_ensemble_kick`)
    Method to initiate a single path ensemble.

fix_path_by_tis (:py:func:`.fix_path_by_tis`)
    Method to fix an initial path that starts and ends at the wrong
    interface via TIS moves.

"""
import copy
import logging
from pyretis.core.common import compute_weight
from pyretis.core.path import paste_paths, Path
from pyretis.core.tis import shoot, time_reversal
from pyretis.inout.screen import REFERENCE
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = [
    'generate_initial_path_kick',
    'initiate_kick',
    'initiate_kicki',
    'initiate_kick_max',
    'initiate_path_ensemble_kick',
    'fix_path_by_tis'
]


def _order_value(order):
    """Return the first order parameter component (scalar or vector)."""
    try:
        return order[0]
    except (TypeError, IndexError):
        return order


def initiate_kick(simulation, settings, cycle):
    """Run the initiation method for several ensembles.

    Please see documentation for :py:func:`.initiate_path_ensemble_kick`.

    """
    init_settings = settings['initial-path']
    start = init_settings.get('kick-from', 'initial').lower()
    if start == 'previous':
        logger.progress('Kick-initiate using previous configuration')
        return initiate_kick_max(simulation, settings, cycle)
    if start == 'initial':
        logger.progress('Kick-initiate using initial configuration')
        return initiate_kicki(simulation, settings, cycle)
    errtxt = f'Unknown argument {start} for kick-from'
    logger.error(errtxt)
    raise ValueError(errtxt)


def initiate_kicki(simulation, settings, cycle):
    """Run the initialization for several ensembles.

    Please see documentation for :py:func:`.initiate_path_ensemble_kick`.

    """
    for ensemble, set_ens in zip(simulation.ensembles, settings['ensemble']):
        ens_name = ensemble['path_ensemble'].ensemble_name
        logger.log(REFERENCE, '\n--- Initiating path ensemble: %s ---',
                   ens_name)
        accept, initial_path, status = initiate_path_ensemble_kick(
            ensemble,
            set_ens['tis'],
            cycle)
        yield accept, initial_path, status, ensemble['path_ensemble']


def initiate_kick_max(simulation, settings, cycle):
    """Run the initiation for several ensembles.

    This method is similar to :py:func:`.initiate_kick`, but here we
    update the initial point for an ensemble to use that of the previous
    path (if this exists).

    Please see documentation for :py:func:`.initiate_path_ensemble_kick`.

    """
    last_paths = []
    for ensemble, set_ens in zip(simulation.ensembles, settings['ensemble']):
        path_ensemble = ensemble['path_ensemble']
        ens_name = path_ensemble.ensemble_name
        logger.log(REFERENCE, '\n--- Initiating path ensemble: %s ---',
                   ens_name)
        if last_paths:
            # Look for the phase point closest to the middle interface.
            middle = path_ensemble.interfaces[1]
            current = None
            min_dist = float('inf')
            for last_path in last_paths:
                for phase_point in last_path.phasepoints:
                    dist = middle - _order_value(phase_point.order)
                    if 0 <= dist <= min_dist:
                        min_dist = dist
                        current = phase_point.copy()
            if current is not None:
                logger.info('Initial state set to: %s', current)
                ensemble['system'] = current
        accept, initial_path, status = initiate_path_ensemble_kick(
            ensemble,
            set_ens['tis'],
            cycle)
        last_paths.append(initial_path)
        yield accept, initial_path, status, path_ensemble


def initiate_path_ensemble_kick(ensemble, tis_settings, cycle):
    """Run the "kick" initiate for a given ensemble.

    Parameters
    ----------
    ensemble: dict
        It contains:

        * `path_ensemble` : object like :py:class:`.PathEnsemble`
          The path ensemble to create an initial path for.
        * `system` : object like :py:class:`.System`
          The system is used here since we need access to the temperature
          and to the particle list.
        * `order_function` : object like :py:class:`.OrderParameter`
          The class used for obtaining the order parameter(s).
        * `engine` : object like :py:class:`.EngineBase`
          The engine to use for propagating a path.
        * `rgen` : object like :py:class:`.RandomGenerator`
          This is the random generator that will be used.

    tis_settings : dict
        This dictionary contains the TIS settings.
    cycle : integer, optional
        The cycle number we are initiating at, typically this will be 0
        which is the default value.

    Returns
    -------
    out[0] : boolean
        True if the initial path was accepted.
    out[1] : object like py:class:`.PathBase`
        The initial path.
    out[2] : string
        The status of the path.

    """
    initial_path = None
    logger.info('Will generate initial path by kicking')
    ensemble['engine'].exe_dir = \
        ensemble['path_ensemble'].directory['generate']

    for attempt in generate_initial_path_kick(ensemble, tis_settings):
        if attempt[0] is True:
            initial_path = attempt[-1]
            break
    accept = True
    status = 'ACC'
    ensemble['path_ensemble'].add_path_data(initial_path, status, cycle)
    # Ask the engine to do clean up after the initialisation:
    ensemble['engine'].clean_up()

    return accept, initial_path, status


def generate_initial_path_kick(ensemble, tis_settings):
    """Generate an initial path with the kicking method.

    This function will generate an initial path by repeatedly kicking a
    phase-space point until the middle interface is crossed.
    The point before and after kicking are stored, so when the
    middle interface is crossed we have two points we can integrate
    forward and backward in time. This function is intended for use with
    TIS. For use with RETIS one should set the appropriate
    `tis_settings` so that the starting conditions are fine (i.e. for
    the [0^-] ensemble it might be different for the other ensembles).

    Parameters
    ----------
    ensemble: dict.
        It contains:

        * `path_ensemble` : object like :py:class:`.PathEnsemble`
          The path ensemble to create an initial path for.
        * `system` : object like :py:class:`.System`
          System is used here since we need access to the temperature
          and to the particle list.
        * `order_function` : object like :py:class:`.OrderParameter`
          The class used for obtaining the order parameter(s).
        * `engine` : object like :py:class:`.EngineBase`
          The engine to use for propagating a path.
        * `rgen` : object like :py:class:`.RandomGenerator`
          This is the random generator that will be used.

    tis_settings : dict
        This dictionary contains settings for TIS. Explicitly used here:

        * `start_cond`: string, starting condition, 'L'eft or 'R'ight.
        * `maxlength`: integer, maximum allowed length of paths.

    Yields
    ------
    out[0] : boolean
        True if we are finished.
    out[1] : string
        The current status of this method.
    out[2] : object like :py:class:`.PathBase` or None
        If we are finished (and successful) this is the generated
        initial path. Otherwise, this is equal to None.

    """
    engine = ensemble['engine']
    interfaces = ensemble['path_ensemble'].interfaces
    rgen = ensemble['path_ensemble'].rgen
    initial_system = ensemble['system'].copy()
    middle = interfaces[1]
    kick_required = True
    if set(ensemble['path_ensemble'].start_condition) == set(['R', 'L']):
        initial_op = _order_value(initial_system.order)
        logger.info("Initial OP: %s", initial_op)

        # If we are allowed to end anywhere and the initial system is between
        # the two interfaces, don't kick.
        if interfaces[0] < initial_op < interfaces[-1]:
            kick_required = False
        # Else kick to the closest interfaces
        elif abs(interfaces[0]-initial_op) < abs(interfaces[-1]-initial_op):
            middle = interfaces[0]
        else:
            middle = interfaces[-1]
    maxlen = tis_settings['maxlength']

    while True:
        # Start from the initial system:
        ensemble['system'] = initial_system.copy()
        ensemble['interfaces'] = interfaces
        if kick_required:
            logger.info('Searching crossing with middle interface')
            left_point, right_point = engine.kick_across_middle(ensemble,
                                                                middle,
                                                                tis_settings)
        else:
            logger.info('Not searching a crossing with the middle interface')
            right_point = initial_system.copy()
            # Set starting point for both propagators to the current system
            left_point = initial_system.copy()

        logger.info('Propagating from crossing points')
        # kick_across_middle will return two points, one immediately
        # left of the interface and one immediately right of the
        # interface. So we have two points (`leftpoint` and the
        # current `system.particles`). We then propagate the current
        # phase point forward:
        logger.info('Propagating forward from crossing points '
                    'for initial path')
        path_forw = Path(rgen=rgen, maxlen=maxlen)
        ensemble['system'] = right_point
        success, msg = engine.propagate(path_forw, ensemble,
                                        reverse=False)

        if 'exp' in tis_settings.get('shooting_move', 'sh'):
            path_forw.status = 'ACC'
            path_forw.generated = ('ki', 0, 0, 0)
            yield True, 'Initial path generated.', path_forw

        if not success:
            logger.warning('Forward path not successful: %s', msg)
            yield False, f'Forward path failed: {msg}', None
            continue
        # And we propagate the `leftpoint` backward:
        path_back = Path(rgen=rgen, maxlen=maxlen)
        ensemble['system'] = left_point
        success, msg = engine.propagate(path_back, ensemble,
                                        reverse=True)
        if not success:
            logger.warning('Backward path not successful: %s', msg)
            yield False, f'Backward path failed: {msg}', None
            continue
        # Merge backward and forward, here we do not set maxlen since
        # both backward and forward may have this length:
        initial_path = paste_paths(path_back, path_forw, overlap=False)
        if initial_path.length >= maxlen:
            logger.warning('Initial path too long (exceeded "MAXLEN")')
            yield False, 'Initial path was too long.', None
            continue
        start, end, _, _ = initial_path.check_interfaces(interfaces)
        # OK, now it is time to check the path:
        # 0) We can start at the starting condition, pass the middle
        #    and continue all the way to the end - perfect!
        # 1) We can start at the starting condition, pass the middle
        #    and return to starting condition - this is perfectly fine.
        # 2) We can start at the wrong interface, pass the middle and
        #    end at the same (wrong) interface - we fix this by
        #    additional TIS moves.
        # 3) We can start at wrong interface and end at the starting
        #    condition. To fix this, we reverse the path.
        if start in set(ensemble['path_ensemble'].start_condition):
            # Case 0) & 1):
            break
        # Now we do the other cases:
        if end in set(ensemble['path_ensemble'].start_condition):
            # Case 3) (and start != start_cond):
            logger.info('Initial path is in the wrong direction')
            initial_path = initial_path.reverse(ensemble['order_function'])
            logger.info('Initial path has been reversed!')
            break
        if end == start:
            # Case 2):
            logger.info('Initial path starts and ends at wrong interface')
            logger.info('Will perform TIS moves to fix it!')
            yield False, 'Trying to fix path by TIS moves', None

            initial_path = fix_path_by_tis(initial_path, ensemble,
                                           tis_settings)
            break
        # This is reached if the start/end conditions
        # has been modified or are inconsistent.
        logger.warning('Could not generate initial path, will retry!')
        yield (False, 'Could not generate initial path will retry!',
               None)
        continue
    initial_path.status = 'ACC'
    initial_path.generated = ('ki', 0, 0, 0)
    if tis_settings.get('high_accept', False):
        move = tis_settings.get('shooting_move', 'sh')
        if move == 'wf':
            int_cap = tis_settings.get('interface_cap', interfaces[2])
            interfaces = [interfaces[0], interfaces[1], int_cap]
        if move in ('ss', 'wf'):
            initial_path.weight = compute_weight(initial_path,
                                                 interfaces,
                                                 move)
    yield True, 'Initial path generated.', initial_path


def _get_help(start_cond, interfaces):
    """Define some helper methods for :py:func:`.fix_path_by_tis`.

    This method returns two methods that :py:func:`.fix_path_by_tis`
    can use to determine if a new path is an improvement compared to
    the current path and if the "fixing" is done.

    Parameters
    ----------
    start_cond : string
        The starting condition (from the TIS settings). Left or Right.
    interfaces : list of floats
        The interfaces, on form ``[left, middle, right]``.

    Returns
    -------
    out[0] : method
        The method which determines if a new path represents an
        improvement over the current path.
    out[1] : method
        The method which determines if we are done, that is if we
        can accept the current path.

    """
    left, middle, right = interfaces
    improved, done = None, None
    if start_cond == 'R' or middle == right:
        def improved_r(newp, current):
            """Return True if new path is an improvement."""
            return (newp.ordermax[0] > current.ordermax[0] and
                    newp.ordermin[0] < middle)

        def done_r(path):
            """Return True if the path can be accepted."""
            return path.ordermax[0] > right

        improved = improved_r
        done = done_r
    elif start_cond == 'L' or middle == left:
        def improved_l(newp, current):
            """Return True if new path is an improvement."""
            return (newp.ordermin[0] < current.ordermin[0] and
                    newp.ordermax[0] > middle)

        def done_l(path):
            """Return True if the path can be accepted."""
            return path.ordermin[0] < left

        improved = improved_l
        done = done_l
    else:
        logger.error('Unknown start condition (should be "R" or "L")')
        raise ValueError('Unknown start condition (should be "R" or "L")')
    return improved, done


def fix_path_by_tis(initial_path, ensemble, tis_settings):
    """Fix a path that starts and ends at the wrong interfaces.

    The given path is amended by performing TIS moves (shooting and
    time reversal). Note that this function is intended to be used
    in connection with the initialisation.

    Parameters
    ----------
    initial_path : object like :py:class:`.PathBase`
        This is the initial path to fix. It starts & ends at the
        wrong interface.
    ensemble: dict.
        It contains:

        * `path_ensemble` : object like :py:class:`.PathEnsemble`
          The path ensemble to create an initial path for.
        * `system` : object like :py:class:`.System`
          System is used here since we need access to the temperature
          and to the particle list.
        * `order_function` : object like :py:class:`.OrderParameter`
          The class used for obtaining the order parameter(s).
        * `engine` : object like :py:class:`.EngineBase`
          The engine to use for propagating a path.
        * `rgen` : object like :py:class:`.RandomGenerator`
          This is the random generator that will be used.

    tis_settings : dict
        Settings for TIS method, here we explicitly use:

        * `start_cond`: string which defines the start condition.
        * `maxlength`: integer which gives the maximum allowed path
          length.

        Note that we here explicitly set some local TIS settings for
        use in the `make_tis_step` function.

    Returns
    -------
    out : object like :py:class:`.PathBase`
        The amended path.

    """
    logger.debug('Attempting to fix path by running TIS moves.')

    path_ensemble = ensemble['path_ensemble']
    local_tis_settings = copy.deepcopy(tis_settings)
    local_tis_settings['allowmaxlength'] = True
    local_tis_settings['aimless'] = True
    improved, check_ok = _get_help(path_ensemble.start_condition,
                                   path_ensemble.interfaces)

    backup_path = True
    path_ok = False

    path_ensemble.last_path = initial_path

    while not path_ok:
        logger.debug('Performing a TIS move to improve the initial path')
        if backup_path:  # Move the initial path to safe place:
            logger.debug('Moving initial_path')
            path_ensemble.move_path_to_generate(initial_path, '_')
            backup_path = False

        accept, trial, _ = time_reversal(initial_path,
                                         ensemble['order_function'],
                                         path_ensemble.interfaces,
                                         path_ensemble.start_condition)

        if not accept:
            target = path_ensemble.interfaces[1]
            best_frame = ensemble['system']

            for phasepoint in initial_path.phasepoints:
                if abs(_order_value(best_frame.order) - target) > \
                        abs(_order_value(phasepoint.order) - target):
                    best_frame = phasepoint.copy()

            accept, trial, _ = shoot(ensemble,
                                     local_tis_settings,
                                     path_ensemble.start_condition,
                                     best_frame)

        if accept:
            if improved(trial, initial_path):
                logger.debug('TIS move improved path.')
                initial_path = trial
                path_ensemble.last_path = initial_path
                backup_path = True
            else:
                logger.debug('TIS move did not improve path')
            path_ok = check_ok(initial_path)

    initial_path.generated = ('ki', 0, 0, 0)
    initial_path.status = 'ACC'
    return initial_path
