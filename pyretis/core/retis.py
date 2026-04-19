# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This module contains functions for RETIS.

This module defines functions that are needed to perform Replica
Exchange Transition Interface Sampling (RETIS). The RETIS algorithm
was first described by van Erp [RETIS]_.


Methods defined here
~~~~~~~~~~~~~~~~~~~~

make_retis_step (:py:func:`.make_retis_step`)
    Function to select and execute the RETIS move.

retis_moves (:py:func:`.retis_moves`)
    Function to perform RETIS swapping moves - it selects what scheme
    to use, i.e. ``[0^-] <-> [0^+], [1^+] <-> [2^+], ...`` or
    ``[0^+] <-> [1^+], [2^+] <-> [3^+], ...``.

retis_swap_wrapper (:py:func:`.retis_swap_wrapper`)
    The function that actually swaps two path ensembles.
    This function decides which swap will be done:
    it can be a typical RETIS swap (retis_swap), or a swap between two PPTIS
    ensembles with limited memory (prretis_swap).

retis_swap (:py:func:`.retis_swap`)
    The function that actually swaps two path ensembles.
    A swap between two typical RETIS ensembles, not PPTIS ensembles.

repptis_swap (:py:func:`.repptis_swap`)
    The function that actually swaps two path ensembles.
    A swap between two PPTIS ensembles with truncated paths,
    so this is not the typical RETIS swap.

retis_swap_zero (:py:func:`.retis_swap_zero`)
    The function that performs the swapping for the
    ``[0^-] <-> [0^+]`` swap.

high_acc_swap (:py:func:`.high_acc_wap`)
    The function coputes if a path generated via SS can be accepted
    for swapping in accordance to super detail balance.

References
~~~~~~~~~~
.. [RETIS] T. S. van Erp,
   Phys. Rev. Lett. 98, 26830 (2007),
   http://dx.doi.org/10.1103/PhysRevLett.98.268301

"""
import copy
import logging
from pyretis.core.common import (null_move,
                                 compute_weight,
                                 priority_checker)
from pyretis.core.tis import make_tis, paste_paths
from pyretis.inout.restart import write_ensemble_restart

logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())

# pylint: disable=W0106

__all__ = ['high_acc_swap',
           'make_retis_step',
           'repptis_swap',
           'retis_moves',
           'retis_swap_wrapper',
           'retis_swap',
           'retis_swap_zero']


def make_retis_step(ensembles, rgen, settings, cycle):
    """Determine and execute the appropriate RETIS move.

    Here we will determine what kind of RETIS moves we should do.
    We have two options:

    1) Do the RETIS swapping moves. This is done by calling
       :py:func:`.retis_moves`.
    2) Do TIS moves, either for all ensembles or for just one, based on
       values of relative shoot frequencies. This is done by calling
       :py:func:`.retis_tis_moves`.

    This function will just determine and execute the appropriate move
    (1 or 2) based on the given swapping frequencies in the `settings`
    and drawing a random number from the random number generator `rgen`.

    Parameters
    ----------
    ensembles : list of dicts of objects
        This is a list of the ensembles we are using in the RETIS method.
        It contains:

        * `path_ensemble`: object like :py:class:`.PathEnsemble`
          This is used for storing results for the simulation. It
          is also used for defining the interfaces for this simulation.
        * `system`: object like :py:class:`.System`
          System is used here since we need access to the temperature
          and to the particle list.
        * `order_function`: object like :py:class:`.OrderParameter`
          The class used for calculating the order parameter(s).
        * `engine`: object like :py:class:`.EngineBase`
          The engine to use for propagating a path.

    rgen : object like :py:class:`.RandomGenerator`
        This is a random generator. Here we assume that we can call
        `rgen.rand()` to draw random uniform numbers.
    settings : dict
        Contains the settings for the simulation.
    cycle : integer
        The current cycle number.

    Returns
    -------
    out : generator
        This generator yields the results after performing the
        RETIS moves.

    """
    prio_skip = priority_checker(ensembles, settings)
    swap_freq = settings['retis']['swapfreq']
    swap = False if True in prio_skip else rgen.rand() < swap_freq

    if swap:
        # Do RETIS moves
        logger.info('Performing RETIS swapping move(s).')
        results = retis_moves(ensembles, rgen, settings, cycle)
    else:
        logger.info('Performing RETIS TIS move(s)')
        results = make_tis(ensembles, rgen, settings, cycle)
    return results


def retis_moves(ensembles, rgen, settings, cycle):
    """Perform RETIS moves on the given ensembles.

    This function will perform RETIS moves on the given ensembles.
    First we have two strategies based on
    `settings['retis']['swapsimul']`:

    1) If `settings['retis']['swapsimul']` is True we will perform
       several swaps, either ``[0^-] <-> [0^+], [1^+] <-> [2^+], ...``
       or ``[0^+] <-> [1^+], [2^+] <-> [3^+], ...``. Which one of these
       two swap options we use is determined randomly and they have
       equal probability.

    2) If `settings['retis']['swapsimul']` is False we will just
       perform one swap for randomly chosen ensembles, i.e. we pick a
       random ensemble and try to swap with the ensemble to the right.
       Here we may also perform null moves if the
       `settings['retis']['nullmove']` specifies so.

    Parameters
    ----------
    ensembles : list of dicts of objects
        Each contains:

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
        This is a random generator. Here we assume that we can call
        `rgen.rand()` to draw random uniform numbers.
    settings : dict
        This dict contains the settings for the RETIS method.
    cycle : integer
        The current cycle number.

    Yields
    ------
    out : dict
        This dictionary contains the result of the RETIS moves.

    """
    if settings['retis']['swapsimul']:
        # Here we have two schemes:
        # 1) scheme == 0: [0^-] <-> [0^+], [1^+] <-> [2^+], ...
        # 2) scheme == 1: [0^+] <-> [1^+], [2^+] <-> [3^+], ...
        if len(ensembles) < 3:
            # Small number of ensembles, can only do the [0^-] <-> [0^+] swap:
            scheme = 0
        else:
            scheme = 0 if rgen.rand() < 0.5 else 1
        for i in range(scheme, len(ensembles) - 1, 2):
            idx = ensembles[i]['path_ensemble'].ensemble_number
            accept, trial, status = retis_swap_wrapper(
                ensembles, idx, settings, cycle
            )
            result = {
                'ensemble_number': idx,
                'mc-move': 'swap',
                'status': status,
                'trial': trial[0],
                'accept': accept,
                'swap-with': idx + 1,
            }
            yield result
            result = {
                'ensemble_number': idx + 1,
                'mc-move': 'swap',
                'status': status,
                'trial': trial[1],
                'accept': accept,
                'swap-with': idx,
            }
            yield result
        # We might have missed some ensembles in the two schemes.
        # Here, we do null moves in these, if requested:
        if settings['retis']['nullmoves']:
            if len(ensembles) % 2 != scheme:  # Missed last ensemble:
                # This is perhaps strange but it is equivalent to:
                # (scheme == 0 and len(ensembles) % 2 != 0) or
                # (scheme == 1 and len(ensembles) % 2 == 0)
                accept, trial, status = null_move(ensembles[-1], cycle)
                result = {
                    'ensemble_number':
                        ensembles[-1]['path_ensemble'].ensemble_number,
                    'mc-move': 'nullmove',
                    'status': status,
                    'trial': trial,
                    'accept': accept,
                }
                yield result
            # We always miss the first ensemble in scheme 1:
            if scheme == 1:
                accept, trial, status = null_move(ensembles[0], cycle)
                result = {
                    'ensemble_number':
                        ensembles[0]['path_ensemble'].ensemble_number,
                    'mc-move': 'nullmove',
                    'status': status,
                    'trial': trial,
                    'accept': accept,
                }
                yield result
    else:  # Just swap two ensembles:
        idx = rgen.random_integers(0, len(ensembles) - 2)
        accept, trial, status = retis_swap_wrapper(
            ensembles, idx, settings, cycle
        )
        result = {
            'ensemble_number': idx,
            'mc-move': 'swap',
            'status': status,
            'trial': trial[0],
            'accept': accept,
            'swap-with': idx + 1,
        }
        yield result
        result = {
            'ensemble_number': idx + 1,
            'mc-move': 'swap',
            'status': status,
            'trial': trial[1],
            'accept': accept,
            'swap-with': idx,
        }
        yield result
        # Do null moves in the other ensembles, if requested:
        if settings['retis']['nullmoves']:
            for ensemble in ensembles:
                idx2 = ensemble['path_ensemble'].ensemble_number
                if idx2 not in (idx, idx + 1):
                    accept, trial, status = null_move(ensemble, cycle)
                    result = {
                        'ensemble_number': idx2,
                        'mc-move': 'nullmove',
                        'status': status,
                        'trial': trial,
                        'accept': accept,
                    }
                    yield result


def retis_swap_wrapper(ensembles, idx, settings, cycle):
    """Swap the last accepted paths of two TIS or two PPTIS ensembles.

    This function selects the correct swap function to use, based on the
    simulation task (retis or repptis).

    Parameters
    ----------
    ensembles : list of dictionaries of objects
        This is a list of the ensembles we are using in the RE(PP)TIS method.
        Each one contains:

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

    idx : integer
        Definition of what path ensembles to swap. We will swap
        `ensembles[idx]` with `ensembles[idx+1]`.
    settings : dict
        This dict contains the settings for the RE(PP)TIS method.
    cycle : integer
        Current cycle number.

    Returns
    -------
    out[0] : boolean
        Should the path be accepted or not?
    out[1] : list of object like :py:class:`.PathBase`
        The trial paths.
    out[2] : string
        The status for the trial paths.

    """
    if settings['simulation']['task'] == "repptis":
        return repptis_swap(ensembles, idx, settings, cycle)

    return retis_swap(ensembles, idx, settings, cycle)


def retis_swap(ensembles, idx, settings, cycle):
    """Perform a RETIS swapping move for two ensembles.

    These ensembles are not PPTIS ensembles.

    The RETIS swapping move will attempt to swap accepted paths between
    two ensembles in the hope that path from [i^+] is an acceptable path
    for [(i+1)^+] as well. We have two cases:

    1) If we try to swap between [0^-] and [0^+] we need to integrate
       the equations of motion.
    2) Otherwise, we can just swap and accept if the path from [i^+] is
       an acceptable path for [(i+1)^+]. The path from [(i+1)^+] is
       always acceptable for [i^+] (by construction).

    Parameters
    ----------
    ensembles : list of dictionaries of objects
        This is a list of the ensembles we are using in the RETIS method.
        Each one contains:

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

    idx : integer
        Definition of what path ensembles to swap. We will swap
        `ensembles[idx]` with `ensembles[idx+1]`. If `idx == 0` we have
        case 1) defined above.
    settings : dict
        This dict contains the settings for the RETIS method.
    cycle : integer
        Current cycle number.

    Returns
    -------
    out[0] : boolean
        Should the path be accepted or not?
    out[1] : list of object like :py:class:`.PathBase`
        The trial paths.
    out[2] : string
        The status for the trial paths.

    Note
    ----
    Note that path.generated is **NOT** updated here. This is because
    we are just swapping references and not the paths. In case the
    swap is rejected updating this would invalidate the last accepted
    path.

    """
    logger.info("Swapping: %s <-> %s",
                ensembles[idx]['path_ensemble'].ensemble_name,
                ensembles[idx+1]['path_ensemble'].ensemble_name)

    if idx == 0:
        return retis_swap_zero(ensembles, settings, cycle)
    path_ensemble1 = ensembles[idx]['path_ensemble']
    path_ensemble2 = ensembles[idx+1]['path_ensemble']

    path1 = path_ensemble1.last_path
    path2 = path_ensemble2.last_path
    ens_moves = [settings['ensemble'][i]['tis'].get('shooting_move', 'sh')
                 for i in [idx, idx+1]]
    intf_w = [list(i) for i in (path_ensemble1.interfaces,
                                path_ensemble2.interfaces)]
    for i, j in enumerate([settings['ensemble'][k] for k in (idx, idx+1)]):
        if ens_moves[i] == 'wf':
            intf_w[i][2] = j['tis'].get('interface_cap', intf_w[i][2])

    # Check if path1 can be accepted in ensemble 2:
    cross = path1.check_interfaces(path_ensemble2.interfaces)[-1]
    accept = cross[1]
    status = 'ACC' if accept else 'NCR'

    if accept and settings['tis'].get('high_accept', False):
        if 'ss' in ens_moves or 'wf' in ens_moves:
            accept, status = high_acc_swap([path1, path2],
                                           ensembles[idx]['rgen'],
                                           intf_w[0], intf_w[1],
                                           ens_moves)
    if accept:  # Accept the swap:
        logger.info('Swap was accepted.')
        # To avoid overwriting files, we move the paths to the
        # generate directory here.
        path_ensemble2.move_path_to_generate(path1)
        path_ensemble1.move_path_to_generate(path2)
        for i, path, path_ensemble, flag in ((0, path2, path_ensemble1, 's+'),
                                             (1, path1, path_ensemble2, 's-')):
            path.status = status
            ens_set = settings['ensemble'][idx + i]
            move = ens_moves[i]
            path.weight = compute_weight(path, intf_w[i], move)\
                if (ens_set['tis'].get('high_accept', False) and
                    move in ('wf', 'ss')) else 1

            path.set_move(flag) if path.get_move() != 'ld' else \
                path.set_move('ld')

            # Then moved into accepted directory by the `add_path_data` below.
            path_ensemble.add_path_data(path, status, cycle=cycle)
            if cycle % ens_set.get('output', {}).get('restart-file', 1) == 0:
                write_ensemble_restart(ensembles[idx+i], ens_set)

        return accept, (path2, path1), status

    logger.info('Swap was rejected. (%s)', status)

    # Make shallow copies:
    trial1 = copy.copy(path2)
    trial2 = copy.copy(path1)
    trial1.set_move('s+')  # Came from right.
    trial2.set_move('s-')  # Came from left.
    # Calculate weights:
    for i, trial in ((0, trial1), (1, trial2)):
        trial.weight = compute_weight(trial, intf_w[i], ens_moves[i])\
            if (settings['ensemble'][idx+i]['tis'].get('high_accept', False)
                and ens_moves[i] in ('wf', 'ss')) else 1.
    path_ensemble1.add_path_data(trial1, status, cycle=cycle)
    path_ensemble2.add_path_data(trial2, status, cycle=cycle)

    return accept, (trial1, trial2), status


def retis_swap_zero(ensembles, settings, cycle):
    """Perform the RETIS swapping for ``[0^-] <-> [0^+]`` swaps.

    The RETIS swapping move for ensembles [0^-] and [0^+] requires some
    extra integration. Here we are generating new paths for [0^-] and
    [0^+] in the following way:

    1) For [0^-] we take the initial point in [0^+] and integrate
       backward in time. This is merged with the second point in [0^+]
       to give the final path. The initial point in [0^+] starts to the
       left of the interface and the second point is on the right
       side - i.e. the path will cross the interface at the end points.
       If we let the last point in [0^+] be called ``A_0`` and the
       second last point ``B``, and we let ``A_1, A_2, ...`` be the
       points on the backward trajectory generated from ``A_0`` then
       the final path will be made up of the points
       ``[..., A_2, A_1, A_0, B]``. Here, ``B`` will be on the right
       side of the interface and the first point of the path will also
       be on the right side.

    2) For [0^+] we take the last point of [0^-] and use that as an
       initial point to generate a new trajectory for [0^+] by
       integration forward in time. We also include the second last
       point of the [0^-] trajectory which is on the left side of the
       interface. We let the second last point be ``B`` (this is on the
       left side of the interface), the last point ``A_0`` and the
       points generated from ``A_0`` we denote by ``A_1, A_2, ...``.
       Then the resulting path will be ``[B, A_0, A_1, A_2, ...]``.
       Here, ``B`` will be on the left side of the interface and the
       last point of the path will also be on the left side of the
       interface.

    Parameters
    ----------
    ensembles : list of dictionaries of objects
        This is a list of the ensembles we are using in the RETIS method.
        It contains:

        * `path_ensemble`: object like :py:class:`.PathEnsemble`
          This is used for storing results for the simulation. It
          is also used for defining the interfaces for this simulation.
        * `system`: object like :py:class:`.System`
          System is used here since we need access to the temperature
          and to the particle list.
        * `order_function`: object like :py:class:`.OrderParameter`
          The class used for calculating the order parameter(s).
        * `engine`: object like :py:class:`.EngineBase`
          The engine to use for propagating a path.

    settings : dict
        This dict contains the settings for the RETIS method.
    cycle : integer
        The current cycle number.

    Returns
    -------
    out : string
        The result of the swapping move.

    """
    path_ensemble0 = ensembles[0]['path_ensemble']
    path_ensemble1 = ensembles[1]['path_ensemble']
    engine0, engine1 = ensembles[0]['engine'], ensembles[1]['engine']
    maxlen0 = settings['ensemble'][0]['tis']['maxlength']
    maxlen1 = settings['ensemble'][1]['tis']['maxlength']

    ens_moves = [settings['ensemble'][i]['tis'].get('shooting_move', 'sh')
                 for i in [0, 1]]
    intf_w = [list(i) for i in (path_ensemble0.interfaces,
                                path_ensemble1.interfaces)]
    for i, j in enumerate([settings['ensemble'][k] for k in (0, 1)]):
        if ens_moves[i] == 'wf':
            intf_w[i][2] = j['tis'].get('interface_cap', intf_w[i][2])

    # 0. check if MD is allowed
    allowed = (path_ensemble0.last_path.get_end_point(
                path_ensemble0.interfaces[0],
                path_ensemble0.interfaces[-1]) == 'R')
    # 1. Generate path for [0^-] from [0^+]:
    # We generate from the first point of the path in [0^+]:
    logger.debug('Swapping [0^-] <-> [0^+]')
    logger.debug('Creating path for [0^-]')
    system = path_ensemble1.last_path.phasepoints[0].copy()
    logger.debug('Initial point is: %s', system)
    ensembles[0]['system'] = system
    # Propagate it backward in time:
    path_tmp = path_ensemble1.last_path.empty_path(maxlen=maxlen1-1)
    if allowed:
        logger.debug('Propagating for [0^-]')
        engine0.propagate(path_tmp, ensembles[0], reverse=True)
    else:
        logger.debug('Not propagating for [0^-]')
        path_tmp.append(system)
    path0 = path_tmp.empty_path(maxlen=maxlen0)
    for phasepoint in reversed(path_tmp.phasepoints):
        path0.append(phasepoint)
    # Add second point from [0^+] at the end:
    logger.debug('Adding second point from [0^+]:')
    # Here we make a copy of the phase point, as we will update
    # the configuration and append it to the new path:
    phase_point = path_ensemble1.last_path.phasepoints[1].copy()
    logger.debug('Point is %s', phase_point)
    engine0.dump_phasepoint(phase_point, 'second')
    path0.append(phase_point)
    if path0.length == maxlen0:
        path0.status = 'BTX'
    elif path0.length < 3:
        path0.status = 'BTS'
    elif ('L' not in set(path_ensemble0.start_condition) and
          'L' in path0.check_interfaces(path_ensemble0.interfaces)[:2]):
        path0.status = '0-L'
    else:
        path0.status = 'ACC'

    # 2. Generate path for [0^+] from [0^-]:
    logger.debug('Creating path for [0^+] from [0^-]')
    # This path will be generated starting from the LAST point of [0^-] which
    # should be on the right side of the interface. We will also add the
    # SECOND LAST point from [0^-] which should be on the left side of the
    # interface, this is added after we have generated the path and we
    # save space for this point by letting maxlen = maxlen1-1 here:
    path_tmp = path0.empty_path(maxlen=maxlen1-1)
    # We start the generation from the LAST point:
    # Again, the copy below is not needed as the propagate
    # method will not alter the initial state.
    system = path_ensemble0.last_path.phasepoints[-1].copy()
    if allowed:
        logger.debug('Initial point is %s', system)
        ensembles[1]['system'] = system
        logger.debug('Propagating for [0^+]')
        engine1.propagate(path_tmp, ensembles[1], reverse=False)
        # Ok, now we need to just add the SECOND LAST point from [0^-] as
        # the first point for the path:
        path1 = path_tmp.empty_path(maxlen=maxlen1)
        phase_point = path_ensemble0.last_path.phasepoints[-2].copy()
        logger.debug('Add second last point: %s', phase_point)
        engine1.dump_phasepoint(phase_point, 'second_last')
        path1.append(phase_point)
        path1 += path_tmp  # Add rest of the path.
    else:
        path1 = path_tmp
        path1.append(system)
        logger.debug('Skipping propagating for [0^+] from L')

    if path_ensemble1.last_path.get_move() != 'ld':
        path0.set_move('s+')
    else:
        path0.set_move('ld')

    if path_ensemble0.last_path.get_move() != 'ld':
        path1.set_move('s-')
    else:
        path1.set_move('ld')
    if path1.length >= maxlen1:
        path1.status = 'FTX'
    elif path1.length < 3:
        path1.status = 'FTS'
    else:
        path1.status = 'ACC'
    logger.debug('Done with swap zero!')

    # Final checks:
    accept = path0.status == 'ACC' and path1.status == 'ACC'
    status = 'ACC' if accept else (path0.status if path0.status != 'ACC' else
                                   path1.status)
    # High Acceptance swap is required when Wire Fencing are used
    if accept and settings['tis'].get('high_accept', False):
        if 'wf' in ens_moves:
            accept, status = high_acc_swap([path1, path_ensemble1.last_path],
                                           ensembles[0]['rgen'],
                                           intf_w[0],
                                           intf_w[1],
                                           ens_moves)

    for i, path, path_ensemble, flag in ((0, path0, path_ensemble0, 's+'),
                                         (1, path1, path_ensemble1, 's-')):
        if not accept and path.status == 'ACC':
            path.status = status

        # These should be 1 unless length of paths equals 3.
        # This technicality is not yet fixed. (An issue is open as a reminder)

        ens_set = settings['ensemble'][i]
        move = ens_moves[i]
        path.weight = compute_weight(path, intf_w[i], move)\
            if (ens_set['tis'].get('high_accept', False) and
                move in ('wf', 'ss')) else 1

        if accept:
            path_ensemble.move_path_to_generate(path)
        else:
            logger.debug("Rejected swap path in [0^%s], %s", flag[:-1], status)
        path_ensemble.add_path_data(path, path.status, cycle=cycle)
        if cycle % ens_set.get('output', {}).get('restart-file', 1) == 0:
            write_ensemble_restart(ensembles[i], settings['ensemble'][i])

    return accept, (path0, path1), status


def repptis_swap(ensembles, idx, settings, cycle, inc_sh=True):
    """Perform a replica exchange move between adjacent PPTIS ensembles.

    First, propagation directions are randomly chosen for both paths, and it is
    checked whether the pathtypes (i.e. LMR, RMR, ...) are compatible. If not,
    the move is rejected. If the proposed directions are compatible with the
    pathtypes, the path segments in the overlapping region are cut, and they
    are propagated in the proposed directions.

    When compatible, the paths are created as follows:
    1) new [i^+-] path from old [(i+1)^+-] path. The path segment in the 'LM'
    region is cut out of the old [(i+1)^+-] path (with one point left of the
    L interface, and one point right of the M interface). This path segment is
    then extended by propagating the phasepoint left of the L interface in the
    proposed direction till a crossing condition of the [i^+-] interface occurs

    2) new [(i+1)^+-] path from old [i^+-] path. The path segment in the 'MR'
    region is cut out of the old [i^+-] path (with one point left of the M
    interface, and one point right of the R interface). This path segment is
    then extended by propagating the phasepoint right of the R interface in the
    proposed direction until a crossing condition of the [(i+1)^+-] interface
    occurs.

    Floating point precision problem (for GRO):
    The phase point to be shot from, is first dumped to a .gro file, where
    precision goes from ~9 decimals to 3 decimals. This impacts the order
    parameter of the dumped phasepoint, and if it was very closely located to
    an interface, it might even switch sides w.r.t. that interface. The problem
    has an easy fix. We include the shooting point in the overlap_path, such
    that the order parameter does not change. When we shoot from this point,
    the order parameter will still change due to precision decrease, and the
    first point of the propagated_path will be slightly different. We discard
    this point in the propagated_path, and keep the one from the overlap_path.
    Whichever shootpoint you keep, there will be a mismatch. The choice we make
    does not cause trouble with the swapping move. The above described behavior
    (i.e. the fix) is enabled with the paramter inc_sh set to True, being
    the standard behavior.

    Parameters
    ----------
    ensembles : list of dictionaries of objects
        This is a list of the ensembles we are using in the RETIS method.
        It contains:

        * `path_ensemble`: object like :py:class:`.PathEnsemble`
        This is used for storing results for the simulation. It
        is also used for defining the interfaces for this simulation.
        * `system`: object like :py:class:`.System`
        System is used here since we need access to the temperature
        and to the particle list.
        * `order_function`: object like :py:class:`.OrderParameter`
        The class used for calculating the order parameter(s).
        * `engine`: object like :py:class:`.EngineBase`
        The engine to use for propagating a path.

    idx : int
        The index of the ensemble that will be swapped.
    settings : dict
        This dict contains the settings for the RETIS method.
    cycle : integer
        The current cycle number.
    inc_sh : bool
        If True, the phasepoint that is used for shooting is included in the
        overlap path. If False, it will be included in the newly propagated
        path (with lower .gro precision).

    Returns
    -------
    out[0] : boolean
        Should the paths be accepted or not?
    out[1] : list of object like :py:class:`.PathBase`
        The trial paths.
    out[2] : string
        The status for the trial paths.

    """
    msg = f"Swapping: {ensembles[idx]['path_ensemble'].ensemble_name} <-> "
    msg += f"{ensembles[idx+1]['path_ensemble'].ensemble_name}"
    logger.info(msg)

    # if idx == 0, we perform a swap_zero move
    if idx == 0:
        # Since we allow RML paths in [0+-'], we need to do an extra check,
        # as this pathtype will result in a bogus swap...
        # NB: if we don't check for this, the move will be rejected anyways,
        # as the first propagated ph will be OOB for [0-], which results in
        # an "FTS" or "BTS" status. However, it's cleaner to reject the move
        # in advance...
        start, end, _, _ =\
            ensembles[1]['path_ensemble'].last_path.check_interfaces(
                ensembles[1]['path_ensemble'].interfaces)
        ptype = ""+start+"M"+end
        if ptype == "RML":
            msg = "Rejecting swap: RML in [0+-'] cannot be swapped with [0-]!"
            logger.info(msg)
            return reject_repptis_swap_wrong_prop(
                ensembles[0]['path_ensemble'], ensembles[1]['path_ensemble'],
                cycle, settings, idx, ensembles)
        return retis_swap_zero(ensembles, settings, cycle)

    # else we perform a repptis swap move
    path_ensemble0 = ensembles[idx]['path_ensemble']
    path_ensemble1 = ensembles[idx+1]['path_ensemble']

    # Choose propagation directions for the two paths, and check if the
    # proposed move is allowed or not
    allowed, propdir0, propdir1, path0_type, path1_type = \
        re_propagation_directions(ensembles, idx)

    # If not allowed, we immediately reject the swap
    if not allowed:
        return reject_repptis_swap_wrong_prop(path_ensemble0, path_ensemble1,
                                              cycle, settings, idx, ensembles)

    # If allowed, the swap move is attempted
    logger.info("Attempting swap.")
    logger.debug("Attempting swap of %s and %s paths with propdirs %s and %s.",
                 path0_type, path1_type, propdir0, propdir1)

    # First create path0_new, then path1_new.
    # We will immediately reject the move if path0_new is not valid.
    # -----------------------------------------------------------------------
    # Create new path for [i+-]
    maxlen0 = settings['ensemble'][idx]['tis']['maxlength']
    engine0 = ensembles[idx]['engine']
    # We distinguish between the part of the path that is copied (overlap,cut),
    # and the part that is propagated (prop).
    # First we cut:
    maxlen_cut0 = maxlen0 if inc_sh else maxlen0 - 1
    path_overlap0_new = path_ensemble0.last_path.empty_path(maxlen=maxlen_cut0)
    path0_new_shootpoint = \
        cut_overlap_phasepoints(path_ensemble1, propdir1, path_overlap0_new,
                                side='right', include_shootpoint=inc_sh)
    # Then create empty path to be filled with propagated points:
    # + 1 because the shootpoint will not be used in the propagated path.
    # (Thus, one point of prop_path will be discarded)
    path0_new = None
    path1_new = None
    maxlen_prop0_new = maxlen0 - path_overlap0_new.length + 1
    path_prop0_new = path_ensemble0.last_path.empty_path(
        maxlen=maxlen_prop0_new)
    # Set the system to the shootpoint:
    ensembles[idx]['system'] = path0_new_shootpoint
    logger.debug("Initial point is %s", ensembles[idx]['system'])
    # Propagate the path:
    engine0.propagate(path_prop0_new, ensembles[idx],
                      reverse=propdir1 == 'backwards')
    if inc_sh:  # We need to remove the shootpoint from the propagated path
        path_prop0_new.phasepoints.pop(0)
    # Append the propagated path to the overlap path:
    if propdir1 == 'backwards':
        path0_new = paste_paths(path_prop0_new, path_overlap0_new,
                                overlap=False, maxlen=maxlen0)
    elif propdir1 == 'forwards':
        path0_new = paste_paths(path_overlap0_new, path_prop0_new,
                                overlap=False, maxlen=maxlen0)
    # Check if the new path is valid:
    path0_new.status = set_swap_status(path0_new, propdir1, maxlen0)

    # If the new path is not valid, we stop and reject the entire swap move
    if path0_new.status != 'ACC':
        return reject_repptis_swap_half(path_ensemble0, path_ensemble1,
                                        path0_new, cycle, settings, idx,
                                        ensembles)
    # -----------------------------------------------------------------------
    # Create new path for [(i+1)+-]
    maxlen1 = settings['ensemble'][idx+1]['tis']['maxlength']
    engine1 = ensembles[idx+1]['engine']
    # First we cut:
    maxlen_cut1 = maxlen1 if inc_sh else maxlen1 - 1
    path_overlap1_new = path_ensemble1.last_path.empty_path(maxlen=maxlen_cut1)
    path1_new_shootpoint = \
        cut_overlap_phasepoints(path_ensemble0, propdir0, path_overlap1_new,
                                side='left', include_shootpoint=inc_sh)
    # Then create empty path to be filled with propagated points:
    # + 1 because the shootpoint will not be used in the propagated path.
    # (Thus, one point of prop_path will be discarded)
    maxlen_prop1_new = maxlen1 - path_overlap1_new.length + 1
    path_prop1_new = path_ensemble1.last_path.empty_path(
        maxlen=maxlen_prop1_new)
    # Set the system to the shootpoint:
    ensembles[idx+1]['system'] = path1_new_shootpoint
    logger.debug("Initial point is %s", ensembles[idx+1]['system'])
    # Propagate the path:
    engine1.propagate(path_prop1_new, ensembles[idx+1],
                      reverse=propdir0 == 'backwards')
    if inc_sh:  # We need to remove the shootpoint from the propagated path
        path_prop1_new.phasepoints.pop(0)
    # Append the propagated path to the overlap path:
    if propdir0 == 'backwards':
        path1_new = paste_paths(path_prop1_new, path_overlap1_new,
                                overlap=False, maxlen=maxlen1)
    elif propdir0 == 'forwards':
        path1_new = paste_paths(path_overlap1_new, path_prop1_new,
                                overlap=False, maxlen=maxlen1)
    # Check if the new path is valid:
    path1_new.status = set_swap_status(path1_new, propdir0, maxlen1)
    # -----------------------------------------------------------------------
    # Bookkeeping
    input_tuple = (path_ensemble0, path_ensemble1, path0_new, path1_new,
                   settings, idx, cycle, ensembles)
    return repptis_swap_bookkeeping(input_tuple)


def repptis_swap_bookkeeping(input_tuple):
    """Bookkeeping for the repptis swap move.

    If you reached here, then the first half of the swapping move was already
    successful. Here it is checked whether the second half of the move is
    successful or not. Afterward the final bookkeeping is performed.
    The input_tuple contains the parameters below, which are assigned at the
    beginning of the function.

    Parameters
    ----------
    input_tuple: tuple, it contains
        path_ensemble0 : object like :py:class:`.PathEnsemble`
            The path ensemble of the first trial path ([i+-]).
        path_ensemble1 : object like :py:class:`.PathEnsemble`
            The path ensemble of the second trial path ([(i+1)+-]).
        path0_new : object like :py:class:`.PathBase`
            The first trial path ([i+-]).
        path1_new : object like :py:class:`.PathBase`
            The second trial path ([(i+1)+-]).
        settings : dict
            This dict contains the settings for the REPPTIS method.
        idx : int
            The index of the ensemble that will be swapped with the next one.
        cycle : integer
            The current cycle number.
        ensembles : list of dictionaries of objects
            This is a list of the ensembles we are using in the REPPTIS method.

    Returns
    -------
    out[0] : boolean
        Should the paths be accepted or not?
    out[1] : list of object like :py:class:`.PathBase`
        The trial paths.
    out[2] : string
        The status for the trial paths.

    """
    # Get the input:
    path_ensemble0, path_ensemble1, path0_new, path1_new, settings, idx, \
        cycle, ensembles = input_tuple

    status = 'ACC'
    if path0_new.status != 'ACC':
        raise RuntimeError("First half of swap should've been ACC.")
    if path1_new.status != 'ACC':  # This can happen (FTX or BTX)
        status = path1_new.status

    for path_new, path_ensemble, flag in ((path0_new, path_ensemble0, 's+'),
                                          (path1_new, path_ensemble1, 's-')):
        path_new.status = status
        path_new.weight = 1  # NB: no high acceptance allowed for repptis
        path_new.set_move('ld' if path_ensemble.last_path.get_move() == 'ld'
                          else flag)

    accept = path0_new.status == 'ACC' and path1_new.status == 'ACC'
    # Here is where you would want to put 'plot_repptis_swap' for debugging
    # plot_repptis_swap(path_ensemble0, path_ensemble1, path0_new, path1_new,
    #                   idx, cycle, accept)
    if accept:
        path0_new = path_ensemble0.copy_path_to_generate(path0_new)
        path1_new = path_ensemble1.copy_path_to_generate(path1_new)

    for i, path_new, path_ensemble in ((0, path0_new, path_ensemble0),
                                       (1, path1_new, path_ensemble1)):
        path_ensemble.add_path_data(path_new, status, cycle=cycle)
        ens_set = settings['ensemble'][idx+i]
        if cycle % ens_set.get('output', {}).get('restart-file', 1) == 0:
            write_ensemble_restart(ensembles[idx+i], ens_set)
    return accept, (path0_new, path1_new), status


def set_swap_status(path, propdir, maxlen):
    """Set the status of a path after a swap-propagation in REPPTIS.

    Normally, the path extension should automatically result in an acceptable
    path of the adjacent ensemble, but there are two exceptions:
    1) The path is too long (exceeds maxlen): breaks detailed balance.
    2) The path is too short (length < 3): Not a path breaks detailed balance?

    Parameters
    ----------
    path : object like :py:class:`.PathBase`
        The path that was extended by a REPPTIS swap move.
    propdir : string
        The propagation direction of the path. Either 'forwards' or 'backwards'
    maxlen : int
        The maximum length of the path.

    Returns
    -------
    status : string
        The status of the path after the REPPTIS swap move.

    """
    # First we check if it's too long
    status_map = {'backwards': {True: 'BTX', False: 'ACC'},
                  'forwards': {True: 'FTX', False: 'ACC'}}
    status = status_map[propdir][path.length >= maxlen]
    # Then we check if it's too short. This is theoretically not possible,
    # as we already start from an overlap path that is by definition 2
    # phasepoints long (and we add another one with shooting). However, we
    # keep the check, just in case something terrible happened.
    if path.length < 3:
        raise RuntimeError('Pathlength < 3 is NOT possible in REPPTIS swap.')

    return status


def reject_repptis_swap_wrong_prop(path_ensemble0, path_ensemble1, cycle,
                                   settings, idx, ensembles):
    """Reject REPPTIS swap for incompatible propagation directions.

    E.g. Trying to extend an LMR path of [i+-] into a [(i+1)+-] path by
    backwards propagation.

    Parameters
    ----------
    path_ensemble0 : object like :py:class:`.PathEnsemble`
        The path ensemble of the first trial path ([i+-]).
    path_ensemble1 : object like :py:class:`.PathEnsemble`
        The path ensemble of the second trial path ([(i+1)+-]).
    cycle : integer
        The current cycle number.
    settings : dict
        This dict contains the settings for the REPPTIS method.
    idx : int
        The index of the ensemble that will be swapped with the next one.
    ensembles : list of dictionaries of objects
        This is a list of the ensembles we are using in the REPPTIS method.

    Returns
    -------
    out[0] : boolean
        Should the paths be accepted or not?
    out[1] : list of object like :py:class:`.PathBase`
        The trial paths.
    out[2] : string
        The status for the trial paths.

    """
    logger.info("Swap rejected: incompatible propagation directions")
    status = 'SWD'  # swap wrong direction
    accept = False
    # Make shallow copies:
    trial0 = copy.copy(path_ensemble1.last_path)
    trial1 = copy.copy(path_ensemble0.last_path)
    for trial, move in ((trial0, 's+'), (trial1, 's-')):
        trial.set_move('ld' if trial.get_move() == 'ld' else move)
        trial.weight = 1.0
    path_ensemble0.add_path_data(trial0, status, cycle=cycle)
    path_ensemble1.add_path_data(trial1, status, cycle=cycle)
    for i in (0, 1):
        ens_set = settings['ensemble'][idx+i]
        if cycle % ens_set.get('output', {}).get('restart-file', 1) == 0:
            write_ensemble_restart(ensembles[idx+i], ens_set)
    return accept, (trial0, trial1), status


def reject_repptis_swap_half(path_ensemble0, path_ensemble1, path0_new, cycle,
                             settings, idx, ensembles):
    """Bookkeeping when the first REPPTIS extension is not acceptable.

    If the first extension of the REPPTIS swap move does not result in an
    acceptable path for the [(i+1)+-] ensemble, we do not waste computational
    time on the second extension. We immediately reject the move.

    Parameters
    ----------
    path_ensemble0 : object like :py:class:`.PathEnsemble`
        The path ensemble of the first trial path ([i+-]).
    path_ensemble1 : object like :py:class:`.PathEnsemble`
        The path ensemble of the second trial path ([(i+1)+-]).
    path0_new : object like :py:class:`.PathBase`
        The first trial path ([i+-]).
    cycle : integer
        The current cycle number.
    settings : dict
        This dict contains the settings for the REPPTIS method.
    idx : int
        The index of the ensemble that will be swapped with the next one.
    ensembles : list of dictionaries of objects
        This is a list of the ensembles we are using in the REPPTIS method.

    Returns
    -------
    out[0] : boolean
        Should the paths be accepted or not?
    out[1] : list of object like :py:class:`.PathBase`
        The trial paths.
    out[2] : string
        The status for the trial paths.

    """
    logger.info("Swap rejected: first path not accepted")
    status = 'SWH'  # swap rejected
    accept = False
    # Make shallow copies:
    trial0 = copy.copy(path0_new)
    trial1 = copy.copy(path_ensemble0.last_path)
    for trial, move in ((trial0, 's+'), (trial1, 's-')):
        trial.set_move('ld' if trial.get_move() == 'ld' else move)
        trial.weight = 1.0
    path_ensemble0.add_path_data(trial0, status, cycle=cycle)
    path_ensemble1.add_path_data(trial1, status, cycle=cycle)
    for i in (0, 1):
        ens_set = settings['ensemble'][idx+i]
        if cycle % ens_set.get('output', {}).get('restart-file', 1) == 0:
            write_ensemble_restart(ensembles[idx+i], ens_set)
    return accept, (trial0, trial1), status


def cut_overlap_phasepoints(path_ensemble, propdir, overlap_path, side,
                            include_shootpoint=True):
    """Cut path in overlapping region of adjacent ensembles.

    Cuts the phasepoints in the overlapping region of two adjacent PPTIS
    ensembles [i^+] and [(i+1)^-]. The overlapping region is [L,M] if
    we are cutting from an [(i+1)^+-] path (denoted by side = 'right'), and
    [M,R] if we are cutting from an [i^+-] path (denoted by side = 'left).
    The propagation direction is the direction in which the cut path is to be
    propagated.

    Parameters
    ----------
    path_ensemble: object like :py:class:`.PathEnsemble`
          This is used for storing results for the simulation. It
          is also used for defining the interfaces for this simulation.
    propdir : string
        This is either 'forwards' or 'backwards', and denotes the direction
        in which the cut path is to be propagated.
    overlap_path : object like :py:class:`.Path`
        The path that is to be cut.
    side : string
        This is either 'left' or 'right', and denotes the side from which the
        path originates. left: [i^+-] path to be propagated into [(i+1)^+-]
        right: [(i+1)^+-] path to be propagated into [i^+-]
    include_shootpoint : bool
        If True, the shootpoint of the path is included in the cut path.

    Returns
    -------
    shoot_point : object like :py:class:`.PhasePoint`
        The phasepoint from which the cut path is to be propagated.

    """
    # intf_M defines when to stop cutting phasepoints
    # For now, this is equivalent for all PPTIS ensembles, but might change
    # in the future.
    if side not in ['right', 'left']:
        raise ValueError(f'Unknown side "{side}"')

    intf_m = path_ensemble.interfaces[0] if path_ensemble.ensemble_number == 1\
        else path_ensemble.interfaces[1]

    stop_idx = 1  # to not include shootpoint in the cut path loop

    # Different behavior depending on the side from which the path originates,
    # and the propagation direction. Note that, for the paths. Note that, if
    # the propagation direction is forwards, we reverse the order in which we
    # cut the phasepoints. We do not re-reverse this in the end, because the
    # function 'paste_paths' will be used. As this cut-out path is the backward
    # part of the merged paths, it will be reversed in that function.
    # For the backwards propagation, the cut-out path is the forwards part of
    # the merged paths, and will not be reversed in 'paste_paths' (or here).
    if side == 'right':
        if propdir == 'backwards':
            shoot_point = path_ensemble.last_path.phasepoints[0].copy()
            if include_shootpoint:
                overlap_path.append(shoot_point.copy())
            phase_points = path_ensemble.last_path.phasepoints[stop_idx:]
            for phase_point in phase_points:
                overlap_path.append(phase_point.copy())
                if phase_point.order[0] >= intf_m:  # if>=M
                    break
        elif propdir == 'forwards':
            shoot_point = path_ensemble.last_path.phasepoints[-1].copy()
            if include_shootpoint:
                overlap_path.append(shoot_point.copy())
            phase_points = path_ensemble.last_path.phasepoints[:-stop_idx]
            for phase_point in reversed(phase_points):
                overlap_path.append(phase_point.copy())
                if phase_point.order[0] >= intf_m:  # if>=M
                    break

    elif side == 'left':
        if propdir == 'backwards':
            shoot_point = path_ensemble.last_path.phasepoints[0].copy()
            if include_shootpoint:
                overlap_path.append(shoot_point.copy())
            phase_points = path_ensemble.last_path.phasepoints[stop_idx:]
            for phase_point in phase_points:
                overlap_path.append(phase_point.copy())
                if phase_point.order[0] <= intf_m:  # if<=M
                    break
        elif propdir == 'forwards':
            shoot_point = path_ensemble.last_path.phasepoints[-1].copy()
            if include_shootpoint:
                overlap_path.append(shoot_point.copy())
            phase_points = path_ensemble.last_path.phasepoints[:-stop_idx]
            for phase_point in reversed(phase_points):
                overlap_path.append(phase_point.copy())
                if phase_point.order[0] <= intf_m:  # if<=M
                    break

    return shoot_point


def fix_directions(ensembles, idx):  # pragma: no cover
    """Propose acceptable propagation direcs for [i^+-] and [(i+1)+-] paths.

    If no acceptable directions are found, a rejection is returned.
    Note that this will **break detailed balance**.

    Parameters
    ----------
    ensembles : list of dictionaries of objects
        This is a list of the ensembles we are using in the RETIS method.
        It contains:

        * `path_ensemble`: object like :py:class:`.PathEnsemble`
          This is used for storing results for the simulation. It
          is also used for defining the interfaces for this simulation.
        * `system`: object like :py:class:`.System`
          System is used here since we need access to the temperature
          and to the particle list.
        * `order_function`: object like :py:class:`.OrderParameter`
          The class used for calculating the order parameter(s).
        * `engine`: object like :py:class:`.EngineBase`
          The engine to use for propagating a path.

    idx : integer
        Index of the ensemble, which will be swapped with the (idx+1) ensemble.

    Returns
    -------
    accepted : bool
        Whether the proposed directions can result in an acceptable swap.
    propdir0 : string
        The proposed propagation direction for the [i^+] path. This is either
        'forwards' or 'backwards', or None if no acceptable direction was found
    propdir1 : string
        The proposed propagation direction for the [(i+1)^-] path. This is
        either 'forwards' or 'backwards', or None if no acceptable direction
        was found.
    path0_type : string
        The type of the [i^+] path. This is either 'LMR', 'RMR', 'LML' or 'RML'
    path1_type : string
        The type of the [(i+1)^-] path. This is either 'LMR', 'RMR', 'LML' or
        'RML'.

    Notes
    -----
    This function does **not** obey detailed balance (DB). Used for debugging.
    It could allow detailed balance, but then the path ensemble definitions
    should be changed (weight dependent on path type!)

    """
    path_ensemble0 = ensembles[idx]['path_ensemble']
    path_ensemble1 = ensembles[idx+1]['path_ensemble']

    # Get the start and end points of the paths (L or R)
    path0_old_start, path0_old_end, _, _ = \
        path_ensemble0.last_path.check_interfaces(path_ensemble0.interfaces)
    path1_old_start, path1_old_end, _, _ = \
        path_ensemble1.last_path.check_interfaces(path_ensemble1.interfaces)

    # Define the path types, for nice output
    path0_type = ""+path0_old_start+"M"+path0_old_end
    path1_type = ""+path1_old_start+"M"+path1_old_end

    # Define the allowed propagation directions
    propdirs = {'leftright': {'LMR': 'forwards',
                              'RMR': 'forwards',  # backwards could work too
                              'LML': None,
                              'RML': 'backwards'},
                'rightleft': {'LMR': 'backwards',
                              'RMR': None,
                              'LML': 'backwards',  # forwards could work too
                              'RML': 'forwards'}}

    # allocate directions
    propdir0 = propdirs['leftright'][path0_type]
    propdir1 = propdirs['rightleft'][path1_type]

    # Check whether it is acceptable
    accepted = propdir0 in ('forwards', 'backwards') and \
        propdir1 in ('forwards', 'backwards')

    # Log the move
    msg = "Proposed swap move:\n"
    msg += f"{path_ensemble0.ensemble_name} <-->" + \
           f"{path_ensemble1.ensemble_name}\n" + \
           f"{path0_type} <-----> {path1_type}\n" + \
           f"{propdir0}<>{propdir1}"
    logger.info(msg)

    return accepted, propdir0, propdir1, path0_type, path1_type


def re_propagation_directions(ensembles, idx):
    """Proposes propagation directions for the [i^+-] and [(i+1)^+-] paths.

    And checks whether this results in an acceptable swapping move.

    Parameters
    ----------
    ensembles : list of dictionaries of objects
        This is a list of the ensembles we are using in the RETIS method.
        It contains:

        * `path_ensemble`: object like :py:class:`.PathEnsemble`
          This is used for storing results for the simulation. It
          is also used for defining the interfaces for this simulation.
        * `system`: object like :py:class:`.System`
          System is used here since we need access to the temperature
          and to the particle list.
        * `order_function`: object like :py:class:`.OrderParameter`
          The class used for calculating the order parameter(s).
        * `engine`: object like :py:class:`.EngineBase`
          The engine to use for propagating a path.

    idx : integer
        Index of the ensemble, which will be swapped with the (idx+1) ensemble.

    Returns
    -------
    allowed : bool
        Whether the proposed directions can result in an acceptable swap.
    propdir0 : string
        The proposed propagation direction for the [i^+-] path. This is either
        'forwards' or 'backwards', or None if no acceptable direction was found
    propdir1 : string
        The proposed propagation direction for the [(i+1)^+-] path. This is
        either 'forwards' or 'backwards', or None if no acceptable direction
        was found.
    path0_type : string
        The type of the [i^+-] path. This is 'LMR', 'RMR', 'LML' or 'RML'.
    path1_type : string
        The type of the [(i+1)^+-] path. This is either 'LMR', 'RMR', 'LML' or
        'RML'.

    """
    # Set path ensembles
    path_ensemble0 = ensembles[idx]['path_ensemble']
    path_ensemble1 = ensembles[idx+1]['path_ensemble']

    # Check the start and end points of the paths (L or R)
    path0_old_start, path0_old_end, _, _ = \
        path_ensemble0.last_path.check_interfaces(path_ensemble0.interfaces)
    path1_old_start, path1_old_end, _, _ = \
        path_ensemble1.last_path.check_interfaces(path_ensemble1.interfaces)

    # Define the path types, for nice output
    path0_type = ""+path0_old_start+"M"+path0_old_end
    path1_type = ""+path1_old_start+"M"+path1_old_end

    # Pick two random numbers rand0 and rand1 between [0,1].
    # If rand0 < 0.5: attempt backward propagation for path0. Otherwise forward
    # If rand1 < 0.5: attempt backward propagation for path1. Otherwise forward

    # For now, we will just pick a random number for each path from the
    # random number generator of the 0 ensemble, because there are seed issues.
    idx1 = int(round((ensembles[0]['rgen'].rand())[0]) + 0.1)
    idx2 = int(round((ensembles[0]['rgen'].rand())[0]) + 0.1)
    propdir0 = ['backwards', 'forwards'][idx1]
    propdir1 = ['backwards', 'forwards'][idx2]

    # Check if the proposed move is allowed
    allow = {'leftright': {'L': False, 'R': True},
             'rightleft': {'L': True, 'R': False}}

    allowed0 = allow['leftright'][path0_old_start if propdir0 == 'backwards'
                                  else path0_old_end]
    allowed1 = allow['rightleft'][path1_old_start if propdir1 == 'backwards'
                                  else path1_old_end]

    # Both paths must be allowed to begin the swapping move
    allowed = allowed0 and allowed1

    # Log the move
    msg = "Proposed swap move:\n"
    msg += f"{path_ensemble0.ensemble_name} <-->" + \
           f"{path_ensemble1.ensemble_name}\n" + \
           f"{path0_type} <-----> {path1_type}\n" + \
           f"{propdir0}<>{propdir1}\n" + \
           f"Allowed: {allowed}"
    logger.info(msg)

    return allowed, propdir0, propdir1, path0_type, path1_type


def high_acc_swap(paths, rgen, intf0, intf1, ens_moves):
    """Accept or Reject a swap move using the High Acceptance weights.

    Parameters
    ----------
    paths: list of object like :py:class:`.PathBase`
        The path in the LOWER and UPPER ensemble to exchange.
    rgen : object like :py:class:`.RandomGenerator`
        This is a random generator.
    intf0: list of float
        The interfaces of the LOWER ensemble.
    intf1: list of float
        The interfaces of the HIGHER ensemble.
    ens_moves: list of string
        The moves used in the two ensembles.

    Returns
    -------
    out[0] : boolean
        True if the move should be accepted.

    Notes
    -----
     -  This function is needed only when paths generated via Wire Fencing or
        Stone Skipping are involved.
      - In the case that a path bears a flag 'ld', the swap is accepted,
            but the flag will be unchanged.

    """
    # Crossing before the move
    c1_old = compute_weight(paths[0], intf0, ens_moves[0])
    c2_old = compute_weight(paths[1], intf1, ens_moves[1])
    # Crossing if the move would be accepted
    c1_new = compute_weight(paths[1], intf0, ens_moves[0])
    c2_new = compute_weight(paths[0], intf1, ens_moves[1])
    if c1_old == 0 or c2_old == 0:
        logger.warning("div_by_zero. c1_old, c2_old, ens_moves: [%i,%i], %s",
                       c1_old, c2_old, str(ens_moves))
        p_swap_acc = 1
    else:
        p_swap_acc = c1_new*c2_new/(c1_old*c2_old)

    # Finally, randomly decide to accept or not:
    if rgen.rand() < p_swap_acc:
        return True, 'ACC'  # Accepted

    return False, 'HAS'  # Rejected
