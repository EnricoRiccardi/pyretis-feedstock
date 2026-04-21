# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Result comparison for the explorer simulation.

We compare an explore simulation of 16 steps — with initial paths loaded
from sparse trajectory files — against stored reference results.

Determinism is guaranteed by:
  - Fixed seeds in retis.rst (engine seed=0, TIS seed=0, velocity seed=0,
    rgen = rgen-borg at every level)
  - Initial paths loaded from committed files in pippo/

Three checks are performed:
  1. pathensemble.txt, order.txt and energy.txt match stored references.
  2. Every accepted path satisfies its ensemble start/end/interface
     conditions (self-consistency, independent of reference data).
  3. Swapped paths carry consistent order-parameter data across the
     two ensembles that participated in the swap.
"""
import logging
import os
import sys
import colorama
import numpy as np
from pyretis.inout.settings import parse_settings_file
from pyretis.core.pathensemble import generate_ensemble_name
from pyretis.setup.createsimulation import create_ensembles
from pyretis.inout.formats.order import OrderPathFile
from pyretis.testing.simulation_comparison import (
    compare_path_ensemble_data,
    compare_numerical_data,
)

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

RESULTS = 'results'
INPUT_FILE = 'retis.rst'


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _section(msg):
    logger.info('\n%s', msg)
    logger.info('=' * len(msg))


def _ensemble_list(settings):
    """Return (ens_obj, ens_dir) pairs for all ensembles in the run."""
    return [
        (e['path_ensemble'],
         generate_ensemble_name(e['path_ensemble'].ensemble_number))
        for e in create_ensembles(settings)
    ]


# ---------------------------------------------------------------------------
# 1. File-level comparison against stored reference data
# ---------------------------------------------------------------------------

def compare_files(settings):
    """Compare pathensemble.txt, order.txt and energy.txt for every ensemble.

    Parameters
    ----------
    settings : dict
        Parsed simulation settings.

    Returns
    -------
    status : int
        0 if all comparisons pass, positive integer otherwise.
    """
    files = ('pathensemble.txt', 'order.txt', 'energy.txt')
    status = 0
    for _, ens_dir in _ensemble_list(settings):
        _section(f'Ensemble {ens_dir}')
        for fname in files:
            ref = os.path.join(RESULTS, ens_dir, fname)
            new = os.path.join(ens_dir, fname)
            logger.info('  Comparing %s ...', fname)
            if fname == 'pathensemble.txt':
                equal, msg = compare_path_ensemble_data(new, ref)
            else:
                equal, msg = compare_numerical_data(new, ref)
            if equal:
                logger.info('    -> OK')
            else:
                logger.error('    -> MISMATCH: %s', msg)
                status += 1
    return status


# ---------------------------------------------------------------------------
# 2. Self-consistency: accepted paths must satisfy ensemble requirements
# ---------------------------------------------------------------------------

def check_path_file(ens, ens_dir):
    """Verify that accepted paths meet the ensemble interface conditions.

    Parameters
    ----------
    ens : object
        PathEnsemble object.
    ens_dir : str
        Directory name for this ensemble (e.g. '001').

    Returns
    -------
    status : int
        0 if all checks pass, 1 otherwise.
    """
    filename = os.path.join(ens_dir, 'pathensemble.txt')
    start = ens.start_condition
    valid_ends = ('R', 'L')
    broken = False

    with open(filename, 'r', encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.strip().split()
            if cols[7] != 'ACC':
                continue
            step = int(cols[0])
            left, middle, right = cols[3], cols[4], cols[5]
            length = int(cols[6])
            mino, maxo = float(cols[9]), float(cols[10])

            if length < 3:
                logger.error('Suspicious length at step %d', step)
                broken = True
            if left != start:
                logger.error('Wrong start at step %d: %s != %s',
                             step, left, start)
                broken = True
            if middle != 'M':
                logger.error('Middle != M at step %d: %s', step, middle)
                broken = True
            if right not in valid_ends:
                logger.error('Wrong end at step %d: %s', step, right)
                broken = True
            cross = [mino < iface < maxo for iface in ens.interfaces]
            if not cross[0] or not cross[1]:
                logger.error(
                    'Interface crossing check failed at step %d '
                    '(cross=%s, mino=%.4f, maxo=%.4f)',
                    step, cross, mino, maxo,
                )
                broken = True

    if not broken:
        logger.info('  Accepted paths OK for %s', ens.ensemble_name)
    return int(broken)


def run_check_path_files(settings):
    """Run check_path_file for all ensembles."""
    return sum(
        check_path_file(ens, ens_dir)
        for ens, ens_dir in _ensemble_list(settings)
    )


# ---------------------------------------------------------------------------
# 3. Swap consistency: swapped paths must share order-parameter data
# ---------------------------------------------------------------------------

def _get_cycle(traj):
    """Extract the cycle index from a trajectory comment."""
    return int(traj['comment'][0].split('Cycle:')[1].split(',')[0])


def _order_swap_ok(data0, data1, special):
    """Check that two order-parameter arrays are compatible after a swap.

    special == 0  : [0-] was generated from [0+]
    special == 1  : [0+] was generated from [0-]
    otherwise     : full array match
    """
    if special == 0:
        return (np.isclose(data0[-2][1], data1[0][1]) and
                np.isclose(data0[-1][1], data1[1][1]))
    if special == 1:
        return (np.isclose(data0[0][1], data1[-2][1]) and
                np.isclose(data0[1][1], data1[-1][1]))
    return bool(np.allclose(data0, data1))


def _read_path_table(ens_dir):
    """Parse pathensemble.txt; return (paths dict, accepted-at dict)."""
    filename = os.path.join(ens_dir, 'pathensemble.txt')
    paths = {}
    path_acc = {}
    current_acc = None
    with open(filename, 'r', encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.strip().split()
            step = int(cols[0])
            status, move = cols[7], cols[8]
            paths[step] = {
                'status': status,
                'move': move,
                'swap-parent': (None, None),
            }
            if status == 'ACC':
                current_acc = step
            path_acc[step] = current_acc
    return paths, path_acc


def _annotate_swap_parents(paths, ens_left_acc, ens_right_acc):
    """Fill in swap-parent information for each path."""
    for key in paths:
        idx = key - 1
        move = paths[key]['move']
        if move == 's+' and ens_right_acc is not None:
            paths[key]['swap-parent'] = (None, ens_right_acc.get(idx))
        elif move == 's-' and ens_left_acc is not None:
            paths[key]['swap-parent'] = (None, ens_left_acc.get(idx))


def _check_swaps_one_direction(paths, path_acc, ens_dir, partner_dir,
                               move_symbol, special):
    """Verify order-parameter consistency for swaps in one direction."""
    ofile0 = OrderPathFile(os.path.join(ens_dir, 'order.txt'), 'r')
    traj0_list = list(ofile0.load())

    ofile1 = OrderPathFile(os.path.join(partner_dir, 'order.txt'), 'r')
    traj1_by_idx = {_get_cycle(t): t for t in ofile1.load()}

    accepted = {s for s, a in path_acc.items() if a == s}
    errors = set()
    ok = True

    for traj0 in traj0_list:
        idx0 = _get_cycle(traj0)
        if idx0 not in accepted or paths[idx0]['move'] != move_symbol:
            continue
        parent_idx = paths[idx0]['swap-parent'][1]
        traj1 = traj1_by_idx.get(parent_idx)
        if traj1 is None:
            logger.warning('Parent %d not found for swap at step %d',
                           parent_idx, idx0)
            ok = False
            continue
        if not _order_swap_ok(traj0['data'], traj1['data'], special):
            logger.error('Order-parameter mismatch for swap at step %d', idx0)
            errors.add(idx0)
            ok = False

    if ok:
        logger.info('  Swap check %s(%s) <-> (%s) OK',
                    move_symbol, ens_dir, partner_dir)
    else:
        logger.error('  Swap errors at steps: %s', errors)
    return 0 if ok else 1


def check_ensemble_swaps(settings):
    """Check swap consistency for all adjacent ensemble pairs."""
    pairs = _ensemble_list(settings)
    n = len(pairs)

    path_info = {}
    path_acc = {}
    dirs = [ens_dir for _, ens_dir in pairs]
    for i, (_, ens_dir) in enumerate(pairs):
        path_info[i], path_acc[i] = _read_path_table(ens_dir)

    for i in range(n):
        left_acc = path_acc[i - 1] if i > 0 else None
        right_acc = path_acc[i + 1] if i < n - 1 else None
        _annotate_swap_parents(path_info[i], left_acc, right_acc)

    status = 0
    for i in range(n):
        if i < n - 1:
            status += _check_swaps_one_direction(
                path_info[i], path_acc[i], dirs[i], dirs[i + 1], 's+', 2
            )
        if i > 0:
            status += _check_swaps_one_direction(
                path_info[i], path_acc[i], dirs[i], dirs[i - 1], 's-', 2
            )
    return status


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    """Run all three checks; return combined exit code."""
    settings = parse_settings_file(INPUT_FILE)

    _section('1. Comparing output files against reference data')
    ret1 = compare_files(settings)

    _section('2. Self-consistency check on accepted paths')
    ret2 = run_check_path_files(settings)

    _section('3. Swap order-parameter consistency')
    ret3 = check_ensemble_swaps(settings)

    total = ret1 + ret2 + ret3
    if total == 0:
        logger.info('\nAll checks passed.')
    else:
        logger.error('\n%d check(s) failed.', total)
    return total


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,
                        format='%(levelname)s: %(message)s')
    colorama.init(autoreset=True)
    sys.exit(main())
