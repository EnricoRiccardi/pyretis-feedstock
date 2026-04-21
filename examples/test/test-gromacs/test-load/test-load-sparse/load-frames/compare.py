# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Structural checks for GROMACS RETIS sparse-load (frames).

Exact file-by-file reference comparisons turned out to be too brittle for
short GROMACS runs across repeated executions, even with fixed seeds. Instead,
we validate invariants that should hold for every successful run:

1. Expected output files exist and are non-empty.
2. Every accepted path satisfies its ensemble start/end/interface conditions.
3. Swapped paths carry consistent order-parameter data across ensembles.
"""
import logging
import os
import sys

import colorama
import numpy as np

from pyretis.core.pathensemble import generate_ensemble_name
from pyretis.inout.formats.order import OrderPathFile
from pyretis.inout.settings import parse_settings_file
from pyretis.setup.createsimulation import create_ensembles

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

INPUT_FILE = 'retis-load-rc-run.rst'


def _section(msg):
    logger.info('\n%s', msg)
    logger.info('=' * len(msg))


def check_output_files(settings):
    """Verify that expected output files exist and are non-empty."""
    files = ('pathensemble.txt', 'order.txt', 'energy.txt')
    status = 0
    for i in range(len(settings['simulation']['interfaces'])):
        ens_dir = generate_ensemble_name(i)
        _section(f'Ensemble {ens_dir}')
        for fname in files:
            path = os.path.join(ens_dir, fname)
            logger.info('  Checking %s ...', fname)
            if not os.path.isfile(path):
                logger.error('    -> MISSING')
                status += 1
                continue
            if os.path.getsize(path) == 0:
                logger.error('    -> EMPTY')
                status += 1
                continue
            logger.info('    -> OK')
    return status


def check_path_file(ens):
    """Verify that accepted paths meet the ensemble interface conditions."""
    filename = os.path.join(
        generate_ensemble_name(ens.ensemble_number), 'pathensemble.txt'
    )
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
            idx1, idx2 = (1, 2) if ens.ensemble_number == 0 else (0, 1)
            cross = [mino < iface < maxo for iface in ens.interfaces]
            if not cross[idx1] or not cross[idx2]:
                logger.error(
                    'Interface crossing check failed at step %d '
                    '(cross=%s, mino=%.6f, maxo=%.6f)',
                    step, cross, mino, maxo,
                )
                broken = True

    if not broken:
        logger.info('  Accepted paths OK for %s', ens.ensemble_name)
    return int(broken)


def run_check_path_files(settings):
    """Run check_path_file for all ensembles."""
    ensembles = create_ensembles(settings)
    return sum(check_path_file(e['path_ensemble']) for e in ensembles)


def _get_cycle(traj):
    """Extract the cycle index from a trajectory comment."""
    return int(traj['comment'][0].split('Cycle:')[1].split(',')[0])


def _order_swap_ok(data0, data1, special):
    """Check that two order-parameter arrays are compatible after a swap."""
    if special == 0:
        return (np.isclose(data0[-2][1], data1[0][1]) and
                np.isclose(data0[-1][1], data1[1][1]))
    if special == 1:
        return (np.isclose(data0[0][1], data1[-2][1]) and
                np.isclose(data0[1][1], data1[-1][1]))
    return bool(np.allclose(data0, data1))


def _read_path_table(ens_number):
    """Parse pathensemble.txt; return (paths dict, accepted-at dict)."""
    filename = os.path.join(
        generate_ensemble_name(ens_number), 'pathensemble.txt'
    )
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


def _check_swaps_one_direction(paths, path_acc, ens_number, partner_number,
                               move_symbol, special):
    """Verify order-parameter consistency for swaps in one direction."""
    ofile0 = OrderPathFile(
        os.path.join(generate_ensemble_name(ens_number), 'order.txt'), 'r'
    )
    traj0_list = list(ofile0.load())

    ofile1 = OrderPathFile(
        os.path.join(generate_ensemble_name(partner_number), 'order.txt'), 'r'
    )
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
        logger.info(
            '  Swap check %s(%d) <-> (%d) OK',
            move_symbol, ens_number, partner_number,
        )
    else:
        logger.error('  Swap errors at steps: %s', errors)
    return 0 if ok else 1


def check_ensemble_swaps(settings):
    """Check swap consistency for all adjacent ensemble pairs."""
    ensembles = create_ensembles(settings)
    n = len(ensembles)

    path_info = {}
    path_acc = {}
    for i in range(n):
        path_info[i], path_acc[i] = _read_path_table(i)

    for i in range(n):
        left_acc = path_acc[i - 1] if i > 0 else None
        right_acc = path_acc[i + 1] if i < n - 1 else None
        _annotate_swap_parents(path_info[i], left_acc, right_acc)

    status = 0
    for i in range(n):
        if i < n - 1:
            special = 0 if i == 0 else 2
            status += _check_swaps_one_direction(
                path_info[i], path_acc[i], i, i + 1, 's+', special
            )
        if i > 0:
            special = 1 if i == 1 else 2
            status += _check_swaps_one_direction(
                path_info[i], path_acc[i], i, i - 1, 's-', special
            )
    return status


def main():
    """Run all checks; return combined exit code."""
    settings = parse_settings_file(INPUT_FILE)

    _section('1. Checking output files')
    ret1 = check_output_files(settings)

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
