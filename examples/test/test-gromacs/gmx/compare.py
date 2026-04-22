# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Shared comparison helpers for GROMACS example simulations."""
import argparse
import logging
import os
import sys

import colorama
import numpy as np

from pyretis.core.pathensemble import generate_ensemble_name
from pyretis.inout.formats.order import OrderPathFile
from pyretis.inout.settings import parse_settings_file
from pyretis.setup.createsimulation import create_ensembles
from pyretis.testing.simulation_comparison import (
    compare_data_by_columns,
    compare_path_ensemble_data,
    compare_simulation_files,
    compare_traj_archive,
)
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


TRAJ = 'traj'
RETIS_RST = 'retis.rst'


def soft_archive_comparison(dir1, dir2):
    """Compare two archives. Here we accept a difference in energies.

    For GROMACS1 and 2, the potential energy output is different.
    This is due to a difference in the way GROMACS adds dispersion
    corrections when continuing a simulation. So, here we expect
    a failure for the energy files. We check these energy files
    manually, by skipping the potential energy terms.

    Parameters
    ----------
    dir1 : string
        The path to the first archive used in the comparison.
    dir2 : string
        The path to the second archive used in the comparison.

    Returns
    -------
    errors : list of tuples
        These are the files which were found to be different.
    """
    errors = []
    archive_diffs = compare_traj_archive(dir1, dir2)
    for (file1, file2) in archive_diffs:
        # Check if it is an energy file:
        if os.path.basename(file1) == 'energy.txt':
            equal, _ = compare_data_by_columns(
                file1, file2, file_type='energy', skip=['vpot']
            )
            if not equal:
                errors.append((file1, file2))
        else:
            # For other files in archive, they must match exactly:
            equal, _ = compare_simulation_files(file1, file2, mode='cmp')
            if not equal:
                errors.append((file1, file2))
    return errors


def compare_ensemble(run1, run2, ensemble,
                     path_skip=None, energy_skip=None, traj_skip=False):
    """Run the comparison for an ensemble.

    Parameters
    ----------
    run1 : string
        The path to the first directory to use for the comparison.
    run2 : string
        The path to the second directory to use for the comparison.
    ensemble : string
        The label for the current ensemble we are investigating, e.g.
        000, 001, or similar.
    path_skip : list of integers, optional
        This selects columns to skip when comparing path ensemble
        data.
    energy_skip : list of strings, optional
        This selects columns to skip in the energy file.
    traj_skip : bool, optional
        If True, we skip comparison of trajectory archives.

    Returns
    -------
    status : bool
        True if the comparison was successful, False otherwise.
    """
    logger.info(f'Comparing for "{ensemble}"')
    for filei in ('energy.txt', 'order.txt', 'pathensemble.txt'):
        logger.info(f'\tComparing {filei} files...')
        file1 = os.path.join(run1, ensemble, filei)
        file2 = os.path.join(run2, ensemble, filei)
        if filei == 'pathensemble.txt':
            equal, msg = compare_path_ensemble_data(
                file1, file2, skip=path_skip
            )
        elif filei == 'energy.txt' and energy_skip:
            equal, msg = compare_data_by_columns(
                file1, file2, file_type='energy', skip=energy_skip
            )
        else:
            equal, msg = compare_simulation_files(
                file1, file2, mode='numerical'
            )
        lvl = 'success' if equal else 'error'
        if lvl == 'error':
            logger.error(f'\t\t-> {msg}')
        else:
            logger.info(f'\t\t-> {msg}')
        if not equal:
            return False

    if not traj_skip:
        logger.info('\tComparing for trajectory archive:')
        # Check folders like 0_traj and 1_traj
        for i in {'0_', '1_'}:
            archive_errors = soft_archive_comparison(
                os.path.join(run1, ensemble, i + TRAJ),
                os.path.join(run2, ensemble, i + TRAJ),
            )
            if archive_errors:
                logger.error('\t\t-> Archives differ')
                return False
        logger.info('\t\t-> Archives are equal')
    else:
        logger.info('\tSkipping comparison for trajectory archive.')

    return True


def _section(msg):
    """Log a short section header."""
    logger.info('\n%s', msg)
    logger.info('=' * len(msg))


def sparse_load_check_output_files(settings):
    """Verify expected sparse-load output files exist and are non-empty."""
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


def _sparse_load_ends(ensemble_number, valid_ends):
    """Return the accepted end states for the current ensemble."""
    if callable(valid_ends):
        return tuple(valid_ends(ensemble_number))
    return tuple(valid_ends)


def sparse_load_check_path_file(ens, valid_ends=('R', 'L')):
    """Verify that accepted sparse-load paths satisfy ensemble conditions."""
    filename = os.path.join(
        generate_ensemble_name(ens.ensemble_number), 'pathensemble.txt'
    )
    start = ens.start_condition
    ends = _sparse_load_ends(ens.ensemble_number, valid_ends)
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
            if right not in ends:
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


def sparse_load_check_path_files(settings, valid_ends=('R', 'L')):
    """Run sparse-load accepted-path checks for every ensemble."""
    ensembles = create_ensembles(settings)
    return sum(
        sparse_load_check_path_file(e['path_ensemble'], valid_ends)
        for e in ensembles
    )


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
    """Fill in swap-parent information for each sparse-load path."""
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
    traj1_by_idx = {_get_cycle(traj): traj for traj in ofile1.load()}

    accepted = {step for step, acc in path_acc.items() if acc == step}
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


def sparse_load_check_ensemble_swaps(settings):
    """Check sparse-load swap consistency for all adjacent ensembles."""
    ensembles = create_ensembles(settings)
    n_ensembles = len(ensembles)

    path_info = {}
    path_acc = {}
    for i in range(n_ensembles):
        path_info[i], path_acc[i] = _read_path_table(i)

    for i in range(n_ensembles):
        left_acc = path_acc[i - 1] if i > 0 else None
        right_acc = path_acc[i + 1] if i < n_ensembles - 1 else None
        _annotate_swap_parents(path_info[i], left_acc, right_acc)

    status = 0
    for i in range(n_ensembles):
        if i < n_ensembles - 1:
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


def main_sparse_load(settings_file='retis-load-rc-run.rst',
                     valid_ends=('R', 'L')):
    """Run the shared sparse-load structural checks."""
    colorama.init(autoreset=True)
    settings = parse_settings_file(settings_file)

    _section('1. Checking output files')
    ret1 = sparse_load_check_output_files(settings)

    _section('2. Self-consistency check on accepted paths')
    ret2 = sparse_load_check_path_files(settings, valid_ends=valid_ends)

    _section('3. Swap order-parameter consistency')
    ret3 = sparse_load_check_ensemble_swaps(settings)

    total = ret1 + ret2 + ret3
    if total == 0:
        logger.info('\nAll checks passed.')
    else:
        logger.error('\n%d check(s) failed.', total)
    return total


def main():
    """Parse arguments and execute the comparison.

    Returns
    -------
    out : int
        0 if the comparison was successful, 1 otherwise.
    """
    colorama.init(autoreset=True)
    parser = argparse.ArgumentParser()
    parser.add_argument('directories', metavar='N', type=str, nargs=2,
                        help='The directories to use for the comparison')
    parser.add_argument('--path_skip', nargs='+', type=int,
                        help='Columns to skip for pathensemble.txt')
    parser.add_argument('--energy_skip', nargs='+', type=str,
                        help='Columns to skip for energy.txt')
    parser.add_argument('--traj_skip', action='store_true',
                        help='Skip comparison for trajectory archives')
    args = parser.parse_args()

    run1 = args.directories[0]
    run2 = args.directories[1]

    settings = parse_settings_file(os.path.join(run1, RETIS_RST))
    inter = settings['simulation']['interfaces']
    for i, _ in enumerate(inter):
        ensemble_dir = generate_ensemble_name(i)
        result = compare_ensemble(
            run1, run2, ensemble_dir,
            path_skip=args.path_skip,
            energy_skip=args.energy_skip,
            traj_skip=args.traj_skip,
        )
        if not result:
            logger.error(f'Comparison failed for {ensemble_dir}. Aborting!')
            return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
