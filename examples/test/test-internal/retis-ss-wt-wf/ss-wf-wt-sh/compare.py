# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Result comparison for RETIS with SS, WT, and WF.

Here we compare a RETIS simulation of 250 steps to known results.
"""
from collections import OrderedDict
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
    compare_data_by_columns,
    compare_numerical_data
)
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

RESULTS = '../results/ss-wf-wt-sh'


def check_path_file(ens):
    """Check that the accepted paths seem ok.

    Parameters
    ----------
    ens : object
        The path ensemble to check.

    Returns
    -------
    status : int
        0 if successful, 1 otherwise.
    """
    logger.info(f'\nReading for {ens.ensemble_name}')
    filename = os.path.join(generate_ensemble_name(ens.ensemble_number),
                            'pathensemble.txt')
    logger.info(f'Reading: {filename}')
    start = ens.start_condition
    end = ('R') if ens.ensemble_number == 0 else ('R', 'L')
    something_weird = False
    with open(filename, 'r', encoding='utf-8') as inputfile:
        for lines in inputfile:
            if lines.startswith('#'):
                continue
            splitline = lines.strip().split()
            status = splitline[7]
            if status != 'ACC':
                continue
            step = int(splitline[0])
            left = splitline[3]
            middle = splitline[4]
            right = splitline[5]
            length = int(splitline[6])
            mino = float(splitline[9])
            maxo = float(splitline[10])

            if length < 3:
                logger.error(f'Suspicious length for path {step}')
                something_weird = True
            if start != left:
                logger.error(
                    f'Inconsistent start: {start} != {left} (step {step})')
                something_weird = True
            if middle != 'M':
                logger.error(f'Middle differ: M != {middle} (step {step})')
                something_weird = True
            if right not in end:
                logger.error(f'Inconsistent end: {right} (step {step})')
                something_weird = True
            cross = [mino < interpos < maxo for interpos in ens.interfaces]
            if ens.ensemble_number == 0:
                idx1, idx2 = 1, 2
            else:
                idx1, idx2 = 0, 1
            if not cross[idx1] or not cross[idx2]:
                something_weird = True
                logger.error(
                    f'Inconsistent crossings: {cross[idx1]} '
                    f'{cross[idx2]} (step {step})')
    if not something_weird:
        logger.info('Accepted paths are OK!')
        return 0
    return 1


def run_check_path_file(settings):
    """Check paths using given simulation settings."""
    ensembles = create_ensembles(settings)
    retval = [check_path_file(ens['path_ensemble']) for ens in ensembles]
    return sum(retval)


def read_path_file(ens):
    """Read information about paths from pathensemble.txt."""
    logger.info(f'\nReading for {ens.ensemble_name}')
    filename = os.path.join(generate_ensemble_name(ens.ensemble_number),
                            'pathensemble.txt')
    logger.info(f'Reading: {filename}')
    paths = OrderedDict()
    path_acc = OrderedDict()
    current_acc = None
    with open(filename, 'r', encoding='utf-8') as inputfile:
        for lines in inputfile:
            if lines.startswith('#'):
                continue
            splitline = lines.strip().split()
            step = int(splitline[0])
            status = splitline[7]
            move = splitline[8]
            paths[step] = {'status': status,
                           'move': move,
                           'parent': current_acc,
                           'swap-parent': (None, None)}
            if status == 'ACC':
                current_acc = step
            path_acc[step] = current_acc
    return paths, path_acc


def get_swap_parent(paths, ensl, ensr, accl, accr):
    """Get the swapping parent for paths."""
    for key in paths:
        idx = key - 1
        if paths[key]['move'] == 's+':
            paths[key]['swap-parent'] = (ensr, accr[idx])
        elif paths[key]['move'] == 's-':
            paths[key]['swap-parent'] = (ensl, accl[idx])


def check_order_swap(data0, data1, special=2):
    """Check that order parameters are consistent for swapping."""
    if special == 0:
        # for [0-] generated from [0+]:
        # - two last of 0- should equal two first in 0+
        ok0 = np.isclose(data0[-2][1], data1[0][1])
        ok1 = np.isclose(data0[-1][1], data1[1][1])
        return ok0 and ok1
    if special == 1:
        # for [0+] generated from [0-]:
        # - two first of 0+ should equal two last from 0-
        ok0 = np.isclose(data0[0][1], data1[-2][1])
        ok1 = np.isclose(data0[1][1], data1[-1][1])
        return ok0 and ok1
    return np.allclose(data0, data1)


def get_index(traj):
    """Just return index from a comment."""
    return int(traj['comment'][0].split('Cycle:')[1].split(',')[0])


def check_swaps(paths, accepted, ens, kind):
    """Check accepted swaps."""
    ofile0 = OrderPathFile(
        os.path.join(generate_ensemble_name(ens), 'order.txt'), 'r'
    )
    traj0 = ofile0.load()
    if kind == 'left':
        special = 1 if ens == 1 else 2
        ens2 = ens - 1
        move = 's-'
    else:
        special = 0 if ens == 0 else 2
        ens2 = ens + 1
        move = 's+'
    ofile1 = OrderPathFile(
        os.path.join(generate_ensemble_name(ens2), 'order.txt'), 'r'
    )
    traj1 = ofile1.load()
    traj1i, idx1 = {}, None
    errors = set()
    everything_is_ok = True
    for traj0i in traj0:
        idx0 = get_index(traj0i)
        if idx0 in accepted and paths[idx0]['move'] == move:
            parent = paths[idx0]['swap-parent']
            swap_ok = False
            found = False
            if traj1i and idx1 == parent[1]:
                found = True
                swap_ok = check_order_swap(traj0i['data'],
                                           traj1i['data'],
                                           special=special)
            if not found:
                for traj1i in traj1:
                    idx1 = get_index(traj1i)
                    if idx1 == parent[1]:
                        found = True
                        swap_ok = check_order_swap(traj0i['data'],
                                                   traj1i['data'],
                                                   special=special)
                        break
            if not found:
                logger.warning(f'Could not find parent for {idx0}')
                everything_is_ok = False
            elif not swap_ok:
                logger.error(f'Comparison failed for {idx0}')
                everything_is_ok = False
                errors.add(idx0)
    if everything_is_ok:
        logger.info('All swaps are ok!')
        return 0
    logger.error(f'Error for some swaps: {errors}')
    return 1


def check_ensemble_swaps(settings):
    """Check swaps for ensembles from settings."""
    ensembles = create_ensembles(settings)
    path_info = {}
    path_acc = {}
    names = []
    for i, ens in enumerate(ensembles):
        pathi, patha = read_path_file(ens['path_ensemble'])
        path_info[i] = pathi
        path_acc[i] = patha
        names.append(ens['path_ensemble'].ensemble_name)
    status = 0
    for i in range(len(ensembles)):
        if i == 0:
            get_swap_parent(path_info[i], None, i+1, None, path_acc[i+1])
            logger.info(f'\nChecking {names[i]} <- {names[i+1]} swaps...')
            status += check_swaps(path_info[i], path_acc[i], i, kind='right')
        elif i == len(ensembles) - 1:
            get_swap_parent(path_info[i], i-1, None, path_acc[i-1], None)
            logger.info(f'\nChecking {names[i]} <- {names[i-1]} swaps...')
            status += check_swaps(path_info[i], path_acc[i], i, kind='left')
        else:
            get_swap_parent(path_info[i], i-1, i+1,
                            path_acc[i-1], path_acc[i+1])
            logger.info(f'\nChecking {names[i]} -> {names[i+1]} swaps...')
            status += check_swaps(path_info[i], path_acc[i], i, kind='right')
            logger.info(f'Checking {names[i]} <- {names[i-1]} swaps...')
            status += check_swaps(path_info[i], path_acc[i], i, kind='left')
    return status


def compare_ens_files(settings, fname, ftype=None):
    """Compare ensemble text files using centralized logic."""
    inter = settings['simulation']['interfaces']
    retval = 0
    for i in range(len(inter)):
        ens_dir = generate_ensemble_name(i)
        fil1 = os.path.join(ens_dir, fname)
        fil2 = os.path.join(RESULTS, ens_dir, fname)
        if fname == 'pathensemble.txt':
            equal, msg = compare_path_ensemble_data(fil1, fil2)
        elif ftype == 'energy':
            equal, msg = compare_data_by_columns(fil1, fil2, ftype)
        else:
            equal, msg = compare_numerical_data(fil1, fil2)
        if not equal:
            logger.error(f'Mismatch in {fil1}: {msg}')
            retval += 1
        else:
            logger.info(f'Files are equal: {fil1}')
    return retval


def main():
    """Run the full comparison."""
    sets = parse_settings_file('retis.rst')
    logger.info('\nComparing pathensemble.txt files')
    logger.info('================================')
    ret1 = compare_ens_files(sets, 'pathensemble.txt')
    logger.info('\nCheck swaps')
    logger.info('===========')
    ret2 = check_ensemble_swaps(sets)
    logger.info('\nCheck accepted paths')
    logger.info('====================')
    ret3 = run_check_path_file(sets)
    logger.info('\nComparing energy.txt files')
    logger.info('================================')
    ret4 = compare_ens_files(sets, 'energy.txt', 'energy')
    logger.info('\nComparing order.txt files')
    logger.info('================================')
    ret5 = compare_ens_files(sets, 'order.txt', 'order')

    retval = ret1 + ret2 + ret3 + ret4 + ret5
    if retval == 0:
        logger.info('\nComparison is successful!')
    else:
        logger.error('\nComparison failed!')
    return retval


if __name__ == '__main__':
    colorama.init(autoreset=True)
    sys.exit(main())
