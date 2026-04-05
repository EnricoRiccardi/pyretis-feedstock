# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Report comparison for TIS multiple interfaces.

Here we compare a TIS simulation of 50 steps to known results.
"""
import os
import sys
import colorama
from pyretis.inout import print_to_screen
from pyretis.inout.settings import parse_settings_file
from pyretis.core.pathensemble import generate_ensemble_name
from pyretis.setup.createsimulation import create_ensembles
from pyretis.testing.simulation_comparison import (
    compare_text_line_by_line,
    compare_reports_normalized,
    compare_path_ensemble_data
)


RESULTS = 'results'


def check_path_file(ens):
    """Check that the accepted paths seem ok.

    Parameters
    ----------
    ens : pyretis.core.pathensemble.PathEnsemble
        The path ensemble to check.

    Returns
    -------
    status : int
        0 if successful, 1 otherwise.
    """
    print_to_screen(f'\nReading for {ens.ensemble_name}')
    filename = os.path.join(generate_ensemble_name(ens.ensemble_number),
                            'pathensemble.txt')
    print_to_screen(f'Reading: {filename}')
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
                print_to_screen(f'Suspicious length for path {step}',
                                level='error')
                something_weird = True
            if start != left:
                print_to_screen(
                    f'Inconsistent start: {start} != {left} (step {step})',
                    level='error')
                something_weird = True
            if middle != 'M':
                print_to_screen(
                    f'Middle differ: M != {middle} (step {step})',
                    level='error')
                something_weird = True
            if right not in end:
                print_to_screen(
                    f'Inconsistent end: {right} (step {step})',
                    level='error')
                something_weird = True
            cross = [mino < interpos < maxo for interpos in ens.interfaces]
            if ens.ensemble_number == 0:
                idx1, idx2 = 1, 2
            else:
                idx1, idx2 = 0, 1
            if not cross[idx1] or not cross[idx2]:
                something_weird = True
                print_to_screen(
                    f'Inconsistent crossings: {cross[idx1]} '
                    f'{cross[idx2]} (step {step})',
                    level='error'
                )
    if not something_weird:
        print_to_screen('Accepted paths are OK!', level='success')
        return 0
    return 1


def run_check_path_file(settings):
    """Check paths using given simulation settings.

    Parameters
    ----------
    settings : dict
        The simulation settings.

    Returns
    -------
    status : int
        Combined status of all checks.
    """
    ensembles = create_ensembles(settings)
    retval = [check_path_file(ens['path_ensemble']) for ens in ensembles]
    return sum(retval)


def compare_path_files(settings):
    """Compare pathensemble.txt files.

    Parameters
    ----------
    settings : dict
        The simulation settings.

    Returns
    -------
    status : int
        Combined status of all comparisons.
    """
    inter = settings['simulation']['interfaces']
    retval = 0
    for i in range(1, len(inter)):
        ens_dir = generate_ensemble_name(i)
        fil1 = os.path.join(ens_dir, 'pathensemble.txt')
        fil2 = os.path.join(RESULTS, ens_dir, 'pathensemble.txt')
        equal, msg = compare_path_ensemble_data(fil1, fil2)
        if not equal:
            print_to_screen(f'Mismatch in {fil1}: {msg}', level='error')
            retval += 1
        else:
            print_to_screen(f'Files are equal: {fil1}', level='success')
    return retval


def compare_rst_files(settings):
    """Compare tis.rst files.

    Parameters
    ----------
    settings : dict
        The simulation settings.

    Returns
    -------
    status : int
        Combined status of all comparisons.
    """
    inter = settings['simulation']['interfaces']
    retval = 0
    for i in range(1, len(inter)):
        fil1 = os.path.join(f"tis-00{i}.rst")
        fil2 = os.path.join(RESULTS, f"tis-00{i}.rst")
        # exe_path differs by run directory; skip those lines.
        equal, msg = compare_text_line_by_line(
            fil1, fil2, skip_keys=['exe_path']
        )
        if not equal:
            print_to_screen(f'Mismatch in {fil1}: {msg}', level='error')
            retval += 1
        else:
            print_to_screen(f'Files are equal: {fil1}', level='success')
    return retval


def compare_prob_files():
    """Compare probabilities files.

    Returns
    -------
    status : int
        Combined status of all report comparisons.
    """
    names = ['tis-multiple_report.html',
             'tis-multiple_report.rst',
             'tis-multiple_report.tex']
    retval = 0
    for name in names:
        fil1 = os.path.join('report', name)
        fil2 = os.path.join(RESULTS, name)
        equal, msg = compare_reports_normalized(fil1, fil2)
        if not equal:
            print_to_screen(
                f'Report mismatch {name}: {msg}', level='error'
            )
            retval += 1
        else:
            print_to_screen(f'Reports are equal: {name}', level='success')
    return retval


def main():
    """Run the full comparison."""
    sets = parse_settings_file('tis-multiple.rst')
    print_to_screen('\nCheck crossing probabilities', level='message')
    print_to_screen('============================', level='message')
    ret4 = compare_prob_files()
    print_to_screen('\nComparing tis.rst files', level='message')
    print_to_screen('=======================', level='message')
    ret1 = compare_rst_files(sets)
    print_to_screen('\nComparing pathensemble.txt files', level='message')
    print_to_screen('================================', level='message')
    ret2 = compare_path_files(sets)
    print_to_screen('\nCheck accepted paths', level='message')
    print_to_screen('====================', level='message')
    ret3 = run_check_path_file(sets)
    return ret1 + ret2 + ret3 + ret4


if __name__ == '__main__':
    colorama.init(autoreset=True)
    sys.exit(main())
