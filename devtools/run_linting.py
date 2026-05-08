# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Run linting for a set of source files.

Here, we will also give some over-all evaluation of the results.
The exit code of this script will be != 0 if one or more files
have a score below a pre-set value or if the changes are "too big".
"""
import argparse
import colorama
import io
import pathlib
import re
import sys
from colorama import Fore
from operator import itemgetter
from pathlib import Path
from pylint import lint
from pylint.reporters.text import TextReporter

SCORE_THRESHOLD = 9.0


def look_for_source_files(rootdir, skip=None):
    """Find .py files by walking the given directory."""
    root = pathlib.Path(rootdir)
    glob = sorted(root.rglob('*.py'))
    if skip is None:
        return glob
    skip_dirs = [pathlib.Path(i) for i in skip]
    return [
        i for i in glob if not any(parent in skip_dirs for parent in i.parents)
    ]


def lint_file(filename, reg_score, reg_warning, reg_error, white_list=None):
    """Run linting for a single file."""
    cmd = [str(Path(filename))]
    if white_list is not None:
        cmd.append(f"--ignored-modules={','.join(white_list)}")

    buf = io.StringIO()

    lint.Run(cmd, reporter=TextReporter(buf), exit=False)

    txt = buf.getvalue()

    warnings = reg_warning.findall(txt)
    errors = reg_error.findall(txt)
    result = reg_score.findall(txt)

    score, delta = process_result(result)

    return score, delta, warnings, errors


def run_linting(files, white_list=None):
    """Run pylint for a set of files, individually."""
    reg_score = re.compile(r'-?\d+\.\d+\/\d+')
    reg_warning = re.compile(r'.+warning \(W.+')
    reg_error = re.compile(r'.+error \(E.+')
    results = {
        'all': [],
        'below': [],
        'big-change': [],
        'change': [],
        'warnings': [],
        'errors': []
        }
    print(Fore.BLUE + 'Running pylint...')
    for i, filei in enumerate(files):
        printed = False
        if i % 10 == 0:
            print(Fore.BLUE + '{} files remaining...'.format(len(files) - i))
        score, delta, warnings, errors = lint_file(
            filei,
            reg_score,
            reg_warning,
            reg_error,
            white_list=white_list,
        )
        if score is None:
            print(Fore.YELLOW + 'Skipping file: {}'.format(filei))
            continue

        if delta is not None:
            results['change'].append((filei, score, delta, abs(delta)))
            if delta < -0.5:
                results['big-change'].append((filei, score, delta))
                if not printed:
                    print(Fore.RED + '{} {} {}'.format(filei, score, delta))
                    printed = True

        if score < SCORE_THRESHOLD:
            results['below'].append((filei, score, delta))
            if not printed:
                print(Fore.RED + '{} {} {}'.format(filei, score, delta))
                printed = True

        if not printed:
            print(filei, score, delta)

        if warnings:
            results['warnings'].append((filei, score, delta))
            print_list(warnings, txt='warnings', color=Fore.YELLOW)
        if errors:
            results['errors'].append((filei, score, delta))
            print_list(errors, txt='errors', color=Fore.RED)

        results['all'].append((filei, score, delta))
    return results


def process_result(result):
    """Do a simple evaluation of the score."""
    try:
        score = float(result[0].split('/')[0])
    except IndexError:
        return None, None
    try:
        delta = score - float(result[1].split('/')[0])
    except IndexError:
        delta = None
    return score, delta


def print_list(list_to_print, txt='warnings', color=Fore.YELLOW):
    """Print out the contents in a list with some colour and text."""
    if list_to_print:
        print(color + '    Found {}:'.format(txt))
        for i in list_to_print:
            print(color + '    {}'.format(i.strip()))


def print_below(results):
    """Print the below results if any."""
    if results['below']:
        fmt = '    {}: {} (delta: {})'
        print()
        print(Fore.RED + 'The following files *MUST BE FIXED*:')
        for i in results['below']:
            print(fmt.format(*i))
        return 1
    return 0


def print_big_change(results):
    """Print big results if any."""
    if results['big-change']:
        fmt = '    {}: {} (delta: {})'
        print()
        print(Fore.RED +
              'The following files had a *TOO LARGE* negative change:')
        for i in results['big-change']:
            print(fmt.format(*i))
        return 1
    return 0


def print_all_warnings(results):
    """Print files with warnings if any."""
    if results['warnings']:
        fmt = '    {}'
        print()
        print(Fore.YELLOW +
              'The following files gave warnings:')
        for i in results['warnings']:
            print(fmt.format(i[0]))
        return 1
    return 0


def print_all_errors(results):
    """Print files with errors if any."""
    if results['errors']:
        fmt = '    {}'
        print()
        print(Fore.RED +
              'The following files have errors:')
        for i in results['errors']:
            print(fmt.format(i[0]))
        return 1
    return 0


def evaluate_results(results):
    """Do a simple evaluation of the results."""
    fmt = '    {}: {} (delta: {})'
    # First print out 10 (if possible) files with lowerest scores:
    maxlen = min(10, len(results['all']))
    txt = []
    for i in sorted(results['all'], key=itemgetter(1))[:maxlen]:
        if i[1] < 10.0:
            txt.append(fmt.format(*i))
    print()
    if txt:
        print(Fore.BLUE + 'The lowest scores:')
        for i in txt:
            print(i)
    else:
        print(Fore.GREEN + 'No scores below 10! Good job :-)')

    # Print out the 10 biggest changes:
    maxlen = min(10, len(results['change']))
    if maxlen > 0:
        txt = []
        for i in sorted(results['change'], key=itemgetter(3),
                        reverse=True)[:maxlen]:
            if abs(i[3]) > 0.0:
                txt.append(fmt.format(*i))
        print()
        if txt:
            print(Fore.BLUE + 'The biggest changes:')
            for i in txt:
                print(i)
        else:
            print()
            print(Fore.BLUE + 'No changes found.')

    fail1 = print_below(results)
    fail2 = print_big_change(results)
    print_all_warnings(results)
    fail3 = print_all_errors(results)
    return fail1 + fail2 + fail3


def create_argument_parser():
    """Create a argument parser for the linting tool."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i',
        '--input',
        nargs='+',
        help='Input directories to search for .py files.',
        required=True,
    )
    parser.add_argument(
        '-s',
        '--skip',
        nargs='+',
        help='Subdirectories to skip.',
        required=False,
    )
    parser.add_argument(
        '-w',
        '--white-list',
        nargs='+',
        help=(
            'Python modules to white-list, please see the pylint '
            'documentation on the option: '
            '"extension-pkg-whitelist"'
        ),
        required=False,
    )
    return parser


def main(args):
    """Run the analysis."""
    return_code = 0
    for startdir in args.input:
        print(Fore.GREEN + '\nRunning in: "{}"'.format(startdir))
        files = look_for_source_files(startdir, skip=args.skip)
        if files:
            lint_result = run_linting(files, white_list=args.white_list)
            result = evaluate_results(lint_result)
            return_code += result
        else:
            print(Fore.RED + 'No files found in {}!'.format(startdir))
            return_code += 1
    sys.exit(return_code)


if __name__ == '__main__':
    colorama.init(autoreset=True)
    PARSER = create_argument_parser()
    main(PARSER.parse_args())
