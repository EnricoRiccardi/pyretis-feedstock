# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Methods for comparing simulation results.

This module defines methods that can be used for comparing results
from different simulations, such as output files, reports, and
path ensembles.
"""
import math
import os
import filecmp
import itertools
import numpy as np
from pyretis.testing.helpers import search_for_files
from pyretis.inout.formats.energy import EnergyPathFile
from pyretis.inout.formats.order import OrderPathFile
from pyretis.inout.formats.path import PathExtFile
from pyretis.inout.formats.cross import CrossFile


# Names of the expected output files in archive directories:
ARCHIVE_FILES = {'energy.txt', 'order.txt', 'traj.txt'}
# Names of other expected output files:
OUTPUT_FILES = {'energy.txt', 'order.txt', 'pathensemble.txt'}
# Define readers for loading data:
READERS = {
    'energy': EnergyPathFile,
    'order': OrderPathFile,
    'traj': PathExtFile,
}


def read_files(*files, read_comments=True):
    """Read files into memory.

    Here, we assume that we are given small files and that we
    can read these into memory.

    Parameters
    ----------
    files : tuple of str
        These are the paths to the files we are to read.
    read_comments : bool, optional
        If False, we skip lines starting with a "#".

    Returns
    -------
    all_data : list of list of str
        The data read from the different files.
    """
    all_data = []
    for filename in files:
        data = []
        with open(filename, 'r', encoding="utf8") as infile:
            for line in infile:
                if not read_comments and line.strip().startswith('#'):
                    continue
                data.append(line)
        all_data.append(data)
    return all_data


def compare_text_line_by_line(file1, file2, skip=None, skip_keys=None):
    """Compare two files, line by line.

    Parameters
    ----------
    file1 : str
        The path to the first file to compare.
    file2 : str
        The path to the second file to compare.
    skip : list of int, optional
        These are 0-indexed line numbers we are to skip.
    skip_keys : list of str, optional
        Lines whose first token matches any key in this list are filtered
        out from both files before comparison. Useful for ignoring settings
        like ``exe_path`` that differ by run directory.

    Returns
    -------
    equal : bool
        True if the files are deemed to be equal.
    msg : str
        A descriptive message of the result of the comparison.
    """
    all_data = read_files(file1, file2, read_comments=True)
    assert len(all_data) == 2
    if skip_keys:
        def keep(line):
            """Return True if line should be kept."""
            token = line.split()[0] if line.split() else ''
            return token not in skip_keys
        data1 = [line for line in all_data[0] if keep(line)]
        data2 = [line for line in all_data[1] if keep(line)]
    else:
        data1, data2 = all_data[0], all_data[1]
    if len(data1) != len(data2):
        return False, 'The number of lines in the files differ'
    for i, (lini, linj) in enumerate(zip(data1, data2)):
        if skip and i in skip:
            continue
        if not lini.rstrip('\n') == linj.rstrip('\n'):
            return False, f'Line {i} differs: {lini.strip()} != {linj.strip()}'

    return True, 'Files are equal'


def compare_data_by_columns(file1, file2, file_type, skip=None):
    """Compare two output PyRETIS data files by columns.

    This method compares files where numbers are stored in columns
    and the columns have specific labels. Here, we also compare
    labels and comments.

    Parameters
    ----------
    file1 : str
        The path to the first file to compare.
    file2 : str
        The path to the second file to compare.
    file_type : str
        A string used to determine the file type (e.g., 'energy').
    skip : list of str, optional
        A list of items from the loaded data we are to skip.
        This can, for instance, be certain energy terms that are
        not absolute and can't easily be compared.

    Returns
    -------
    equal : bool
        True if the files are deemed to be equal.
    msg : str
        A descriptive message of the result of the comparison.
    """
    reader = READERS[file_type]
    data1 = reader(file1, 'r').load()
    data2 = reader(file2, 'r').load()
    # Compare the files by compare the block found in the file:
    for block1, block2 in zip(data1, data2):
        # Compare block comments, tolerating last-digit float differences.
        # Comments may embed order-parameter floats that differ by 1 ULP
        # between Python/numpy versions.
        if block1['comment'] != block2['comment']:
            if len(block1['comment']) != len(block2['comment']):
                return False, 'Block comment differs'
            for c1, c2 in zip(block1['comment'], block2['comment']):
                if c1 == c2:
                    continue
                # Try token-by-token float comparison
                t1, t2 = c1.split(), c2.split()
                if len(t1) != len(t2):
                    return False, 'Block comment differs'
                for tok1, tok2 in zip(t1, t2):
                    if tok1 == tok2:
                        continue
                    # Strip surrounding punctuation before float parse
                    stripped1 = tok1.strip("(),';")
                    stripped2 = tok2.strip("(),';")
                    try:
                        if not math.isclose(float(stripped1), float(stripped2),
                                            rel_tol=1e-12):
                            return False, 'Block comment differs'
                    except ValueError:
                        return False, 'Block comment differs'
        # Compare terms found in the blocks:
        if sorted(block1['data'].keys()) != sorted(block2['data'].keys()):
            return False, 'Different items in block data'
        # Compare numerical data:
        for key, val in block1['data'].items():
            if skip and key in skip:
                continue
            if not np.allclose(val, block2['data'][key]):
                return False, 'Block terms differ'
    return True, 'Files are equal'


def compare_numerical_data(file1, file2, rel_tol=1e-5):
    """Compare two files containing numerical data.

    Here, we compare files that contain numerical data. We don't
    care about comments here, we just compare the actual numerical data.

    Parameters
    ----------
    file1 : str
        The path to the first file to compare.
    file2 : str
        The path to the second file to compare.
    rel_tol : float, optional
        Relative tolerance for the comparison.

    Returns
    -------
    equal : bool
        True if the files are deemed to be equal.
    msg : str
        A descriptive message of the result of the comparison.
    """
    data1 = np.loadtxt(file1)
    data2 = np.loadtxt(file2)
    if not np.allclose(data1, data2, rtol=rel_tol, equal_nan=True):
        return False, 'Numerical data differ'
    return True, 'Files are equal'


def compare_numerical_mse(file1, file2, tol=1e-12):
    """Compare two numerical files using mean squared error.

    Parameters
    ----------
    file1 : str
        The path to the first file to compare.
    file2 : str
        The path to the second file to compare.
    tol : float, optional
        Tolerance for the mean squared error.

    Returns
    -------
    equal : bool
        True if the MSE is below the tolerance.
    msg : str
        A descriptive message with the MSE value.
    """
    data1 = np.loadtxt(file1)
    data2 = np.loadtxt(file2)
    if data1.shape != data2.shape:
        return False, f'Shapes differ: {data1.shape} != {data2.shape}'
    mse = np.mean((data1 - data2)**2)
    if mse > tol:
        return False, f'MSE {mse} > {tol}'
    return True, f'MSE {mse} is within tolerance'


def compare_restarted_text_files(file11, file12, file2):
    """Check if file2 is equal to file11 + file12 minus one overlapping line.

    We handle headers (lines starting with '#') by skipping them in the
    second file part.

    Parameters
    ----------
    file11 : str
        Path to the first part of the restarted simulation output.
    file12 : str
        Path to the second part of the restarted simulation output.
    file2 : str
        Path to the full (continuous) simulation output.

    Returns
    -------
    equal : bool
        True if the files match the pattern.
    msg : str
        A descriptive message of the result.
    """
    with open(file11, 'r', encoding='utf-8') as f11, \
         open(file12, 'r', encoding='utf-8') as f12, \
         open(file2, 'r', encoding='utf-8') as f2:

        f11_lines = f11.readlines()
        f12_lines = f12.readlines()
        f2_lines = f2.readlines()

        # Find first non-comment line in f12
        idx12 = 0
        while idx12 < len(f12_lines) and f12_lines[idx12].startswith('#'):
            idx12 += 1

        if idx12 >= len(f12_lines):
            return False, 'Part 2 of restarted file contains no data'

        # Check overlap
        line11_last = f11_lines[-1] if f11_lines else None
        line12_first_data = f12_lines[idx12]

        if line11_last != line12_first_data:
            msg = ('Overlapping lines differ between part 1 and '
                   'part 2')
            return False, msg

        # Combined = Part 1 + Part 2
        # (skipping its header and its first data line)
        combined = f11_lines + f12_lines[idx12 + 1:]

        if len(combined) != len(f2_lines):
            msg = f'Line count mismatch: {len(combined)} != {len(f2_lines)}'
            return False, msg

        for i, (l1, l2) in enumerate(zip(combined, f2_lines)):
            if l1 != l2:
                return False, f'Mismatch at line {i}'

    return True, 'Restarted files match the full simulation'


def compare_simulation_files(file1, file2, skip=None, mode='line'):
    """Top-level function to compare two simulation output files.

    Parameters
    ----------
    file1 : str
        The path to the first file to compare.
    file2 : str
        The path to the second file to compare.
    skip : list of str or list of int, optional
        A list of items that are to be skipped in the comparison.
    mode : str, optional
        A string used to determine how we do the comparison:
        'numerical' will select a comparison of numerical blocks;
        'line' will select a line-by-line text comparison;
        anything else will perform a literal file comparison.

    Returns
    -------
    equal : bool
        True if the files were found to be equal, False otherwise.
    msg : str
        A string with information about the comparison result.
    """
    if mode == 'numerical':
        return compare_numerical_data(file1, file2)
    if mode == 'line':
        return compare_text_line_by_line(file1, file2, skip=skip)
    equal = filecmp.cmp(file1, file2, shallow=False)
    msg = 'Files are equal' if equal else 'Files are not equal'
    return equal, msg


def compare_traj_archive(dir1, dir2):
    """Compare archived trajectories between two directories.

    These archives consist of trajectory information such as
    energies, order parameters and positions. Here, we verify that
    the output written by PyRETIS is identical in the two cases.

    Parameters
    ----------
    dir1 : str
        The path to the first directory to use in the comparison.
    dir2 : str
        The path to the second directory to use in the comparison.

    Returns
    -------
    errors : list of tuple
        This list contains the files which differed, if any.
    """
    errors = []
    files1 = sorted(search_for_files(dir1))
    files2 = sorted(search_for_files(dir2))
    # Are the number of files equal:
    if len(files1) != len(files2):
        errors.append((dir1, dir2))
        return errors
    # Compare the files that are written by PyRETIS:
    for file1, file2 in zip(files1, files2):
        basename1 = os.path.basename(file1)
        basename2 = os.path.basename(file2)
        if basename1 != basename2:
            errors.append((file1, file2))
            continue
        if basename1 in ARCHIVE_FILES:
            equal, _ = compare_simulation_files(file1, file2, mode='cmp')
            if not equal:
                errors.append((file1, file2))
    return errors


def compare_path_ensemble_data(file1, file2, rel_tol=1e-5, skip=None):
    """Compare two path ensemble files.

    We compare line-by-line, but skip comments and we check that
    numbers are close, as judged by the given relative tolarance.

    Parameters
    ----------
    file1 : str
        The path to the first file to consider in the comparison.
    file2 : str
        The path to the second file to consider in the comparison.
    rel_tol : float, optional
        A relative tolerance used to determine if numbers are equal.
    skip : list of int, optional
        These are columns we are to skip in the comparison.

    Returns
    -------
    equal : bool
        True if the files are equal, False otherwise.
    msg : str
        A message describing the result of the comparison.
    """
    all_data = read_files(file1, file2, read_comments=False)
    assert len(all_data) == 2
    if not len(all_data[0]) == len(all_data[1]):
        return False, 'The number of lines in the files differ'
    # Define the expected data types for the columns in the path
    # ensemble files:
    data_types = {
        0: int, 1: int, 2: int, 3: str, 4: str, 5: str, 6: int, 7: str, 8: str,
        9: float, 10: float, 11: int, 12: int, 13: float, 14: int, 15: int,
    }
    for i, (line1, line2) in enumerate(zip(*all_data)):
        stuff1 = line1.split()
        stuff2 = line2.split()
        for col, func in data_types.items():
            if skip and col in skip:
                continue
            if func == str:
                check = func(stuff1[col]) == func(stuff2[col])
            else:
                check = math.isclose(
                    func(stuff1[col]), func(stuff2[col]), rel_tol=rel_tol
                )
            if not check:
                return False, f'Files differ on line {i}, column {col}'

    return True, 'Files are equal'


def compare_reports_normalized(fil1, fil2):
    """Compare two reports, normalizing common version/time differences.

    This function ignores Docutils version meta-data, timestamps, and
    common spelling variations (grey/gray) in CSS to remain robust against
    environment differences.

    Parameters
    ----------
    fil1 : str
        The path to the first report to compare.
    fil2 : str
        The path to the second report to compare.

    Returns
    -------
    equal : bool
        True if reports are essentially equal.
    msg : str
        Description of mismatch if found.
    """
    def get_clean_lines(filepath):
        """Read lines and skip non-essential content."""
        clean_lines = []
        with open(filepath, 'r', encoding='utf-8') as infile:
            for line in infile:
                # Skip Docutils generator line
                if 'meta name="generator"' in line:
                    continue
                # Skip timestamps
                if 'generated by PyRETIS' in line or \
                   (line.strip().startswith('on ') and
                    (line.strip().endswith('.') or
                     line.strip().endswith('.</p>'))):
                    continue
                # Normalize gray/grey for CSS
                line = line.replace('color: gray', 'color: grey')
                # Skip specific ID lines in CSS that might shift
                if ':Id: $Id: html4css1.css' in line:
                    continue
                clean_lines.append(line)
        return clean_lines

    clean1 = get_clean_lines(fil1)
    clean2 = get_clean_lines(fil2)

    if len(clean1) != len(clean2):
        return False, 'Reports differ in number of contentful lines'

    for i, (l1, l2) in enumerate(zip(clean1, clean2)):
        if l1 != l2:
            return False, f'Mismatch at cleaned line {i}: {l1.strip()}'

    return True, 'Reports are essentially equal'


def compare_restarted_cross_files(file11, file12, file2):
    """Compare CrossFile data from a restarted simulation.

    Parameters
    ----------
    file11 : str
        Path to the first part of the crossing data.
    file12 : str
        Path to the second part of the crossing data.
    file2 : str
        Path to the full (continuous) crossing data.

    Returns
    -------
    equal : bool
        True if the crossing data matches.
    msg : str
        A descriptive message.
    """
    def load_flattened(fpath):
        blocks = list(CrossFile(fpath, 'r').load())
        flat_data = []
        for block in blocks:
            flat_data.extend(block['data'])
        return np.array(flat_data)

    data2 = load_flattened(file2)
    data11 = load_flattened(file11)
    data12 = load_flattened(file12)

    if data11.size > 0 and data12.size > 0:
        if np.array_equal(data11[-1], data12[0]):
            combined_data = np.vstack((data11, data12[1:]))
        else:
            combined_data = np.vstack((data11, data12))
    elif data11.size > 0:
        combined_data = data11
    else:
        combined_data = data12

    if combined_data.shape != data2.shape:
        msg = (f'Data shape mismatch: {combined_data.shape} != '
               f'{data2.shape}')
        return False, msg

    if not np.array_equal(combined_data, data2):
        return False, 'Crossing data mismatch'

    return True, 'Crossing data matches'
