# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Result comparison for GROMACS RETIS simulation.

Here we compare a RETIS simulation of 250 steps to known results
stored in a gzipped tar archive.
"""
import os
import sys
import tarfile
import colorama
from pyretis.inout.settings import parse_settings_file
from pyretis.core.pathensemble import generate_ensemble_name
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


RESULTS = 'results'
RESULTS_TGZ = {
    '5.1.4': 'results-5.1.4.tgz',
}


def read_tarfile_content(tar_path, inner_path):
    """Extract and read a file's content from a .tgz archive.

    Parameters
    ----------
    tar_path : string
        Path to the .tgz file.
    inner_path : string
        The internal path of the file to extract.

    Returns
    -------
    content : bytes
        The raw content of the extracted file, or None if not found.
    """
    with tarfile.open(tar_path, 'r:gz') as tar:
        for member in tar.getmembers():
            if member.isreg() and member.name == inner_path:
                with tar.extractfile(member) as tfile:
                    return tfile.read()
    return None


def find_sorted_traj_files(rootdir):
    """Find and sort all traj.txt files in a directory tree.

    Parameters
    ----------
    rootdir : string
        The root directory to search in.

    Returns
    -------
    sorted_files : list of tuples
        A list of (filepath, cycle_number) sorted by cycle.
    """
    found_files = []
    for root, _, files in os.walk(rootdir):
        if 'traj.txt' in files:
            filepath = os.path.join(root, 'traj.txt')
            with open(filepath, 'r', encoding='utf-8') as infile:
                line = infile.readline()
                # Expected format: "# Cycle:1, Name: ..."
                try:
                    cycle = int(line.split(':')[1].split(',')[0])
                    found_files.append((filepath, cycle))
                except (IndexError, ValueError):
                    continue
    return sorted(found_files, key=lambda x: x[1])


def recreate_traj_txt(root, ensemble):
    """Recreate a consolidated traj.txt from individual step files.

    Parameters
    ----------
    root : string
        Root directory of the run.
    ensemble : int
        The ensemble index.
    """
    ens_name = generate_ensemble_name(ensemble)
    traj_dir = os.path.join(root, ens_name, 'traj')
    sorted_trajs = find_sorted_traj_files(traj_dir)
    target_path = os.path.join(root, ens_name, 'traj.txt')
    with open(target_path, 'w', encoding='utf-8') as outfile:
        for filepath, _ in sorted_trajs:
            with open(filepath, 'r', encoding='utf-8') as infile:
                outfile.write(infile.read())


def compare_ensemble_results(settings, root, results_tgz):
    """Compare all expected output files for a simulation.

    Parameters
    ----------
    settings : dict
        Simulation settings.
    root : string
        The directory containing current results.
    results_tgz : string
        The archive containing reference results.

    Returns
    -------
    status : bool
        True if all files match, False otherwise.
    """
    inter = settings['simulation']['interfaces']
    for i in range(len(inter)):
        ens_name = generate_ensemble_name(i)
        for fname in ('pathensemble.txt', 'order.txt', 'traj.txt'):
            if fname == 'traj.txt':
                recreate_traj_txt(root, i)

            # Compare against the reference in the tarball
            ref_path = os.path.join(RESULTS, ens_name, fname)
            tar_path = os.path.join(root, results_tgz)
            reference_data = read_tarfile_content(tar_path, ref_path)

            current_path = os.path.join(root, ens_name, fname)
            logger.info(f'Comparing: {current_path}')

            if reference_data is None:
                logger.error(f'\t-> *Reference missing: {ref_path}*')
                return False

            with open(current_path, 'rb') as current_file:
                if current_file.read() != reference_data:
                    logger.error('\t-> *Files differ!*')
                    return False
            logger.info('\t-> Files are equal!')

    logger.info('All files are equal!')
    return True


def main(argv):
    """Parse GROMACS version and run comparisons.

    Parameters
    ----------
    argv : list
        Command line arguments.
    """
    try:
        gmx_version = argv[1].split()[-1]
    except (IndexError, AttributeError):
        gmx_version = '5.1.4'

    logger.info(f'GROMACS version family: {gmx_version}')
    correct_tgz = RESULTS_TGZ.get(gmx_version)

    if not correct_tgz:
        logger.info(f'GROMACS version {gmx_version} is not supported '
                    'for result comparison.')
        logger.info('Supported versions:')
        for ver in RESULTS_TGZ:
            logger.info(f'- {ver}')
        sys.exit(0)

    logger.info(f'Using results from: {correct_tgz}')

    for dirname in ('gromacs1', 'gromacs2'):
        if not os.path.isdir(dirname):
            continue
        sets = parse_settings_file(os.path.join(dirname, 'retis.rst'))
        header = f'Comparing files for: {dirname}'
        logger.info(f'\n{header}')
        logger.info('=' * len(header))

        if not compare_ensemble_results(sets, dirname, correct_tgz):
            sys.exit(1)


if __name__ == '__main__':
    colorama.init(autoreset=True)
    main(sys.argv)
