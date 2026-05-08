# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This script will help to copy trjectory files needed for restarting.

Here, we copy the last accepted path(s) from a previous simulation.
"""
import os
import sys
import shutil
import colorama
from pyretis.inout.archive import PathStorage
from pyretis.inout.settings import parse_settings_file
from pyretis.core.pathensemble import generate_ensemble_name
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def read_path_file(filename):
    """Get index for last accepted path."""
    last_idx = None
    with open(filename, 'r') as fileh:
        for i, lines in enumerate(fileh):
            if i == 0:
                continue
            split = lines.strip().split()
            idx = int(split[0])
            status = split[7]
            if status == 'ACC':
                last_idx = idx
    return last_idx


def get_files_from_directory(source, target):
    """Investigate and copy an accepted path for the given ensemble."""
    logger.info(f'Checking directory:\xa0{source}')
    path_file = os.path.join(source, 'pathensemble.txt')
    last_one = read_path_file(path_file)
    logger.info(f'Will use path no. {last_one}')

    # Copy the accepted last path from the traj archive.
    # This is the path where it is stored:
    traj_archive = os.path.join(
        source,
        'traj',  # Note, this is hardcoded in PyRETIS
        PathStorage.archive_acc,
        PathStorage.out_dir_fmt.format(last_one)
    )
    # We first copy the order, energy and traj files:
    for filei in ('order.txt', 'energy.txt', 'traj.txt'):
        src = os.path.join(traj_archive, filei)
        dest = os.path.join(target, filei)
        shutil.copy(src, dest)
    # Next, copy the actual trajectory files. Here we copy the files
    # corresponding to the last accepted trajectory found in the given
    # source. This trajectory is copied into a sub-directory called "accepted"
    # of the target directory.
    traj_archive_raw = os.path.join(traj_archive, PathStorage.archive_traj)
    for i in os.scandir(traj_archive_raw):
        if i.is_file:
            src = os.path.join(traj_archive_raw, i.name)
            target_dir = os.path.join(target, 'accepted')
            try:
                os.makedirs(target_dir)
            except FileExistsError:
                pass
            dest = os.path.join(target_dir, i.name)
            shutil.copy(src, dest)


def main(source, target):
    """Copy the files."""
    settings = parse_settings_file(os.path.join(source, 'retis.rst'))
    interf = settings['simulation']['interfaces']
    for i, _ in enumerate(interf):
        ens = generate_ensemble_name(i)
        target_ens = os.path.join(target, ens)
        source_ens = os.path.join(source, ens)
        try:
            os.makedirs(target_ens)
        except FileExistsError:
            pass
        for filei in ('pathensemble.txt', 'ensemble.restart'):
            src = os.path.join(source_ens, filei)
            dest = os.path.join(target_ens, filei)
            shutil.copy(src, dest)
        get_files_from_directory(source_ens, target_ens)


if __name__ == '__main__':
    colorama.init(autoreset=True)
    main(sys.argv[1], sys.argv[2])
