# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This script will help to copy files for restarting.

Here, we copy the last accepted path for each ensemble.
"""
import os
import shutil
import sys
import pickle
from pyretis.inout.settings import parse_settings_file
from pyretis.core.pathensemble import generate_ensemble_name


def main(source_dir, target_dir):
    """Copy the files.

    Parameters
    ----------
    source_dir : string
        Path to the directory containing the ensemble directories we
        are copying from.
    target_dir : string
        Path to the directory where we will be copying the accepted
        trajectories to.

    """
    settings = parse_settings_file(os.path.join(source_dir, 'retis.rst'))
    nint = len(settings['simulation']['interfaces'])
    for i in range(nint):
        ens = generate_ensemble_name(i)
        target = os.path.join(target_dir, ens)
        source = os.path.join(source_dir, ens)
        shutil.copytree(source, target)
        # Read restart file and update paths:
        restart = os.path.join(target, 'ensemble.restart')
        info = {}
        with open(restart, 'rb') as infile:
            info = pickle.load(infile)
        newpos = []
        for phasepoint in info['path_ensemble']['last_path']['phasepoints']:
            filename, idx = phasepoint['particles']['config']
            name = os.path.basename(filename)
            abs_path = os.path.abspath(os.path.join(target, 'accepted', name))
            newpos.append((abs_path, idx))
        for phasepoint, pos in zip(
                info['path_ensemble']['last_path']['phasepoints'], newpos):
            phasepoint['particles']['config'] = pos
        with open(restart, 'wb') as outfile:
            pickle.dump(info, outfile)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
