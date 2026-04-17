# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A test using the LAMMPS engine."""
from operator import itemgetter
import os
import colorama
import numpy as np
from pyretis.inout.common import make_dirs
from pyretis.engines.lammps import LAMMPSEngine
from pyretis.testing.helpers import clean_dir
from pyretis.testing.systemhelp import create_system_ext
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


HERE = os.path.abspath(os.path.dirname(__file__))
TRAJ = os.path.join(HERE, 'lammps_input', 'traj.lammpstrj')


def _convert_snapshot(snapshot):
    """Convert a LAMMPS text snapshot to numbers."""
    snapshot_new = {
        'timestep': int(snapshot['timestep'][0]),
        'number': int(snapshot['number'][0]),
    }
    box = []
    for line in snapshot['box']:
        box.append([float(i) for i in line.split()])
    snapshot_new['box'] = np.array(box)
    # Convert and sort atoms:
    atoms = []
    for line in snapshot['atoms']:
        atoms.append([float(i) for i in line.split()])
    snapshot_new['atoms'] = np.array(sorted(atoms, key=itemgetter(0)))
    return snapshot_new


def read_lammpstrj(filename):
    """Read frames in a LAMMPS trajectory."""
    snap = {}
    key = None
    data = []
    with open(filename) as infile:
        for lines in infile:
            strip_line = lines.strip()
            if strip_line.startswith('ITEM'):
                # This is the start of a new data block in a snapshot.
                # Add this data block to the current snapshot
                if key is not None:
                    snap[key] = data
                key = strip_line.split()[1].lower()
                data = []
                if key == 'timestep':
                    # This is the start of a new frame. Return the current
                    # snapshot, if any, and empty the data blocks:
                    if snap:
                        yield _convert_snapshot(snap)
                    snap = {}
                continue
            else:
                data.append(strip_line)
    if snap:
        if key is not None:
            snap[key] = data
        yield _convert_snapshot(snap)


def dump_phasepoint():
    """Use the LAMMPS engine to run a MD simulation forward in time."""
    logger.info('\nTesting that we can dump phase points.\n')
    engine = LAMMPSEngine('lmp_serial', 'lammps_input', 2,
                          extra_files=['dw-wca.in'])
    # Create a dummy system:
    system = create_system_ext(pos=(TRAJ, 0))
    exe_dir = os.path.join(HERE, 'dump')
    # Set up some directories:
    make_dirs(exe_dir)
    clean_dir(exe_dir)
    engine.exe_dir = exe_dir
    dumped_files = []
    for i in (0, 2, 4, 6, 8, 10):
        newpos = (TRAJ, i)
        logger.info(f'\t-> Dumping: {newpos}')
        system.particles.set_pos(newpos)
        engine.dump_phasepoint(system, f'dump-{i:02d}')
        pos = system.particles.get_pos()
        assert pos[1] == 0
        dumped_files.append(pos[0])
    logger.info('\nRunning some comparisons:\n')
    for i, frame in enumerate(read_lammpstrj(TRAJ)):
        traj = dumped_files[i]
        logger.info(f'\t-> Comparing for dumped frame: {traj}')
        frame_dump = [i for i in read_lammpstrj(traj)]
        assert len(frame_dump) == 1
        assert frame['number'] == frame_dump[0]['number']
        for key in ('box', 'atoms'):
            assert np.allclose(frame[key], frame_dump[0][key])
        assert frame_dump[0]['timestep'] == 0
        assert frame['timestep'] == engine.subcycles * i
        logger.info('\t\t->Frame was ok!')
    logger.info('\nConclusion: We can dump phase points.')


def main():
    """Run the comparisons."""
    dump_phasepoint()


if __name__ == '__main__':
    colorama.init(autoreset=True)
    main()
