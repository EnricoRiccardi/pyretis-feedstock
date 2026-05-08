#! /usr/bin/env python3
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This is a mock engine for external testing.

This script is inteneded to be used as an external engine and it
will simply just create the output files.
"""
import sys
import os
import numpy as np
from pyretis.inout.formats.xyz import write_xyz_trajectory


def write_energy(outfile, steps=10):
    """Store CP2K style energies."""
    header = ' '.join(['# Step Nr.', 'Time[fs]', 'Kin.[a.u.]',
                       'Temp[K]', 'Pot.[a.u.]', 'Cons Qty[a.u.]',
                       'UsedTime[s]'])
    fmt = '{:10d} {:8.6f} {:8.6f} {:8.6f} {:8.6f} {:8.6f} {:8.6f}'
    with open(outfile, 'w', encoding='utf-8') as output:
        output.write(f'{header}\n')
        for i in range(steps):
            time = i * 0.5
            kin = i * 0.1
            temp = i * 100
            pot = -0.1 * i
            txt = fmt.format(i, time, kin, temp, pot, i, i)
            output.write(f'{txt}\n')


def write_wfn_file(outfile):
    """Write some gibberish data to represent a wfn file."""
    if os.path.isfile(outfile):
        # Just keep one backup file. This is useful for the testing.
        backup = f'{outfile}.bak-1'
        os.rename(outfile, backup)
    with open(outfile, 'w', encoding='utf-8') as output:
        output.write('Some gibberish')


def write_xyz_file(outfile, factor, steps=10):
    """Write some mock xyz-data."""
    names = ['H', 'H']
    box = [1., 2., 3.]
    fac = np.array([factor, np.ones(3) + factor])
    for i in range(steps):
        pos = np.ones((2, 3)) * i * fac
        vel = np.zeros((2, 3))
        write_xyz_trajectory(outfile, pos, vel, names, box, step=i,
                             append=True)


def read_project_name(filename):
    """Get the CP2K project name."""
    with open(filename, 'r', encoding='utf-8') as infile:
        for lines in infile:
            if lines.find('PROJECT') != -1:
                projname = lines.split()[-1]
                return projname
    return None


def write_cp2k_restart(filename, factorx, factorv, steps=10):
    """Write a CP2K-like restart file."""
    xyz = (steps - 1) * factorx
    vel = (steps - 1) * factorv
    if os.path.isfile(filename):
        # Just keep one backup file. This is useful for the testing.
        backup = f'{filename}.bak-1'
        os.rename(filename, backup)
    with open(filename, 'w', encoding='utf-8') as outfile:
        outfile.write('&FORCE_EVAL\n')
        outfile.write('&SUBSYS\n')
        outfile.write('&COORD\n')
        for i in range(2):
            xyzi = xyz + np.ones(3) * i
            outfile.write(
                '{:s} {:8.6f} {:8.6f} {:8.6f}\n'.format('H', *xyzi)
            )
        outfile.write('&END COORD\n')
        outfile.write('&CELL\n')
        outfile.write(
            '{:s} {:8.6f} {:8.6f} {:8.6f}\n'.format('A', 1., 0., 0.)
        )
        outfile.write(
            '{:s} {:8.6f} {:8.6f} {:8.6f}\n'.format('B', 0., 2., 0.)
        )
        outfile.write(
            '{:s} {:8.6f} {:8.6f} {:8.6f}\n'.format('C', 0., 0., 3.)
        )
        outfile.write('&END CELL\n')
        outfile.write('&VELOCITY\n')
        for i in range(2):
            veli = vel + np.ones(3) * i
            outfile.write(
                '{:8.6f} {:8.6f} {:8.6f}\n'.format(*veli)
            )
        outfile.write('&END VELOCITY\n')
        outfile.write('&END SUBSYS\n')
        outfile.write('&END FORCE_EVAL\n')


if __name__ == '__main__':
    # pylint: disable=invalid-name
    print('This is the mock CP2K command')
    args = sys.argv
    name = read_project_name(args[-1])
    if not name:
        sys.exit(1)

    efile = f'{name}-1.ener'
    write_energy(efile, steps=10)

    wfile = f'{name}-RESTART.wfn'
    write_wfn_file(wfile)

    xfile = f'{name}-pos-1.xyz'
    xfac = np.array([0.1, 0.2, 0.3])
    write_xyz_file(xfile, xfac, steps=10)

    vfile = f'{name}-vel-1.xyz'
    vfac = np.array([1.1, 1.2, 1.3])
    write_xyz_file(vfile, vfac, steps=10)

    rfile = f'{name}-1.restart'
    write_cp2k_restart(rfile, xfac, vfac, steps=10)
