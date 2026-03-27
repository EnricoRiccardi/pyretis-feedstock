#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This is a mock engine for testing external GROMACS engines.

This script is inteneded to be used as a spoof of the mdrun command
for GROMACS.
"""
import sys
import os
from time import sleep
import numpy as np
from pyretis.inout.formats.gromacs import (
    read_gromos96_file,
    write_gromos96_file,
    write_trr_frame,
)


def simple_parser(args, need_args):
    """Do some simple parsing of arguments."""
    for i, arg in enumerate(args):
        if arg in need_args:
            need_args[arg] = args[i+1]
    # Check that all options are set:
    for key, val in need_args.items():
        if val is None:
            print(f'Missing {key}', file=sys.stderr, end='\n')
            sys.exit(1)


def check_that_files_exist(keys, args):
    """Check that files in input arguments actually exist."""
    for key in keys:
        if not os.path.isfile(args[key]):
            print(f'Missing file {args[key]}',
                  file=sys.stderr, end='\n')
            sys.exit(1)


def read_mock_tpr(filename):
    """Read info from the mock tpr file."""
    config = None
    steps = 0
    nstxout = 0
    gen = False
    with open(filename, 'r', encoding='utf-8') as inputfile:
        for lines in inputfile:
            if lines.startswith('-c'):
                config = lines.split()[-1]
            elif lines.startswith('nsteps ='):
                steps = int(lines.split()[-1])
            elif lines.startswith('gen_vel ='):
                gen = lines.strip().split()[-1] == 'yes'
            elif lines.startswith('nstxout ='):
                nstxout = int(lines.strip().split()[-1])
    frame, xyz, vel, _ = read_gromos96_file(config)
    config = [frame, xyz, vel]
    return config, steps, nstxout, gen


def mock_mdrun(args, crash=False, sleep_time=None):
    """Fake the gmx mdrun command."""
    print('Running gmx mdrun, well the mock version', file=sys.stdout)
    need_args = {'-s': None, '-deffnm': None, '-c': None}
    simple_parser(args, need_args)

    check_that_files_exist(('-s',), need_args)
    config, steps, nstxout, gen = read_mock_tpr(need_args['-s'])

    print('Writing energy data...', file=sys.stdout)
    write_mock_edr(f"{need_args['-deffnm']}.edr", steps,
                   nstxout, gen=gen)

    if crash:
        print('Crash error for testing.', file=sys.stderr)
        sys.exit(1)

    print('Writing trajectory data...', file=sys.stdout)
    print(config[2][0][0], file=sys.stdout)

    write_mock_trr(f"{need_args['-deffnm']}.trr", steps, nstxout,
                   config[1], config[2], pause=sleep_time)

    if config[2][0][0] > 0:
        xyz = config[1] + steps * np.ones_like(config[1])
    else:
        xyz = config[1] - steps * np.ones_like(config[1])

    print('Writing cpt file...', file=sys.stdout)
    write_mock_cpt(f"{need_args['-deffnm']}.cpt", steps, xyz)

    print('Writing final configuration...', file=sys.stdout)
    write_gromos96_file(need_args['-c'], config[0],
                        xyz, config[2])


def write_mock_trr(filename, steps, nstxout, xyz, vel, pause=None):
    # pylint: disable=too-many-arguments
    """Write a mock TRR file."""
    xyzc = np.copy(xyz)
    velc = np.copy(vel)
    append = False
    natoms = len(xyzc)
    for i in range(0, steps + 1):
        if i % nstxout == 0:
            data = {'box': np.eye(3), 'x': xyzc, 'v': velc,
                    'natoms': natoms, 'time': i*0.002, 'step': i,
                    'lambda': 0.0}
            write_trr_frame(filename, data, append=append)
            if not append:
                append = True
            if pause:
                sleep(pause)
        if vel[0][0] > 0:
            xyzc += np.ones_like(xyz)
        else:
            xyzc -= np.ones_like(xyz)
        velc -= np.ones_like(velc)


def write_mock_cpt(filename, steps, xyz):
    """Write a mock cpt file."""
    with open(filename, 'w', encoding='utf-8') as outfile:
        outfile.write(f'Last step: {steps}\n')
        for j in xyz:
            outfile.write('{:12.7f} {:12.7f} {:12.7f}\n'.format(*j))


def read_mock_cpt(filename):
    """Read coordinates."""
    xyz = []
    steps = -1
    with open(filename, 'r', encoding='utf-8') as infile:
        for i, lines in enumerate(infile):
            if i == 0:
                steps = int(lines.strip().split()[-1])
            else:
                xyz.append([float(i) for i in lines.strip().split()])
    return steps, np.array(xyz)


def write_mock_edr(filename, steps, nstxout, gen=False, start=0):
    """Write a mock edr file."""
    header = [
        '# This file was created by the MOCK gmx program.',
        '@    title "GROMACS Energies',
        '@    xaxis  label "Time (ps)"',
        '@    yaxis  label "(kJ/mol)"',
        '@TYPE xy',
        '@ view 0.15, 0.15, 0.75, 0.85',
        '@ legend on',
        '@ legend box on',
        '@ legend loctype view',
        '@ legend 0.78, 0.8',
        '@ legend length 2',
        '@ s0 legend "Potential"',
        '@ s1 legend "Total Energy"',
        '@ s2 legend "Kinetic En."',
    ]
    fmt = '{:10.6f}  {:10.6f}  {:10.6f}  {:10.6f}\n'
    mode = 'a' if start > 0 else 'w'
    with open(filename, mode, encoding='utf-8') as output:
        if mode == 'w':
            for line in header:
                output.write(f'{line}\n')
        for i in range(0, steps + 1):
            j = float(i)
            if not gen:
                if i % nstxout == 0:
                    output.write(fmt.format(j, -j, j * 2, j))
            else:
                output.write(fmt.format(j, -j, j * 2, j + 1234.5678))


if __name__ == '__main__':
    # pylint: disable=invalid-name
    print('This is the mock mdrun command.')
    if '-crash' in sys.argv:
        mock_mdrun(sys.argv[1:], crash=True)
    else:
        if '-sleep' in sys.argv:
            mock_mdrun(sys.argv[1:], crash=False, sleep_time=0.15)
        else:
            mock_mdrun(sys.argv[1:], crash=False)
