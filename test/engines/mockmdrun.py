#! /usr/bin/env python3
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This is a mock engine for testing external GROMACS engines.

This script is inteneded to be used as a spoof of the mdrun command
for GROMACS.
"""
import sys
import os
import numpy as np
from pyretis.inout.formats.gromacs import (
    read_gromacs_gro_file,
    write_gromacs_gro_file,
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
    gen = False
    with open(filename, 'r', encoding='utf-8') as inputfile:
        for lines in inputfile:
            if lines.startswith('-c'):
                config = lines.split()[-1]
            elif lines.startswith('nsteps ='):
                steps = int(lines.split()[-1])
            elif lines.startswith('gen_vel ='):
                gen = lines.strip().split()[-1] == 'yes'
    frame, xyz, vel, _ = read_gromacs_gro_file(config)
    vel = np.ones_like(xyz)
    config = [frame, xyz, vel]
    return config, steps, gen


def mock_mdrun(args):
    """Fake the gmx mdrun command."""
    print('Running the mock gmx mdrun program...', file=sys.stdout)
    need_args = {'-s': None, '-deffnm': None, '-c': None}
    simple_parser(args, need_args)
    check_that_files_exist(('-s',), need_args)
    config, steps, gen = read_mock_tpr(need_args['-s'])
    print('Writing log file...', file=sys.stdout)
    logfile = f"{need_args['-deffnm']}.log"
    with open(logfile, 'w', encoding='utf-8') as output:
        output.write('Mock GROMACS log file\n')
    print('Writing energy file...', file=sys.stdout)
    write_mock_edr(f"{need_args['-deffnm']}.edr", steps, gen=gen)
    print('Writing trr file...', file=sys.stdout)
    write_mock_trr(f"{need_args['-deffnm']}.trr", steps, config[1],
                   start=0)
    xyz = config[1] + steps * np.ones_like(config[1])
    print('Writing cpt file...', file=sys.stdout)
    write_mock_cpt(f"{need_args['-deffnm']}.cpt", steps, xyz)
    print('Writing final configuration file...', file=sys.stdout)
    write_gromacs_gro_file(need_args['-c'], config[0],
                           xyz, config[2])


def mock_mdrun_continue(args):
    """Fake the gmx mdrun command."""
    print('Running the mock gmx mdrun command with -append', file=sys.stdout)
    need_args = {'-s': None, '-deffnm': None, '-c': None,
                 '-cpi': None}
    simple_parser(args, need_args)
    check_that_files_exist(('-s', '-cpi'), need_args)
    config, steps, gen = read_mock_tpr(need_args['-s'])
    stepsp, xyz = read_mock_cpt(need_args['-cpi'])
    print('Writing log file...', file=sys.stdout)
    logfile = f"{need_args['-deffnm']}.log"
    with open(logfile, 'a', encoding='utf-8') as output:
        output.write('Mock GROMACS log file\n')
    print('Writing edr file...', file=sys.stdout)
    write_mock_edr(f"{need_args['-deffnm']}.edr", steps, gen=gen,
                   start=stepsp)
    print('Writing trr file...', file=sys.stdout)
    write_mock_trr(f"{need_args['-deffnm']}.trr", steps, xyz,
                   start=stepsp)
    xyzf = config[1] + (steps + stepsp) * np.ones_like(config[1])
    print('Writing cpt file...', file=sys.stdout)
    write_mock_cpt(f"{need_args['-deffnm']}.cpt",
                   steps + stepsp, xyzf)
    print('Writing final configuration file...', file=sys.stdout)
    write_gromacs_gro_file(need_args['-c'], config[0],
                           xyzf, config[2])
    print('Updating log...', file=sys.stdout)
    backuplog = f"#{need_args['-deffnm']}.log"
    with open(backuplog, 'a', encoding='utf-8') as output:
        output.write('Mock GROMACS log file. Backup!\n')


def write_mock_trr(filename, steps, xyz, start=0):
    """Write a mock TRR file."""
    xyzc = np.copy(xyz)
    mode = 'a' if start > 0 else 'w'
    with open(filename, mode, encoding='utf-8') as outfile:
        for i in range(start, start + steps + 1):
            write = True
            if steps != 0 and i % steps != 0:
                write = False
            if mode == 'a' and i == start:
                write = False
            if write:
                outfile.write(f'Step: {i}\n')
                for j in xyzc:
                    outfile.write('{:12.7f} {:12.7f} {:12.7f}\n'.format(*j))
            xyzc += np.ones_like(xyz)


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


def write_mock_edr(filename, steps, gen=False, start=0):
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
        for i in range(start, start + steps + 1):
            j = float(i)
            if not gen:
                write = True
                if steps != 0 and i % steps != 0:
                    write = False
                if mode == 'a' and i == start:
                    write = False
                if write:
                    output.write(fmt.format(j, -j, j * 2, j))
            else:
                output.write(fmt.format(j, -j, j * 2, j + 1234.5678))


if __name__ == '__main__':
    # pylint: disable=invalid-name
    print('This is the mock mdrun command')
    cont = '-append' in sys.argv
    if cont:
        mock_mdrun_continue(sys.argv[1:])
    else:
        mock_mdrun(sys.argv[1:])
