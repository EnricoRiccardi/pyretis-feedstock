#! /usr/bin/env python3
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This is a mock engine for testing of an external GROMACS engine.

This script is intended to spoof the behaviour of the grompp program.
It does not actually do the things grompp does, it is just here so we
can have a script for testing the GROMACS engine.
"""
import sys
import os
from pyretis.inout.formats.gromacs import (
    write_gromacs_gro_file,
    read_gromacs_gro_file,
    write_gromos96_file,
    read_gromos96_file,
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


def fake_gmx_energy(args):
    """Fake the gmx energy command."""
    print('Running gmx energy...', file=sys.stdout)
    need_args = {'-f': None}
    simple_parser(args, need_args)
    check_that_files_exist(('-f',), need_args)
    with open('energy.xvg', 'w', encoding="utf8") as outfile:
        with open(need_args['-f'], 'r', encoding="utf8") as infile:
            for lines in infile:
                outfile.write(lines)


def fake_gmx_grompp(args):
    """Fake the gmx grompp command."""
    print('Running gmx grompp...', file=sys.stdout)
    need_args = {'-f': None, '-c': None, '-p': None, '-o': None}
    simple_parser(args, need_args)
    check_that_files_exist(('-f', '-c', '-p'), need_args)
    with open(need_args['-o'], 'w', encoding="utf8") as fileh:
        fileh.write('Mock GROMACS .tpr file. Input files were:\n')
        for key in ('-f', '-c', '-p'):
            fileh.write(f'{key} = {need_args[key]}\n')
        fileh.write('Input settings (.mdp) were:\n')
        with open(need_args['-f'], 'r', encoding="utf8") as infile:
            for lines in infile:
                fileh.write(lines)


def fake_gmx_trjconv(args):
    """Fake the gmx trjconv command."""
    print('Running gmx trjconv...', file=sys.stdout)
    need_args = {'-f': None, '-s': None, '-o': None, '-b': None,
                 '-dump': None}
    simple_parser(args, need_args)
    check_that_files_exist(('-f', '-s'), need_args)
    read_write_gromacs(need_args['-f'], need_args['-o'])


def read_write_gromacs(infile, outfile):
    """Basically just copy a file, but ensure we have velocities."""
    if infile.endswith('gro'):
        snapshot, xyz, vel, _ = read_gromacs_gro_file(infile)
        write_gromacs_gro_file(outfile, snapshot, xyz, vel)
    else:
        snapshot, xyz, vel, _ = read_gromos96_file(infile)
        write_gromos96_file(outfile, snapshot, xyz, vel)


def fake_gmx_converttpr(args):
    """Fake the gmx convert-tpr command."""
    print('Running gmx convert-tpr...', file=sys.stdout)
    need_args = {'-extend': None, '-s': None, '-o': None}
    simple_parser(args, need_args)
    check_that_files_exist(('-s',), need_args)
    with open(need_args['-o'], 'w', encoding="utf8") as outfile:
        with open(need_args['-s'], 'r', encoding="utf8") as infile:
            for lines in infile:
                outfile.write(lines)
        txt = need_args['-extend']
        outfile.write(f'Extend = {txt}\n')


if __name__ == '__main__':
    # pylint: disable=invalid-name
    known = {
        'energy': fake_gmx_energy,
        'grompp': fake_gmx_grompp,
        'convert-tpr': fake_gmx_converttpr,
        'trjconv': fake_gmx_trjconv,
    }
    args1 = None
    try:
        args1 = sys.argv[1]
    except IndexError:
        sys.exit(1)
    if args1 not in known:
        sys.exit(1)
    func = known[args1]
    func(sys.argv[1:])
