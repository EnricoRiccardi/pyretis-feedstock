#! /usr/bin/env python3
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This is a mock engine for external testing.

This script is inteneded to be used as an external engine and will
create some LAMMPS-like output.
"""
import sys
import os
import tempfile
from pyretis.engines.lammps import read_lammps_input


def get_arguments(lammps_args):
    """Read arguments given to this program.

    Note, we here assuming that the arguments we are served
    are on the form: -in <inputfile> -l <logfile> -screen <screenfile>

    """
    looking_for = {'-in', '-l', '-screen'}
    found = {}
    for i, arg in enumerate(lammps_args):
        if arg in looking_for:
            var = lammps_args[i+1]
            found[arg] = var
    errors = False
    for key in looking_for:
        if key not in found:
            errors = True
            print(f'ERROR: Missing input argument "{key}"')
    if errors:
        sys.exit(1)
    return found


def write_fake_log(logfile, steps, freq):
    """Create a fake log file for LAMMPS."""
    with open(logfile, 'w', encoding='utf-8') as output:
        output.write('This is not LAMMPS.\n')
        output.write('But we are pretending to be.\n\n')
        output.write('Please have some fake thermo data:\n')
        output.write('Step Temp Press PotEng KinEng TotEng\n')
        for i in range(steps + 1):
            if i % freq == 0:
                j = 5 * i * freq
                output.write(
                    f'{i} {j + 1} {j + 2} {j + 3} {j + 4} {j + 5}\n'
                )
        output.write('Loop time of 0 on 1 procs for 12 steps with 0 atoms\n')
        output.write('\n')
        output.write('Total wall time: 0:00:00')


def write_fake_screen(screenfile, settings):
    """Create a fake output file for LAMMPS.

    Here, we just repeat the input settings.

    """
    with open(screenfile, 'w', encoding='utf-8') as output:
        for key, val in settings:
            output.write(f'{key} {val}\n')


def make_dump_file(filename, steps, freq):
    """Create a fake dump file."""
    with open(filename, 'w', encoding='utf-8') as output:
        output.write('# Fake dump file.')
        for i in range(steps + 1):
            if i % freq == 0:
                output.write('ITEM: TIMESTEP\n')
                output.write(f'{i}\n')


def make_restart_file(filename):
    """Create a fake restart file."""
    with open(filename, 'w', encoding='utf-8') as output:
        output.write('Fake restart file.')


def make_order_file(filename, steps, freq):
    """Create a fake order parameter file."""
    with open(filename, 'w', encoding='utf-8') as output:
        output.write('# Fake order parameter file.\n')
        for i in range(steps + 1):
            if i % freq == 0:
                output.write(
                    f'{i} {i + 1} {-(i + 1)}\n'
                )


def make_fake_output_files(settings, steps):
    """Create output files as requested by the settings."""
    for key, val in settings:
        if key == 'dump':
            # Assume that the setting is: dump ID group-ID style N filename
            filename = val.split()[4].strip()
            freq = int(val.split()[3])
            print(f'Dumping to: "{filename}"')
            make_dump_file(filename, steps, freq)
        elif key == 'write_restart':
            filename = val.split()[0].strip()
            print(f'Writing restart to: "{filename}"')
            make_restart_file(filename)
        elif key == 'fix':
            # Check if we are to write order parameters.
            if val.startswith('order_output'):
                val_split = [i.strip() for i in val.split()]
                idx = val_split.index('append')
                filename = val_split[idx+1]
                idx = val_split.index('print')
                freq = int(val_split[idx+1])
                print(f'Writing order parameters to: "{filename}"')
                make_order_file(filename, steps, freq)


def main(lammps_args):
    """Run the mock LAMMPS engine."""
    arguments = get_arguments(lammps_args)
    settings = read_lammps_input(arguments['-in'])
    # First just figure out the number of steps:
    steps = 0
    freq = 1
    for key, val in settings:
        if key == 'run' and len(val.split()) == 1:
            steps = int(val)
        if key == 'thermo' and len(val.split()) == 1:
            freq = int(val)
    # Create the output files
    write_fake_log(arguments['-l'], steps, freq)
    write_fake_screen(arguments['-screen'], settings)
    make_fake_output_files(settings, steps)


if __name__ == '__main__':
    print('This is the mock LAMMPS command')
    try:
        main(sys.argv)
    except BaseException as err:
        import traceback
        tmpdir = tempfile.gettempdir()
        errlog = os.path.join(tmpdir, 'mocklammps_error.log')
        with open(errlog, 'w') as f:
            f.write(f"ERROR in mocklammps.py: {err}\n")
            traceback.print_exc(file=f)
        sys.exit(1)
