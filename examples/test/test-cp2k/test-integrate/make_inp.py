# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Update the number of steps for a CP2K input file."""
import sys
from pyretis.inout.settings import parse_settings_file


def main(settings_file, inp_file, output_file):
    """Update number of steps in the CP2K input file."""
    settings = parse_settings_file(settings_file)
    steps = settings['simulation']['steps']
    subcycles = settings['engine']['subcycles']
    gmx_steps = steps * subcycles
    with open(inp_file, 'r') as infile:
        with open(output_file, 'w') as output:
            for lines in infile:
                if lines.strip().startswith('STEPS'):
                    output.write(f'STEPS = {gmx_steps}\n')
                else:
                    output.write(lines)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
