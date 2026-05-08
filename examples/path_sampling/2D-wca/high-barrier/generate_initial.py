# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Here, we just generate the initial config using PyRETIS tools."""
from pyretis.tools.lattice import generate_lattice
from pyretis.inout.formats.xyz import write_xyz_file

# pylint: disable=invalid-name
filename = 'initial.xyz'
density = 0.7
xyz, size = generate_lattice('sq', [5, 5], density=density)
names = ['B', 'B'] + ['A'] * 23
print(f'Generating for density: {density}')
print(f'Box size: {size}')
print(f'Writing file: {filename}')
write_xyz_file(filename, xyz, names=names, header='Initial config')
