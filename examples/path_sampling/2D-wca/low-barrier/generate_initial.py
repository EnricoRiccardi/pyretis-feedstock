# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Here, we just generate the initial config using PyRETIS tools."""
from pyretis.tools.lattice import generate_lattice
from pyretis.inout.format.xyz import write_xyz_file


# pylint: disable=invalid-name
filename = 'initial.xyz'
density = 9.0 / 15.0
xyz, size = generate_lattice('sq', [3, 3], density=density)
names = ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'B', 'B']
print(f'Generating for density: {density}')
print(f'Box size: {size}')
print(f'Writing file: {filename}')
write_xyz_file(filename, xyz, names=names, header='Initial config')
