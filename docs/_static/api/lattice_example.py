# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Example of creating a lattice."""
from pyretis.tools import generate_lattice
from pyretis.inout.formats.xyz import write_xyz_file
xyz, size = generate_lattice('diamond', [3, 3, 3], lcon=3.57)
write_xyz_file('diamond.xyz', xyz, names=['C']*len(xyz), header='Diamond')

xyz, size = generate_lattice('hcp', [3, 3, 3], lcon=2.5)
name = ['A', 'B'] * (len(xyz) // 2)
write_xyz_file('hcp.xyz', xyz, names=name, header='HCP test')

xyz, size = generate_lattice('sq2', [3, 3], lcon=1.0)
write_xyz_file('sq2.xyz', xyz, header='sq2 test')
