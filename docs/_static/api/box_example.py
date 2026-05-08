# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Example of creating a lattice and making a box."""
from pyretis.core import create_box
from pyretis.tools import generate_lattice
from pyretis.inout.formats.xyz import write_xyz_file

xyz, size = generate_lattice('diamond', [2, 2, 2], lcon=3.567)
write_xyz_file('diamond_small.xyz', xyz, names=['C']*len(xyz))

xyz, size = generate_lattice('diamond', [10, 10, 10], lcon=3.567)
write_xyz_file('diamond.xyz', xyz, names=['C']*len(xyz))

low = [i[0] for i in size]
high = [i[1] for i in size]
box = create_box(low=low, high=high)
print(box)
