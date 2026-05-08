# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Example of creating a lattice and making a box."""
from pyretis.core import create_box
box1 = create_box()
print(box1)
box2 = create_box(cell=[10, 10, 10])
print(box2)
box3 = create_box(low=[0, -10, 10], high=[10, 10, 20],
                  periodic=[True, True, False])
print(box3)
