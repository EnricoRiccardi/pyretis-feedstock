# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Example of creating a simple order parameter."""
from pyretis.core import System, Particles
from pyretis.orderparameter import Position

position = Position(0, dim='x')

# Define a simple system for testing:
system = System()
system.particles = Particles(dim=3)
system.add_particle([1.0, 2.0, 3.0], mass=1.0, name='Ar', ptype=0)
print('Order =', position.calculate(system))
