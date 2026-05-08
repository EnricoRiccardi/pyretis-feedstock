# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Example of interaction with a Composite order parameter."""
import numpy as np
from pyretis.core import System, Particles
from pyretis.orderparameter import (
    OrderParameter,
    Position,
    CompositeOrderParameter,
)

# Create a new, empty, order parameter:
order_parameter = CompositeOrderParameter()

# Add a position order parameter:
order_parameter.add_orderparameter(Position(0, dim='x'))

# Add a custom order parameter:
position_y = OrderParameter(description='Position along y-axis')


def collective_y(system):
    """Position along y-axis."""
    return [system.particles.pos[0][1]]


position_y.calculate = collective_y
order_parameter.add_orderparameter(position_y)

# Add another order parameter:
position_cos = OrderParameter(description='Cosine of z-coordinate')


def collective_cos(system):
    """Additional collective variable."""
    return [np.cos(np.pi * system.particles.pos[0][2])]


position_cos.calculate = collective_cos
order_parameter.add_orderparameter(position_cos)

# Create a system for testing the new order parameter:
i_sys = System()
i_sys.particles = Particles(dim=3)
i_sys.add_particle([1.0, 2.0, 3.0], mass=1.0, name='Ar', ptype=0)
print('Order =', order_parameter.calculate(i_sys))
