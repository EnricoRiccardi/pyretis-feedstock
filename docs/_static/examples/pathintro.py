# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Example of interaction with a Path."""
import numpy as np
from pyretis.core.system import System
from pyretis.core.particles import Particles
from pyretis.core.path import Path
from pyretis.core.random_gen import RandomGenerator

path = Path(rgen=RandomGenerator(seed=0))  # Create empty path.

# Add some phase points to the path:
for i in range(10):
    phasepoint = System()
    phasepoint.order = [i]
    phasepoint.particles = Particles(dim=3)
    phasepoint.add_particle(np.zeros(3), vel=np.zeros(3))
    phasepoint.vpot = i
    phasepoint.ekin = i
    path.append(phasepoint)

# Loop over the phase points in the path:
print('Looping forward:')
for i, phasepoint in enumerate(path.phasepoints):
    print('Point {}. Order parameter = {}'.format(i, phasepoint.order))

# Loop over the phase points in the path:
print('Looping backward:')
for phasepoint in reversed(path.phasepoints):
    print('Order parameter = {}'.format(phasepoint.order))

# Get some randomly chosen shooting points:
print('Generating shooting points:')
for i in range(10):
    point, idx = path.get_shooting_point()
    print(
        'Shooting point {}: index = {}, order = {}'.format(i, idx, point.order)
    )
