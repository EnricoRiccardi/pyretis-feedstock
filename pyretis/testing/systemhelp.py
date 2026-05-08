# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Methods that might be useful for testing.

This module defines methods that are useful in connection with
systems.

"""
from pyretis.core.system import System
from pyretis.core.particles import ParticlesExt


def create_system_ext(pos=None, vel=False):
    """Create an external system with given positions and velocities."""
    system = System()
    system.particles = ParticlesExt(dim=3)
    if pos is None:
        system.particles.set_pos((None, None))
    else:
        system.particles.set_pos(pos)
    system.particles.set_vel(vel)
    return system
