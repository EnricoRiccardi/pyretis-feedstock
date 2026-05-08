# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the C implementation of the Velocity Verlet integrator.

This test is comparing:

1) The Python (numpy) implementation

2) The C implementation.
"""
# pylint: disable=invalid-name
import unittest
import numpy as np
from pyretis.core import System, create_box, Particles
from pyretis.simulation import Simulation
from pyretis.core.units import create_conversion_factors
from pyretis.forcefield import ForceField
from pyretis.forcefield.potentials import PairLennardJonesCutnp
from pyretis.tools import generate_lattice
from pyretis.engines import VelocityVerlet
from vvintegratorc import VelocityVerletC


def create_positions():
    """Create particles for the test."""
    create_conversion_factors('lj')
    lattice, size = generate_lattice('fcc', [3, 3, 3], density=0.9)
    size = np.array(size)
    box = create_box(low=size[:, 0], high=size[:, 1],
                     periodic=[True, True, True])
    npart = len(lattice)
    lattice += np.random.randn(npart, 3) * 0.05
    particles = Particles(dim=3)
    for pos in lattice:
        particles.add_particle(pos, np.zeros_like(pos), np.zeros_like(pos),
                               mass=1.0, name='Ar', ptype=0)
    msg = 'Created lattice with {} atoms.'
    print(msg.format(particles.npart))
    return particles, box


def run_test(steps, integrator, system=None):
    """Execute a test MD simulation."""
    if system is None:
        # create system
        particles, box = create_positions()
        system = System(temperature=1.0, units='lj', box=box)
        system.particles = particles
        initial = system.copy()
    else:
        initial = None
    parameters = {0: {'sigma': 1.0, 'epsilon': 1.0, 'rcut': 2.5}}
    potentialnp = PairLennardJonesCutnp(dim=3, shift=True)
    forcefieldnp = ForceField('Lennard-Jones force field',
                              potential=[potentialnp], params=[parameters])
    system.forcefield = forcefieldnp
    controls = {'steps': steps}
    simulation = Simulation(settings=None, controls=controls)
    task_integrate = {'func': integrator.integration_step,
                      'args': [system]}
    simulation.add_task(task_integrate)
    traj = []
    for _ in simulation.run():
        traj.append(
            {
                'pos': system.particles.get_pos(),
                'vel': system.particles.get_vel(),
                'force': system.particles.get_force()
            }
        )
    return traj, initial, system


class VVIntegratorTest(unittest.TestCase):
    """Run the tests for the C integrator."""

    def test_integrator(self):
        """Test by integrating the equations of motion."""
        integrator = VelocityVerlet(0.0025)
        traj, initial, system = run_test(20, integrator)
        # reset to initial state
        system = initial
        # repeat with external integrator:
        integratorc = VelocityVerletC(0.0025)
        traj2, _, _ = run_test(20, integratorc, system)
        for trj1, trj2 in zip(traj, traj2):
            posok = np.allclose(trj1['pos'], trj2['pos'])
            self.assertTrue(posok)
            velok = np.allclose(trj1['vel'], trj2['vel'])
            self.assertTrue(velok)
            forceok = np.allclose(trj1['force'], trj2['force'])
            self.assertTrue(forceok)


if __name__ == '__main__':
    unittest.main()
