# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A 1D harmonic potential."""
from matplotlib import pyplot as plt
import numpy as np
from pyretis.core import System, Particles
from pyretis.forcefield import ForceField
from pyretis.forcefield.potential import PotentialFunction


class Harmonic1D(PotentialFunction):
    """A 1D harmonic potential function."""

    def __init__(self):
        """Set up the potential."""
        super().__init__(dim=1, desc='1D harmonic potential')
        self.eq_pos = 0.0  # equilibrium position
        self.k_force = 1.0  # force constant

    def potential(self, system):
        """Calculate potential energy."""
        pos = system.particles.pos  # Get positions from the particle list.
        vpot = 0.5 * self.k_force * (pos - self.eq_pos)**2
        return vpot.sum()

    def force(self, system):
        """Calculate the force."""
        pos = system.particles.pos  # Get positions from the particle list.
        virial = None  # We are lazy and do not calculate the virial here.
        forces = -self.k_force * (pos - self.eq_pos)
        return forces, virial

    def potential_and_force(self, system):
        """Calculate force and potential."""
        pot = self.potential(system)
        force, virial = self.force(system)
        return pot, force, virial


if __name__ == '__main__':
    # Create empty force field:
    forcefield = ForceField(desc='1D Harmonic potential')
    pot = Harmonic1D()  # create potential
    # Add the potential function to the force field:
    forcefield.add_potential(pot)
    system = System()
    system.particles = Particles(dim=1)
    system.add_particle(pos=np.zeros(1))
    # Do some plotting, first calculate the potential as some locations:
    vpotentials = []
    positions = np.linspace(-2.5, 2.5, 100)
    for xi in positions:
        system.particles.pos = xi
        vpotentials.append(forcefield.evaluate_potential(system))
    # And plot it using matplotlib:
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(positions, vpotentials, lw=3)
    plt.show()
