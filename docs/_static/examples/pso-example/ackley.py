# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""The Ackley function, used as a potential."""
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # pylint: disable=unused-import
from pyretis.forcefield import PotentialFunction


TWO_PI = np.pi * 2.0
EXP = np.exp(1)


@np.vectorize
def ackley_potential(x, y):  # pylint: disable=invalid-name
    """Evaluate the Ackley function."""
    return (-20.0 * np.exp(-0.2*np.sqrt(0.5*(x**2 + y**2))) -
            np.exp(0.5 * (np.cos(TWO_PI * x) + np.cos(TWO_PI * y))) +
            EXP + 20)


class Ackley(PotentialFunction):
    """A implementation of the Ackley function.

    Note that the usage of this potential function differs from
    the usual usage for force fields.

    """

    def __init__(self):
        """Set up the function."""
        super().__init__(dim=2, desc='The Ackley function')

    def potential(self, system):
        """Evaluate the potential, note that we return all values."""
        xpos = system.particles.pos[:, 0]
        ypos = system.particles.pos[:, 1]
        pot = ackley_potential(xpos, ypos)
        return pot


def main():
    """Plot the Ackley function."""
    xgrid, ygrid = np.meshgrid(np.linspace(-5, 5, 100),
                               np.linspace(-5, 5, 100))
    zgrid = ackley_potential(xgrid, ygrid)

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, projection='3d')
    ax1.contourf(xgrid, ygrid, zgrid)
    ax2.plot_surface(xgrid, ygrid, zgrid, cmap=plt.get_cmap('viridis'))
    plt.show()


if __name__ == '__main__':
    main()
