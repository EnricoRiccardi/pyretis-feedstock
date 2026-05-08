"""Time the FORTRAN implementation of the Lennard-Jones potential.

This timing is simply done by evaluating the Lennard-Jones forces
(and potential) for different system sizes.
"""
# pylint: disable=invalid-name
import timeit
import numpy as np
from pyretis.core import System, create_box, Particles
from pyretis.core.units import create_conversion_factors
from pyretis.tools import generate_lattice
from ljpotentialfp import PairLennardJonesCutFp


def set_up_initial_state(nlattice=5):
    """Create particles for the test.

    This will set up a 3D lattice with 4*nlattice**3 particles.
    """
    create_conversion_factors('lj')
    lattice, size = generate_lattice('fcc', [nlattice] * 3, density=0.9)
    npart = len(lattice)
    size = np.array(size)
    lattice += np.random.randn(npart, 3) * 0.05
    box = create_box(low=size[:, 0], high=size[:, 1],
                     periodic=[True, True, True])
    sys = System(temperature=1.0, units='lj', box=box)
    sys.particles = Particles(dim=3)
    for pos in lattice:
        sys.add_particle(name='Ar', pos=pos, mass=1.0, ptype=0)
    print(f'Created lattice with {sys.particles.npart} atoms.')
    return sys


def test_wrapper(func, *args, **kwargs):
    """A simple wrapper for calling functions."""
    def wrapped():
        return func(*args, **kwargs)
    return wrapped


def test_function(function, system, repeat=3, number=5):
    """Run the test for a function."""
    print(f'Testing function: {function.__name__}')
    wrapped = test_wrapper(function, system)
    res = timeit.repeat(wrapped, repeat=repeat, number=number)
    best = min(res) / float(number)
    avg = np.average([resi / float(number) for resi in res])
    std = np.std([resi / float(number) for resi in res])
    print(f'Best: {best}')
    print(f'Average: {avg} +- {std}')
    return best, avg, std


def main():
    """Run the timings."""
    parameters = {0: {'sigma': 1.0, 'epsilon': 1.0, 'rcut': 2.5}}
    # set up potentials:
    potential = PairLennardJonesCutFp(dim=3, shift=True, mixing='geometric')
    potential.set_parameters(parameters)

    results = []

    for i in range(3, 16):
        system = set_up_initial_state(nlattice=i)
        print('Testing FORTRAN implementation')
        time1 = test_function(potential.potential,
                              system,
                              number=10, repeat=3)
        results.append((system.particles.npart, time1[0], time1[1], time1[2]))
    results = np.array(results)
    np.savetxt('timings.txt', results, fmt='%i %.9e %.9e %.9e',
               header='N best avg std')


if __name__ == '__main__':
    main()
