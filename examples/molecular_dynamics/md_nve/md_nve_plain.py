# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Example of running a MD simulation using the PyRETIS library.

The system considered is a simple Lennard-Jones fluid.

"""
import sys
import numpy as np
from matplotlib import pyplot as plt
from pyretis.core import create_box, Particles, System
from pyretis.simulation import SimulationNVE
from pyretis.engines import VelocityVerlet
from pyretis.tools import generate_lattice
from pyretis.forcefield import ForceField
from pyretis.forcefield.potentials import PairLennardJonesCutnp


# If you have a recent version of matplotlib you can make
# the plot nicer by loading a style, e.g.:
plt.style.use('seaborn-colorblind')


def set_up_system():
    """Set up a system using the PyRETIS library."""
    print('Creating box:')
    xyz, size = generate_lattice('fcc', [3, 3, 3], density=0.9)
    size = np.array(size)
    box = create_box(low=size[:, 0], high=size[:, 1])
    print(box)

    print('Creating system:')
    system = System(units='lj', box=box, temperature=2.0)
    system.particles = Particles(dim=3)
    for pos in xyz:
        system.add_particle(pos, vel=np.zeros_like(pos),
                            force=np.zeros_like(pos),
                            mass=1.0, name='Ar', ptype=0)

    gen_settings = {'distribution': 'maxwell', 'momentum': True}
    system.generate_velocities(**gen_settings)
    print(system.particles)

    print('Creating force field:')
    potentials = [PairLennardJonesCutnp(dim=3, shift=True, mixing='geometric')]
    parameters = [{0: {'sigma': 1, 'epsilon': 1, 'rcut': 2.5}}]
    ffield = ForceField('Lennard Jones force field',
                        potential=potentials, params=parameters)
    system.forcefield = ffield
    print(system.forcefield)
    return system


def set_up_simulation(system):
    """Set up the simulation."""
    print('Creating simulation:')
    engine = VelocityVerlet(0.002)
    simulation = SimulationNVE(ensemble={'system': system,
                                         'engine': engine},
                               settings={},
                               controls={'steps': 200})
    return simulation


def run_simulation(simulation):
    """Run the simulation and collect some outputs."""
    ekin = []
    vpot = []
    etot = []
    step = []
    for result in simulation.run():
        if result['cycle']['step'] % 10 == 0:
            print('Step:', result['cycle']['step'])
        step.append(result['cycle']['step'])
        ekin.append(result['thermo']['ekin'])
        vpot.append(result['thermo']['vpot'])
        etot.append(result['thermo']['etot'])
    return step, vpot, ekin, etot


def plot_results(step, vpot, ekin, etot, show=True):
    """Plot some energies."""
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(step, ekin, label='Kinetic energy')
    ax1.plot(step, etot, label='Total energy')
    ax1.plot(step, vpot, label='Potential energy')
    ax1.set_xlabel('Step no.')
    ax1.set_ylabel('Energy / reduced units')
    ax1.legend()
    fig1.tight_layout()
    fig1.savefig('out.png', dpi=150)
    if show:
        plt.show()


def main(args):
    """Set up, run simulation and plot some energies."""
    system = set_up_system()
    simulation = set_up_simulation(system)
    step, vpot, ekin, etot = run_simulation(simulation)
    plot_results(step, vpot, ekin, etot, show='noplot' not in args)


if __name__ == '__main__':
    main(sys.argv[1:])
