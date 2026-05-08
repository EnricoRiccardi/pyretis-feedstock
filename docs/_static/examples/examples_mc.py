# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This is an example of running an Umbrella Window simulation."""
import numpy as np
from matplotlib import pyplot as plt
from pyretis.core import System, RandomGenerator, Particles
from pyretis.simulation.mc_simulation import UmbrellaWindowSimulation
from pyretis.forcefield import ForceField
from pyretis.forcefield.potentials import DoubleWell, RectangularWell
from pyretis.analysis import histogram, match_all_histograms


mysystem = System(temperature=500, units='eV/K')
mysystem.particles = Particles(dim=1)
mysystem.add_particle(name='X', pos=np.array([-0.7]))

potential_dw = DoubleWell(a=1, b=1, c=0.02)
potential_rw = RectangularWell()
forcefield = ForceField(desc='Double well', potential=[potential_dw])
forcefield_bias = ForceField(desc='Double well with rectangular bias',
                             potential=[potential_dw, potential_rw])
mysystem.forcefield = forcefield_bias

umbrellas = [[-1.0, -0.4], [-0.5, -0.2], [-0.3, 0.0],
             [-0.1, 0.2], [0.1, 0.4], [0.3, 0.6], [0.5, 1.0]]
n_umb = len(umbrellas)
MINCYCLES = 1e4  # Number of MC steps to perform.
MAXDX = 0.1  # Maximum allowed displacement in the MC step(s).

RANDSEED = 1  # Seed for the random number generator:
RGEN = RandomGenerator(seed=RANDSEED)

trajectory, energy = [], []  # For storing trajectories & energies.
for i, umbrella in enumerate(umbrellas):
    print('Running umbrealla no: {} of {}. Location: {}'.format(i + 1, n_umb,
                                                                umbrella))
    # Move rectangular potential to correct place:
    params = {'left': umbrella[0], 'right': umbrella[1]}
    mysystem.forcefield.update_potential_parameters(potential_rw, params)
    mysystem.potential()  # Re-calculate potential energy.
    over = umbrellas[min(i + 1, n_umb - 1)][0]  # Position we must cross.
    simulation = UmbrellaWindowSimulation(
        ensemble={
            'system': mysystem,
            'rgen': RGEN
        },
        settings={
            'simulation': {
                'task': 'umbrella',
                'umbrella': umbrella,
                'overlap': over,
                'maxdx': MAXDX
            }
        },
        controls={'steps': MINCYCLES}
    )
    # Also create empy list for storing some data:
    traj, ener = [], []
    for result in simulation.run():
        for pos in mysystem.particles.pos:
            traj.append(pos)
            ener.append(mysystem.particles.vpot)
    trajectory.append(np.array(traj))
    energy.append(np.array(ener))
    print('Done. Cycles: {}'.format(simulation.cycle['stepno']))

# We can now post-process the simulation output:
BINS = 100
LIM = (-1.1, 1.1)
histograms = [histogram(traj, bins=BINS, limits=LIM) for traj in trajectory]
# Extract the bins (the midpoints) and the bin-width:
bin_x = histograms[0][-1]
dbin = bin_x[1] - bin_x[0]
# We are going to match these histograms:
print('Matching histograms...')
histograms_s, _, hist_avg = match_all_histograms(histograms, umbrellas)

print('Plotting matched histograms')
fig = plt.figure()
axs = fig.add_subplot(111)
axs.set_yscale('log')
axs.set_xlabel('Position ($x$)', fontsize='large')
axs.set_ylabel('Matched histograms', fontsize='large')
colors = ['blue', 'green', 'darkviolet', 'brown', 'gray', 'crimson', 'cyan']
for i, histo in enumerate(histograms_s):
    axs.bar(bin_x - 0.5 * dbin, histo, dbin, color=colors[i],
            alpha=0.6, log=True)
axs.plot(bin_x, hist_avg, lw=7, color='orangered', alpha=0.6,
         label='Average after matching')
axs.legend()
plt.xlim((-1.1, 1.1))

print('Plotting the free energy')
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
XPOT = np.linspace(-2, 2, 1000)
free = -np.log(hist_avg) / mysystem.temperature['beta']  # The free energy.
# Set up unbiased potential for plotting:
VPOT = []
for xi in XPOT:
    mysystem.particles.pos = xi
    VPOT.append(forcefield.evaluate_potential(mysystem))
VPOT = np.array(VPOT)
free += (VPOT.min() - free.min())
ax2.plot(XPOT, VPOT, 'blue', lw=3, label='Unbiased potential', alpha=0.5)
ax2.plot(bin_x, free, lw=7, alpha=0.5, color='green', label='Free energy')
ax2.set_xlabel('Position ($x$)', fontsize='large')
ax2.set_ylabel('Potential energy ($V(x)$) / eV', fontsize='large')
ax2.legend()
plt.xlim((-1.1, 1.1))
plt.ylim((-0.3, 0.05))
plt.show()
