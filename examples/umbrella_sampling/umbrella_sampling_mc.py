# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Example showing a simple umbrella simulation with PyRETIS.

In this simulation, we study a particle moving in a one-dimensional
potential energy landscape and the goal is to determine this
landscape by performing umbrella simulations.

"""

import sys
import colorama
from tqdm import tqdm
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.cm import get_cmap
from pyretis.core import RandomGenerator, System, Particles
from pyretis.setup import create_simulation, create_force_field
from pyretis.analysis import histogram, match_all_histograms
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


UMBRELLA_WINDOWS = [
    [-1.0, -0.4],
    [-0.5, -0.2],
    [-0.3, 0.0],
    [-0.1, 0.2],
    [0.1, 0.4],
    [0.3, 0.6],
    [0.5, 1.0]
]


DOUBLE_WELL = {
    'class': 'DoubleWell',
    'parameter': {'a': 1.0, 'b': 1.0, 'c': 0.02},
}

DEFAULT_SETTINGS = {
    'system': {
        'dimensions': 1,
        'units': 'metal',
        'temperature': 500,
    },
    'simulation': {
        'task': 'umbrellawindow',
        'rgen': 0,
        'mincycle': 10000,
        'maxdx': 0.1,
        'overlap': 1.0,
        'umbrella': 1,
    },
    'particles': {
        'position': {'input_file': 'position.xyz'},
        'name': ['X'],
        'ptype': [0],
    },
}

plt.style.use('seaborn')


def set_up_simulation(umbrella, over, seed, initial_positions=None):
    """Set up a single umbrella window simulation.

    Parameters
    ----------
    umbrella : list of floats
        The umbrella window we are investigating.
    over : float
        The coordinate we need to cross for the
        current window.
    seed : integer
        A seed for the random number generator.
    initial_positions : numpy.array, optional
        The initial position(s) for the particle(s). If this
        is not provided the default settings are used.

    """
    settings = DEFAULT_SETTINGS.copy()
    settings['simulation']['overlap'] = over
    settings['simulation']['umbrella'] = umbrella
    settings['simulation']['rgen'] = RandomGenerator(seed=seed)
    settings['potential'] = [
        DOUBLE_WELL,
        {
            'class': 'RectangularWell',
            'parameter': {'left': umbrella[0], 'right': umbrella[1]},
        }
    ]
    # Create the umbrella simulation:
    simulation = create_simulation(settings)
    if initial_positions is not None:
        simulation.system.particles.pos = np.copy(initial_positions)
    # Make sure initial potential energy is set:
    simulation.system.potential()
    return simulation


def run_umbrellas(windows):
    """Run the sampling for the given set of windows.

    Parameters
    ----------
    windows : list of floats
        The umbrella windows to investigate.

    """
    msg = '\nRunning umbrella no: {} of {}. Location: {}'
    n_umb = len(windows)
    logger.info('Starting simulations:')
    trajectories = []
    energies = []
    system = None
    for i, window in enumerate(windows):
        logger.info(msg.format(i + 1, n_umb, window))
        over = windows[min(i + 1, n_umb - 1)][0]
        initial_positions = None if system is None else system.particles.pos
        simulation = set_up_simulation(
            window,
            over,
            1,
            initial_positions=initial_positions,
        )
        system = simulation.system
        traj, ener = [], []
        for _ in tqdm(simulation.run()):
            for pos in system.particles.pos:
                traj.append(pos)
                ener.append(system.particles.vpot)
        trajectories.append(traj)
        energies.append(ener)
        nstep = simulation.cycle['step'] - simulation.cycle['startcycle']
        logger.info(f'Done. Cycles: {nstep}')
    return system, trajectories, energies


def analysis_and_plot(system, trajectory, windows):
    """Plot some results from the simulation."""
    # We can now post-process the simulation output.
    bins = 100
    lim = (-1.0, 1.0)
    histograms = [histogram(trj, bins=bins, limits=lim) for trj in trajectory]
    # Extract the bins (the midpoints) and the bin-width:
    bin_x = histograms[0][-1]
    dbin = bin_x[1] - bin_x[0]
    # We are going to match these histograms:
    logger.info('Matching histograms...')
    histograms_s, _, hist_avg = match_all_histograms(histograms, windows)
    # Let us create some simple plots using matplotlib:
    plot_histograms(histograms_s, hist_avg, bin_x, dbin)
    plot_free_energy(hist_avg, bin_x, system)


def plot_histograms(histograms_s, hist_avg, bin_x, dbin):
    """Plot matched histograms.

    Parameters
    ----------
    histograms_s : list of numpy.arrays
        The scaled histograms.
    hist_avg : numpy.array
        The average histogram.
    bin_x : numpy.array
        The midpoints for the histograms.
    dbin : float
        The histogram width.

    """
    logger.info('Plotting matched histograms')
    fig = plt.figure()
    axs = fig.add_subplot(111)
    axs.set_yscale('log')
    axs.set_xlabel('Position ($x$)', fontsize='large')
    axs.set_ylabel('Matched histograms', fontsize='large')
    colors = []
    for maps in ('viridis', 'Spectral', None):
        try:
            cmap = get_cmap(name=maps)
            colors = cmap(np.linspace(0, 1, len(histograms_s)))
            break
        except ValueError:
            continue

    for i, histo in enumerate(histograms_s):
        axs.bar(bin_x - 0.5 * dbin, histo, dbin, color=colors[i],
                alpha=0.8, log=True, edgecolor='#262626')
    axs.plot(bin_x, hist_avg, lw=7, color='orangered', alpha=0.8,
             label='Average after matching')
    axs.legend()
    axs.set_xlim((-1.1, 1.1))
    axs.set_ylim((0.1, hist_avg.max() * 1.25))
    fig.tight_layout()


def calculate_unbiased_potential(xpos):
    """Calculate the unbiased potential for the given locations."""
    # Set up and plot the unbiased potential:
    unbiased = {
        'forcefield': {
            'description': 'Double well',
        },
        'potential': [DOUBLE_WELL]
    }
    forcefield = create_force_field(unbiased)
    # Create a fake system for calculating the potential:
    system = System()
    system.particles = Particles()
    vpot = []
    for i in xpos:
        system.particles.pos = i
        vpot.append(forcefield.evaluate_potential(system))
    return np.array(vpot)


def plot_free_energy(hist_avg, bin_x, system):
    """Plot the free energy obtained.

    Parameters
    ----------
    hist_avg : numpy.array
        The average histogram.
    bin_x : numpy.array
        Midpoints for the bins.
    system : object like :py:class:`.System`
        The system we have been investigating, we are here
        using it to plot the unbiased potential we have been
        sampling.

    """
    logger.info('Plotting the free energy')
    fig = plt.figure()
    axi = fig.add_subplot(111)
    xpos = np.linspace(-2, 2, 1000)
    free = -np.log(hist_avg) / system.temperature['beta']  # Free energy.
    vpot = calculate_unbiased_potential(xpos)
    free += (vpot.min() - free.min())

    axi.plot(xpos, vpot, lw=3, label='Unbiased potential', alpha=0.5)
    axi.plot(bin_x, free, lw=7, alpha=0.5,
             label='Free energy from umbrella simulations')
    axi.set_xlabel('Position ($x$)', fontsize='large')
    axi.set_ylabel('Potential energy ($V(x)$) / eV', fontsize='large')
    axi.legend()
    axi.set_xlim((-1.1, 1.1))
    axi.set_ylim((-0.3, 0.05))
    fig.tight_layout()


if __name__ == '__main__':
    colorama.init(autoreset=True)
    SYS, TRAJ, _ = run_umbrellas(UMBRELLA_WINDOWS)
    analysis_and_plot(SYS, TRAJ, UMBRELLA_WINDOWS)
    if 'noplot' not in sys.argv[1:]:
        plt.show()
