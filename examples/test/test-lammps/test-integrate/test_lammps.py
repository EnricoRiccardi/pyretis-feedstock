# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A test using the LAMMPS engine."""
import os
import sys
import colorama
from matplotlib import pyplot as plt
from pyretis.engines.lammps import LAMMPSEngine
from pyretis.inout.common import make_dirs
from pyretis.testing.systemhelp import create_system_ext
from pyretis.testing.helpers import clean_dir
from plotting import plot_compare, plot_xy
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


HERE = os.path.abspath(os.path.dirname(__file__))
STEPS = 50
SUBCYCLES = 2


def run_in_folder(system, engine, exe_dir, steps):
    """Run the engine in the given folder."""
    make_dirs(exe_dir)
    clean_dir(exe_dir)
    engine.exe_dir = exe_dir
    engine.integrate(system, steps, order_function=None, thermo='full')
    energy = engine.read_energies('pyretis_md')
    order = engine.read_order_parameters('order_pyretis_md.txt')
    return energy, order


def run_forward():
    """Use the LAMMPS engine to run a MD simulation forward in time."""
    logger.info('\nTesting the LAMMPS engine by running forward.\n',)
    engine = LAMMPSEngine('lmp_serial', 'lammps_input', SUBCYCLES)
    logger.info(f'Running forward for {2 * STEPS} steps...')
    # Create a dummy system:
    system = create_system_ext(pos=('system.data', 0))
    exe_dir = os.path.join(HERE, 'run-forward1')
    energy_f1, order_f1 = run_in_folder(system, engine, exe_dir, 2*STEPS)
    # Run forward, for the initial state for 100 steps:
    logger.info('-> Forward done! Resetting system.')
    logger.info(f'\nRunning forward for {STEPS} steps...')
    system = create_system_ext(pos=('system.data', 0))
    exe_dir = os.path.join(HERE, 'run-forward2')
    energy_f2, order_f2 = run_in_folder(system, engine, exe_dir, STEPS)
    # Run for 50 more steps:
    logger.info('-> Forward done!')
    logger.info(f'\nContinuing run for {STEPS} steps...')
    config = system.particles.get_pos()
    traj_file = os.path.join(exe_dir, config[0])
    system.particles.set_pos((traj_file, config[1]))
    exe_dir = os.path.join(HERE, 'run-forward3')
    energy_f3, order_f3 = run_in_folder(system, engine, exe_dir, STEPS)
    # Update step for part 3:
    energy_f3['Step'] += energy_f2['Step'][-1]
    logger.info('-> Forward done! Comparing energies:')

    data_sets = []

    for key in ['PotEng', 'KinEng', 'Press', 'Temp']:
        data_set = []
        for i, energy in enumerate((energy_f1, energy_f2, energy_f3)):
            if i == 0:
                lab = 'Full simulation'
            else:
                lab = f'Simulation, part {i}'
            data_set.append(
                (energy['Step'], energy[key], lab)
            )
        data_sets.append(data_set)

    plot_compare(data_sets, ['PotEng', 'KinEng', 'Press', 'Temp'])

    data_sets_xy = []
    for key in ('PotEng', 'KinEng', 'Press', 'Temp'):
        energy_parts = [i for i in energy_f2[key]]
        # Remove the data point that is common to part 1 and part 2:
        energy_parts.pop()
        # Add data points from part 2:
        energy_parts += [i for i in energy_f3[key]]
        data_sets_xy.append(
            (
                energy_f1[key], energy_parts,
                f'{key}, full simulation',
                f'{key}, part1 & part2'
            )
        )
    plot_xy(data_sets_xy)

    logger.info('-> Will also compare order parameters:')
    data_sets = [
        [
            (energy_f1['Step'], order_f1[:, 0], 'Full simulation'),
            (energy_f2['Step'], order_f2[:, 0], 'Simulation, part 1'),
            (energy_f3['Step'], order_f3[:, 0], 'Simulation, part 2'),
        ],
        [
            (energy_f1['Step'], order_f1[:, 1], 'Full simulation'),
            (energy_f2['Step'], order_f2[:, 1], 'Simulation, part 1'),
            (energy_f3['Step'], order_f3[:, 1], 'Simulation, part 2'),
        ],
    ]
    plot_compare(data_sets, ['OP1', 'OP2'])

    data_sets_xy = []
    for i in (0, 1):
        order_parts = [i for i in order_f2[:, i]]
        # Remove the data point that is common to part 1 and part 2:
        order_parts.pop()
        # Add data points from part 2:
        order_parts += [i for i in order_f3[:, i]]
        data_sets_xy.append(
            (
                order_f1[:, i], order_parts,
                f'OP{i + 1}, full simulation',
                f'OP{i + 1}, part1 & part2'
            )
        )
    plot_xy(data_sets_xy)


def run_forward_backward():
    """Use the LAMMPS engine to run a MD simulation forward in time."""
    logger.info(
        '\nTesting the LAMMPS engine by running forward and backward.\n')
    # Set up for running forward:
    logger.info(f'-> Running forward for {STEPS} steps...')
    engine = LAMMPSEngine('lmp_serial', 'lammps_input', SUBCYCLES)
    system = create_system_ext(pos=('system.data', 0))
    exe_dir = os.path.join(HERE, 'run-forward')
    energy_f, order_f = run_in_folder(system, engine, exe_dir, STEPS)
    logger.info('-> Forward done! Prepare backward run.')
    # Reverse the velocities and run:
    config = system.particles.get_pos()
    traj_file = os.path.join(exe_dir, config[0])
    logger.info('-> Setting configuration.')
    system.particles.set_pos((traj_file, config[1]))
    logger.info('-> Reversing velocities.')
    system.particles.set_vel(True)

    logger.info(f'\n-> Starting backward run ({STEPS}) steps!')
    exe_dir = os.path.join(HERE, 'run-backward')
    energy_b, order_b = run_in_folder(system, engine, exe_dir, STEPS)
    logger.info('\n-> Backward done! Will compare energies:')

    step = energy_b['Step'] - energy_b['Step'][0]
    data_sets = [
        [
            (energy_f['Step'], energy_f['PotEng'],
             'Forward'),
            (step, energy_b['PotEng'][::-1], 'Backward'),
        ],
        [
            (energy_f['Step'], energy_f['KinEng'],
             'Forward'),
            (step, energy_b['KinEng'][::-1], 'Backward'),
        ],
        [
            (energy_f['Step'], energy_f['Press'],
             'Forward'),
            (step, energy_b['Press'][::-1], 'Backward'),
        ],
        [
            (energy_f['Step'], energy_f['Temp'],
             'Forward'),
            (step, energy_b['Temp'][::-1], 'Backward'),
        ],
    ]
    plot_compare(data_sets, ['PotEng', 'KinEng', 'Press', 'Temp'])

    data_sets_xy = [
        (energy_f['PotEng'], energy_b['PotEng'][::-1],
         'PotEng, forward', 'PotEng, backward'),
        (energy_f['KinEng'], energy_b['KinEng'][::-1],
         'KinEng, forward', 'KinEng, backward'),
        (energy_f['Press'], energy_b['Press'][::-1],
         'Press, forward', 'Press, backward'),
        (energy_f['Temp'], energy_b['Temp'][::-1],
         'Temp, forward', 'Temp, backward'),
    ]
    plot_xy(data_sets_xy)

    logger.info('\n-> Will also compare order parameters:')
    logger.warning('NOTE: OP2 IS HERE MULTIPLIED WITH "-1"')
    logger.warning('NOTE: THIS, SINCE IT IS VELOCITY DEPENDENT!')
    data_sets = [
        [
            (step, order_f[:, 0], 'Forward'),
            (step, order_b[:, 0][::-1], 'Backward'),
        ],
        [
            (step, order_f[:, 1], 'Forward'),
            (step, -1 * order_b[:, 1][::-1], 'Backward'),
        ],
    ]
    plot_compare(data_sets, ['OP1', 'OP2'])

    data_sets_xy = [
        (order_f[:, 0], order_b[:, 0][::-1], 'OP1, forward', 'OP1, backward'),
        (order_f[:, 1], -1 * order_b[:, 1][::-1],
         'OP2, forward', 'OP2, backward'),
    ]
    plot_xy(data_sets_xy)


def main(plot=False):
    """Run the comparisons."""
    run_forward()
    run_forward_backward()
    if plot:
        plt.show()


if __name__ == '__main__':
    colorama.init(autoreset=True)
    main(plot=len(sys.argv) >= 2)
