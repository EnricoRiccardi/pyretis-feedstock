# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""
Here we test the basic functionality of the CP2K engine.

1) We generate random velocities.

2) We integrate forward in time.

3) We reverse the velocities.

4) We integrate backward in time from 3.
"""
import os
import sys
import shutil
import time
import itertools
import colorama
import numpy as np
from matplotlib import pyplot as plt
from pyretis.core import System, create_box, ParticlesExt, Path
from pyretis.orderparameter import PositionVelocity
from pyretis.inout.common import make_dirs
from pyretis.inout.settings import parse_settings_file
from pyretis.inout.formats.xyz import read_xyz_file, write_xyz_trajectory
from pyretis.engines.cp2k import (
    CP2KEngine,
    write_for_step_vel,
    convert_snapshot
)
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


plt.style.use('seaborn-v0_8-deep')


def clean_dir(dirname):
    """Remove ALL files in the given directory."""
    for files in os.listdir(dirname):
        filename = os.path.join(dirname, files)
        if os.path.isfile(filename):
            os.remove(filename)


def run_in_steps(engine, system, order_parameter, interfaces,
                 steps=25, exe_dir='forward-step', reverse=False):
    """Run the engine forward in time, in steps.

    Parameters
    ----------
    engine : object like :py:class:`.ExternalMDEngine
        Engine to use for propagation.
    system : object like :py:class:`.System`
        The system we are propagation.
    order_parameter : object like :py:class:`.OrderParameter`
        An order parameter to calculate.
    interfaces : list of floats
        Interfaces to consider, here typically just set to
        ``[-float('inf'), float('inf'), float('inf')]``.
    steps : integer
        The number of steps we will do.
    exe_dir : string
        The folder to use for the execution.
    reverse : boolean
        Selects the time direction.

    """
    logger.info(f'\nRunning {steps} steps in "{exe_dir}"')
    logger.info(f'(Reverse = {reverse})')
    make_dirs(exe_dir)
    folder = os.path.abspath(exe_dir)
    clean_dir(folder)
    engine.exe_dir = folder
    path = Path(None, maxlen=steps)
    start = time.perf_counter()
    engine.propagate(path,
                     {'system': system, 'order_function': order_parameter,
                      'interfaces': interfaces},
                     reverse=reverse)
    end = time.perf_counter()
    logger.info('Propagation done!')
    logger.info(f'Time spent: {end - start}')
    return path


def run_plain_cp2k(engine, system, order_parameter, input_conf,
                   steps=25, exe_dir='forward-plain', reverse=False):
    """Run a plain CP2K simulation.

    Parameters
    ----------
    engine : object like :py:class:`.CP2KEngine`
        This is used to obtain the commands for executing CP2K.
    system : object like :py:class:`.System`
        The system we are propagation.
    order_parameter : object like :py:class:`.OrderParameter`
        An order parameter to calculate.
    input_conf : string
        The input configuration to use for the CP2K run.
    steps : integer
        The number of steps to run.
    exe_dir : string
        The path to where we should execute CP2K.
    reverse : boolean
        Selects the time direction.

    """
    logger.info(f'\nRunning {steps} plain CP2K steps in "{exe_dir}"')
    logger.info(f'(Reverse = {reverse})')
    make_dirs(exe_dir)
    folder = os.path.abspath(exe_dir)
    clean_dir(folder)
    engine.exe_dir = exe_dir
    shutil.copy(input_conf, folder)
    _, _, vel, _ = engine._read_configuration(input_conf)
    input_file = os.path.join(exe_dir, 'md.inp')
    write_for_step_vel(
        engine.input_files['template'],
        input_file,
        engine.timestep,
        (steps - 1) * engine.subcycles,
        os.path.basename(input_conf),
        vel,
        name='md',
        print_freq=1
    )
    engine.add_input_files(exe_dir)
    logger.info('Running CP2K...')
    start = time.perf_counter()
    engine.run_cp2k('md.inp', 'md')
    end = time.perf_counter()
    logger.info('Done!')
    logger.info(f'Time spent: {end - start}')
    energy_file = os.path.join(exe_dir, 'md-1.ener')
    energy = np.loadtxt(energy_file)
    pos_file = os.path.join(exe_dir, 'md-pos-1.xyz')
    vel_file = os.path.join(exe_dir, 'md-vel-1.xyz')
    order = []
    box = np.array([5.0, 5.0, 5.0])
    out_file = os.path.join(exe_dir, 'rev_last.xyz')
    out_file2 = os.path.join(exe_dir, 'rev_last2.xyz')
    for snapshotp, snapshotv in zip(read_xyz_file(pos_file),
                                    read_xyz_file(vel_file)):
        _, xyz, _, names = convert_snapshot(snapshotp)
        _, vel, _, _ = convert_snapshot(snapshotv)  # interpret pos as vel
        system.particles.pos = xyz
        system.particles.vel = vel
        if reverse:
            system.particles.vel *= -1
        system.box.update_size(box)
        order.append(order_parameter.calculate(system))
        write_xyz_trajectory(out_file, xyz, -vel, names, box,
                             append=False)
        write_xyz_trajectory(out_file2, xyz, -vel, names, box,
                             append=False)
    order = np.array(order)
    return energy, order, out_file


def test_genvel(engine, input_file, exe_dir='genvel'):
    """Test generation of velocities.

    Parameters
    ----------
    engine : object like :py:class:`.ExternalMDEngine
        Engine to use for the generation.
    input_file : string
        A file containing the input configuration for CP2K.
    exe_dir : string, optional
        The folder where we will be executing CP2K.

    Returns
    -------
    out_conf : string
        The path to the output configuration with randomized
        velocities.

    """
    logger.info(f'\nRunning CP2K genvel step in "{exe_dir}"')
    make_dirs(exe_dir)
    folder = os.path.abspath(exe_dir)
    clean_dir(folder)
    engine.exe_dir = folder
    start = time.perf_counter()
    out_conf, _ = engine._prepare_shooting_point(input_file)
    end = time.perf_counter()
    logger.info('Generation of velocity done!')
    logger.info(f'Time spent: {end - start}')
    return out_conf


def main(plot=False):
    """Run the tests."""
    settings = parse_settings_file('engine.rst')
    steps = settings['simulation']['steps']
    engine_settings = settings['engine']
    engine = CP2KEngine(
        engine_settings['cp2k'],
        engine_settings['input_path'],
        engine_settings['timestep'],
        engine_settings['subcycles'],
        extra_files=engine_settings.get('extra_files', []),
        conf=engine_settings.get('conf', 'h2.xyz')
    )
    logger.info(f'Testing engine: {engine}')
    logger.info(f'Time step: {engine.timestep}')
    logger.info(f'Subcycles: {engine.subcycles}')
    system = System(units='gromacs',
                    box=create_box(cell=[100, 100, 100]),
                    temperature=500)
    system.particles = ParticlesExt(dim=3)
    initial_conf = engine.input_files['conf']
    interfaces = [-float('inf'), float('inf'), float('inf')]
    order_parameter = PositionVelocity(0, dim='x', periodic=True)

    # First generate some random velocities:
    initial_conf = test_genvel(engine, initial_conf, exe_dir='genvel')
    system.particles.set_pos((initial_conf, None))
    system.particles.set_vel(False)

    pathf = run_in_steps(engine, system, order_parameter, interfaces,
                         steps=steps, exe_dir='forward', reverse=False)

    # set state to last point in trajectory:
    phase_point = pathf.phasepoints[-1]
    system = phase_point.copy()

    pathb = run_in_steps(engine, system, order_parameter, interfaces,
                         steps=steps, exe_dir='backward', reverse=True)

    plainf = run_plain_cp2k(engine, system, order_parameter, initial_conf,
                            steps=steps, exe_dir='forward-plain',
                            reverse=False)

    plainb = run_plain_cp2k(engine, system, order_parameter, plainf[-1],
                            steps=steps, exe_dir='backward-plain',
                            reverse=True)
    obtain_mses(pathf, pathb, plainf, plainb, engine.subcycles)
    if plot:
        steps = np.array([i * engine.subcycles for i in range(steps)])
        plot_path_comparison(steps, pathf, pathb, plainf, plainb,
                             engine.subcycles)


def mse_combinations(text, var, tol=None):
    """Calculate mse for several combinations."""
    for comb in itertools.combinations(var, 2):
        mse = ((comb[0][0] - comb[1][0])**2).mean(axis=0)
        level = 'info'
        tol_ok = True
        if tol:
            try:
                tol_ok = all([abs(i) < tol for i in mse])
            except TypeError:
                tol_ok = abs(mse) < tol
            if not tol_ok:
                level = 'error'
        if level == 'error':
            logger.error(f'MSE {text}: {comb[0][1]} vs {comb[1][1]} = {mse}')
        else:
            logger.info(f'MSE {text}: {comb[0][1]} vs {comb[1][1]} = {mse}')
        if not tol_ok:
            sys.exit(1)


def obtain_mses(pathf, pathb, plainf, plainb, subcycles):
    """Calculate mean square errors."""
    orderf = [i.order for i in pathf.phasepoints]
    orderb = [i.order for i in reversed(pathb.phasepoints)]
    mses = [(np.array(orderf), 'step-forward'),
            (np.array(orderb), 'step-back'),
            (plainf[1][::subcycles, :], 'plain-forward'),
            (plainb[1][::subcycles, :][::-1], 'plain-back')]
    mse_combinations('order parameters', mses, tol=1e-7)

    ekinf = [i.particles.ekin for i in pathf.phasepoints]
    ekinb = [i.particles.ekin for i in reversed(pathb.phasepoints)]
    mses = [(np.array(ekinf), 'step-forward'),
            (np.array(ekinb), 'step-back'),
            (plainf[0][::subcycles, 2], 'plain-forward'),
            (plainb[0][::subcycles, 2][::-1], 'plain-back')]
    mse_combinations('kinetic energy', mses, tol=1e-7)

    vpotf = [i.particles.vpot for i in pathf.phasepoints]
    vpotb = [i.particles.vpot for i in reversed(pathb.phasepoints)]
    mses = [(np.array(vpotf), 'step-forward'),
            (np.array(vpotb), 'step-back'),
            (plainf[0][::subcycles, 4], 'plain-forward'),
            (plainb[0][::subcycles, 4][::-1], 'plain-back')]
    mse_combinations('potential energy', mses, tol=1e-4)


def plot_path_comparison(steps, pathf, pathb, plainf, plainb, subcycles):
    """Just plot some properties for the paths."""
    orderf = np.array([i.order for i in pathf.phasepoints])
    orderb = np.array([i.order for i in reversed(pathb.phasepoints)])
    fig1 = plt.figure(figsize=(12, 6))
    ax11 = fig1.add_subplot(121)
    ax12 = fig1.add_subplot(122)
    ax11.plot(steps, orderf[:, 0], lw=2, ls='-', marker='o',
              label='Forward', ms=12, alpha=0.8)
    ax11.plot(steps, orderb[:, 0], lw=2, ls='--', marker='^',
              label='Backward', ms=12, alpha=0.8)
    ax11.plot(plainf[1][:, 0], lw=2, ls='--', marker='^',
              label='Plain forward', alpha=0.8)
    ax11.plot(plainb[1][:, 0][::-1], lw=2, ls='-.', marker='x',
              label='Plain backward', alpha=0.8)
    ax11.legend()
    ax11.set_title('Order param 1')
    ax12.plot(steps, orderf[:, 1], lw=2, ls='-', marker='o', ms=12,
              alpha=0.8)
    ax12.plot(steps, orderb[:, 1], lw=2, ls='--', marker='^', ms=12,
              alpha=0.8)
    ax12.plot(plainf[1][:, 1], lw=2, ls=':', marker='s', alpha=0.8)
    ax12.plot(plainb[1][:, 1][::-1], lw=2, ls='-.', marker='x',
              alpha=0.8)
    ax12.set_title('Order param 2')

    ekinf = np.array([i.particles.ekin for i in pathf.phasepoints])
    ekinb = np.array(
        [i.particles.ekin for i in reversed(pathb.phasepoints)]
    )
    fig2 = plt.figure(figsize=(12, 6))
    ax21 = fig2.add_subplot(121)
    ax22 = fig2.add_subplot(122)
    ax21.plot(steps, ekinf, lw=2, ls='-', marker='o', label='Forward', ms=12,
              alpha=0.8)
    ax21.plot(steps, ekinb, lw=2, ls='--', marker='^', label='Backward',
              ms=9, alpha=0.8)
    ax21.plot(plainf[0][:, 0], plainf[0][:, 2], lw=2, ls=':', marker='s',
              label='Plain forward', alpha=0.8)
    ax21.plot(plainb[0][:, 0], plainb[0][:, 2][::-1], lw=2, ls='-.',
              marker='x', label='Plain backward', alpha=0.8)
    ax21.set_title('Kinetic energy')
    ax21.legend()

    vpotf = np.array([i.particles.vpot for i in pathf.phasepoints])
    vpotb = np.array(
        [i.particles.vpot for i in reversed(pathb.phasepoints)]
    )
    ax22.plot(steps, vpotf, lw=2, ls='-', marker='o', label='Forward',
              ms=12, alpha=0.8)
    ax22.plot(steps, vpotb, lw=2, ls='--', marker='^', label='Backward',
              ms=12, alpha=0.8)
    ax22.plot(plainf[0][:, 0], plainf[0][:, 4], lw=2, ls=':', marker='s',
              label='Plain forward', alpha=0.8)
    ax22.plot(plainb[0][:, 0], plainb[0][:, 4][::-1], lw=2, ls='-.',
              marker='x', label='Plain backward', alpha=0.8)
    ax22.set_title('Potential energy')
    ax22.legend()

    fig3 = plt.figure(figsize=(12, 6))
    ax31 = fig3.add_subplot(221)
    ax32 = fig3.add_subplot(222)
    ax33 = fig3.add_subplot(223)
    ax34 = fig3.add_subplot(224)

    ax31.scatter(orderf[:, 0], orderb[:, 0], marker='o', label='Backward',
                 alpha=0.8)
    ax31.scatter(orderf[:, 0], plainf[1][::subcycles, 0], marker='s',
                 label='Forward-plain', alpha=0.8)
    ax31.scatter(orderf[:, 0], plainb[1][::subcycles, 0][::-1], marker='^',
                 label='Backward-plain', alpha=0.8)
    minx = min(orderf[:, 0])
    maxx = max(orderf[:, 0])
    ax31.plot([minx, maxx], [minx, maxx], ls=':',
              c='#262626', alpha=0.5, lw=2)
    ax31.set_xlabel('Order 1 Forward')
    ax31.set_ylabel('Order 1')
    ax31.legend()

    ax32.scatter(orderf[:, 1], orderb[:, 1], marker='o', label='Backward',
                 alpha=0.8)
    ax32.scatter(orderf[:, 1], plainf[1][::subcycles, 1], marker='s',
                 label='Forward-plain', alpha=0.8)
    ax32.scatter(orderf[:, 1], plainb[1][::subcycles, 1][::-1], marker='^',
                 label='Backward-plain', alpha=0.8)
    ax32.set_xlabel('Order 2 Forward')
    ax32.set_ylabel('Order 2')
    minx = min(orderf[:, 1])
    maxx = max(orderf[:, 1])
    ax32.plot([minx, maxx], [minx, maxx], ls=':',
              c='#262626', alpha=0.5, lw=2)
    ax32.legend()

    ax33.scatter(ekinf, ekinb, marker='o', label='Backward',
                 alpha=0.8)
    ax33.scatter(ekinf, plainf[0][::subcycles, 2], marker='s',
                 label='Forward-plain', alpha=0.8)
    ax33.scatter(ekinf, plainb[0][::subcycles, 2][::-1], marker='^',
                 label='Backward-plain', alpha=0.8)
    ax33.set_xlabel('Ekin Forward')
    ax33.set_ylabel('Ekin')
    minx = min(ekinf)
    maxx = max(ekinf)
    ax33.plot([minx, maxx], [minx, maxx], ls=':',
              c='#262626', alpha=0.5, lw=2)
    ax33.legend()

    ax34.scatter(vpotf, vpotb, marker='o', label='Backward',
                 alpha=0.8)
    ax34.scatter(vpotf, plainf[0][::subcycles, 4], marker='s',
                 label='Forward-plain', alpha=0.8)
    ax34.scatter(vpotf, plainb[0][::subcycles, 4][::-1], marker='^',
                 label='Backward-plain', alpha=0.8)
    ax34.set_xlabel('Vpot Forward')
    ax34.set_ylabel('Vpot')
    minx = min(vpotf)
    maxx = max(vpotf)
    ax34.plot([minx, maxx], [minx, maxx], ls=':',
              c='#262626', alpha=0.5, lw=2)
    ax34.legend()

    fig1.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()

    plt.show()


if __name__ == '__main__':
    colorama.init(autoreset=True)
    PLOT = len(sys.argv) > 1
    main(plot=PLOT)
