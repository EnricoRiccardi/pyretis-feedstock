# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""
Here we test the basic functionality of the GROMACS engine.

1) We generate random velocities.

2) We integrate forward in time.

3) We reverse the velocities.

4) We integrate backward in time from 3.
"""
import itertools
import os
import shutil
import sys
import time
import colorama
from matplotlib import pyplot as plt
import numpy as np
from pyretis.core import System, create_box, ParticlesExt, Path
from pyretis.orderparameter.orderparameter import PositionVelocity
from pyretis.inout.common import make_dirs
from pyretis.inout.settings import parse_settings_file
from pyretis.inout.formats.gromacs import read_trr_file
from pyretis.engines import GromacsEngine, GromacsEngine2
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
    engine.propagate(path,
                     {'system': system, 'order_function': order_parameter,
                      'interfaces': interfaces},
                     reverse=reverse)
    logger.info('Propagation done!')
    return path


def run_plain_gromacs(engine, system, order_parameter, input_conf,
                      steps=25, exe_dir='forward-plain', reverse=False):
    """Run plain GROMACS.

    Parameters
    ----------
    engine : object like :py:class:`.GromacsEngine`
        This is used to obtain the commands for executing GROMACS.
    system : object like :py:class:`.System`
        The system we are propagation.
    order_parameter : object like :py:class:`.OrderParameter`
        An order parameter to calculate.
    input_conf : string
        The input configuration to use for the GROMACS run.
    steps : integer
        The number of steps to run.
    exe_dir : string
        The path to where we should execute GROMACS.
    reverse : boolean
        Selects the time direction.

    """
    logger.info(f'\nRunning {steps} plain GROMACS steps in "{exe_dir}"')
    make_dirs(exe_dir)
    folder = os.path.abspath(exe_dir)
    clean_dir(folder)
    shutil.copy(input_conf, folder)
    shutil.copy(engine.input_files['input'], folder)
    mdp = os.path.join(folder, 'grompp.mdp')
    engine._modify_input(
        os.path.join(folder, os.path.basename(engine.input_files['input'])),
        mdp,
        {'nsteps': (steps - 1) * engine.subcycles},
        delim='='
    )
    input_file = os.path.join(folder, os.path.basename(input_conf))
    tpr_file = os.path.join(folder, 'topol.tpr')
    grompp = [engine.gmx, 'grompp', '-c', input_file, '-f', mdp,
              '-p', engine.input_files['topology'], '-o', tpr_file]
    logger.info(f'Running grompp in {exe_dir}')
    engine.execute_command(grompp, cwd=exe_dir)
    conf_out = f'run.{engine.ext}'
    exe = engine.mdrun.format('topol.tpr', 'run', conf_out).split()
    logger.info(f"Running \"{' '.join(exe)}\"")
    engine.execute_command(exe, cwd=exe_dir)
    energy_file = os.path.join(folder, 'run.edr')
    energy = engine.get_energies(energy_file)
    trr_file = os.path.join(folder, 'run.trr')
    order = []
    for _, data in read_trr_file(trr_file):
        system.particles.pos = data['x']
        system.particles.vel = data['v']
        if reverse:
            system.particles.vel *= -1
        system.box.update_size(np.diagonal(data['box']))
        order.append(order_parameter.calculate(system))
    order = np.array(order)
    energym = np.zeros((len(order), 2))
    energym[:, 0] = energy['kinetic en.']
    energym[:, 1] = energy['potential']
    return energym, order, trr_file, os.path.join(folder, conf_out)


def main(select=1, plot=False):
    """Perform the test."""
    settings = parse_settings_file('engine-run.rst')
    settings['engine']['input_path'] = os.path.abspath(settings['engine']
                                                       ['input_path'])
    steps = settings['simulation']['steps']
    engine = settings['engine']
    if select == 2:
        gro = GromacsEngine2(engine['gmx'], engine['mdrun'],
                             engine['input_path'], engine['timestep'],
                             engine['subcycles'],
                             maxwarn=0,
                             gmx_format=engine.get('gmx_format', 'g96'),
                             write_vel=True, write_force=True)
    else:
        gro = GromacsEngine(engine['gmx'], engine['mdrun'],
                            engine['input_path'], engine['timestep'],
                            engine['subcycles'],
                            maxwarn=0,
                            gmx_format=engine.get('gmx_format', 'g96'),
                            write_vel=True, write_force=True)
    logger.info(f'Testing for: {gro}')
    logger.info(f'Time step: {gro.timestep}')
    logger.info(f'Subcycles: {gro.subcycles}')
    logger.info(f'GMX format: {gro.ext}')
    # Create dummy variables for the test:
    system = System(units='gromacs',
                    box=create_box(cell=[100, 100, 100]),
                    temperature=200)
    system.particles = ParticlesExt(dim=3)
    initial_conf = gro.input_files['conf']
    system.particles.set_pos((initial_conf, None))
    system.particles.set_vel(False)
    interfaces = [-float('inf'), float('inf'), float('inf')]
    order_parameter = PositionVelocity(1472, dim='z', periodic=True)
    start = time.perf_counter()
    pathf = run_in_steps(gro, system, order_parameter, interfaces,
                         steps=steps,
                         exe_dir=f'{select}-forward-step',
                         reverse=False)
    end = time.perf_counter()
    logger.info(f'Time spent: {end - start}')

    # Set state to last point in trajectory:
    phase_point = pathf.phasepoints[-1]
    system = phase_point.copy()
    start = time.perf_counter()
    pathb = run_in_steps(gro, system, order_parameter, interfaces,
                         steps=steps,
                         exe_dir=f'{select}-backward-step',
                         reverse=True)
    end = time.perf_counter()
    logger.info(f'Time spent: {end - start}')

    # Run plain GROMACS:
    start = time.perf_counter()
    plainf = run_plain_gromacs(gro, system, order_parameter, initial_conf,
                               steps=steps,
                               exe_dir=f'{select}-forward-plain',
                               reverse=False)
    end = time.perf_counter()
    logger.info(f'Time spent: {end - start}')
    last_c = plainf[-1]
    last_r = os.path.join(
        os.path.dirname(last_c),
        f'r_{os.path.basename(last_c)}'
    )
    gro._reverse_velocities(last_c, last_r)
    start = time.perf_counter()
    plainb = run_plain_gromacs(gro, system, order_parameter, last_r,
                               steps=steps,
                               exe_dir=f'{select}-backward-plain',
                               reverse=True)
    end = time.perf_counter()
    logger.info(f'Time spent: {end - start}')

    mse_ok = obtain_mses(pathf, pathb, plainf, plainb)

    if plot:
        logger.info('\nPlotting for comparison')
        plot_path_comparison(pathf, pathb, plainf, plainb)

    if not mse_ok:
        logger.error('\nComparison failed!')
        sys.exit(1)


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
            return False
    return True


def obtain_mses(pathf, pathb, plainf, plainb):
    """Obtain some mean squared errors."""
    orderf = [i.order for i in pathf.phasepoints]
    orderb = [i.order for i in reversed(pathb.phasepoints)]
    mses = [(np.array(orderf), 'step-forward'),
            (np.array(orderb), 'step-back'),
            (plainf[1], 'plain-forward'),
            (plainb[1][::-1, :], 'plain-back')]
    mse_ok = mse_combinations('order parameters', mses, tol=1.0e-9)
    if not mse_ok:
        return mse_ok

    ekinf = [i.particles.ekin for i in pathf.phasepoints]
    ekinb = [i.particles.ekin for i in reversed(pathb.phasepoints)]
    mses = [(np.array(ekinf), 'step-forward'),
            (np.array(ekinb), 'step-back'),
            (plainf[0][:, 0], 'plain-forward'),
            (plainb[0][:, 0][::-1], 'plain-back')]
    mse_ok = mse_combinations('kinetic energy', mses, tol=1.0e-4)
    if not mse_ok:
        return mse_ok

    vpotf = [i.particles.vpot for i in pathf.phasepoints]
    vpotb = [i.particles.vpot for i in reversed(pathb.phasepoints)]
    mses = [(np.array(vpotf), 'step-forward'),
            (np.array(vpotb), 'step-back'),
            (plainf[0][:, 1], 'plain-forward'),
            (plainb[0][:, 1][::-1], 'plain-back')]
    mse_ok = mse_combinations('potential energy', mses, tol=1.0e-4)
    return mse_ok


def plot_path_comparison(pathf, pathb, plainf, plainb):
    """Just plot some properties for the paths."""
    orderf = np.array(
        [i.order for i in pathf.phasepoints]
    )
    orderb = np.array(
        [i.order for i in reversed(pathb.phasepoints)]
    )
    fig1 = plt.figure(figsize=(12, 6))
    ax11 = fig1.add_subplot(121)
    ax12 = fig1.add_subplot(122)
    ax11.plot(orderf[:, 0], lw=2, ls='-', marker='o', label='Forward')
    ax11.plot(orderb[:, 0], lw=2, ls='--', marker='^', label='Backward')
    ax11.plot(plainf[1][:, 0], lw=2, ls=':', marker='s',
              label='Forward-plain')
    ax11.plot(plainb[1][:, 0][::-1], lw=2, ls='--', marker='x',
              label='Backward-plain')
    ax11.legend()
    ax11.set_title('Order param 1')
    ax12.plot(orderf[:, 1], lw=2, ls='-', marker='o')
    ax12.plot(orderb[:, 1], lw=2, ls='--', marker='^')
    ax12.plot(plainf[1][:, 1], lw=2, ls=':', marker='s')
    ax12.plot(plainb[1][:, 1][::-1], lw=2, ls='-.', marker='x')
    ax12.set_title('Order param 2')

    ekinf = [i.particles.ekin for i in pathf.phasepoints]
    ekinb = [i.particles.ekin for i in reversed(pathb.phasepoints)]
    fig2 = plt.figure(figsize=(12, 6))
    ax21 = fig2.add_subplot(131)
    ax22 = fig2.add_subplot(132)
    ax23 = fig2.add_subplot(133)
    ax21.plot(ekinf, lw=2, ls='-', marker='o', label='Forward')
    ax21.plot(ekinb, lw=2, ls='--', marker='^', label='Backward')
    ax21.plot(plainf[0][:, 0], lw=2, ls=':', marker='s',
              label='Forward-plain')
    ax21.plot(plainb[0][:, 0][::-1], lw=2, ls='-.', marker='x',
              label='Backward-plain')
    ax21.set_title('Kinetic energy')
    ax21.legend()

    vpotf = [i.particles.vpot for i in pathf.phasepoints]
    vpotb = [i.particles.vpot for i in reversed(pathb.phasepoints)]
    ax22.plot(vpotf, lw=2, ls='-', marker='o', label='Forward')
    ax22.plot(vpotb, lw=2, ls='--', marker='^', label='Backward')
    ax23.plot(plainf[0][:, 1], lw=2, ls=':', marker='s',
              label='Forward-plain')
    ax23.plot(plainb[0][:, 1][::-1], lw=2, ls='-.', marker='x',
              label='Backward-plain')
    ax22.set_title('Potential energy')
    ax23.set_title('Potential energy')
    ax22.legend()
    ax23.legend()
    fig3 = plt.figure(figsize=(12, 6))
    ax31 = fig3.add_subplot(221)
    ax32 = fig3.add_subplot(222)
    ax33 = fig3.add_subplot(223)
    ax34 = fig3.add_subplot(224)
    ax31.scatter(orderf[:, 0], orderb[:, 0], marker='o', label='Backward',
                 alpha=0.8)
    ax31.scatter(orderf[:, 0], plainf[1][:, 0], marker='s',
                 label='Forward-plain', alpha=0.8)
    ax31.scatter(orderf[:, 0], plainb[1][:, 0][::-1], marker='^',
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
    ax32.scatter(orderf[:, 1], plainf[1][:, 1], marker='s',
                 label='Forward-plain', alpha=0.8)
    ax32.scatter(orderf[:, 1], plainb[1][:, 1][::-1], marker='^',
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
    ax33.scatter(ekinf, plainf[0][:, 0], marker='s',
                 label='Forward-plain', alpha=0.8)
    ax33.scatter(ekinf, plainb[0][:, 0][::-1], marker='^',
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
    ax34.scatter(vpotf, plainf[0][:, 1], marker='s',
                 label='Forward-plain', alpha=0.8)
    ax34.scatter(vpotf, plainb[0][:, 1][::-1], marker='^',
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
    PLOT = len(sys.argv) > 2
    main(select=int(sys.argv[1]), plot=PLOT)
