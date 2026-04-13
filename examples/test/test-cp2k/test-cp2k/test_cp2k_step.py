# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test a single step with the CP2K external engine."""
import os
import time
import colorama
from pyretis.core import System, create_box, ParticlesExt
from pyretis.inout.common import make_dirs
from pyretis.inout.settings import parse_settings_file
from pyretis.engines.cp2k import (
    CP2KEngine,
)
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def clean_dir(dirname):
    """Remove ALL files in the given directory!"""
    for files in os.listdir(dirname):
        filename = os.path.join(dirname, files)
        if os.path.isfile(filename):
            os.remove(filename)


def run_step(engine, system, exe_dir='forward-single-step'):
    """Run the engine forward in time, in steps.

    Parameters
    ----------
    engine : object like :py:class:`.ExternalMDEngine
        Engine to use for propagation.
    system : object like :py:class:`.System`
        The system we are propagation.
    exe_dir : string
        The folder to use for the execution.
    """
    logger.info(f'\nRunning a single CP2K step in "{exe_dir}"')
    make_dirs(exe_dir)
    folder = os.path.abspath(exe_dir)
    clean_dir(folder)
    engine.exe_dir = folder
    engine.step(system, 'single')
    logger.info('Propagation done!')


def test_genvel(engine, input_file, exe_dir='genvel'):
    """Test generation of velocities.

    Parameters
    ----------
    engine : object like :py:class:`.ExternalMDEngine
        Engine to use for the generation.
    input_file : string
        The input configuration used for CP2K.
    exe_dir : string, optional
        The directory where we will be running CP2K.
    """
    logger.info(f'\nRunning CP2K genvel step in "{exe_dir}"')
    make_dirs(exe_dir)
    folder = os.path.abspath(exe_dir)
    clean_dir(folder)
    engine.exe_dir = folder
    engine._prepare_shooting_point(input_file)
    logger.info('Propagation done!')


def main():
    """Execute the test."""
    settings = parse_settings_file('engine.rst')
    engine_settings = settings['engine']
    engine = CP2KEngine(
        engine_settings['cp2k'],
        engine_settings['input_path'],
        engine_settings['timestep'],
        engine_settings['subcycles'],
        engine_settings.get('extra_files', [])
    )
    logger.info(f'Testing engine: {engine}')
    logger.info(f'Time step: {engine.timestep}')
    logger.info(f'Subcycles: {engine.subcycles}')
    system = System(units='gromacs',
                    box=create_box(cell=[100, 100, 100]),
                    temperature=500)
    system.particles = ParticlesExt(dim=3)
    initial_conf = engine.input_files['conf']
    system.particles.set_pos((initial_conf, None))
    system.particles.set_vel(False)

    start = time.perf_counter()
    run_step(engine, system, exe_dir='forward-single-step')
    end = time.perf_counter()
    logger.info(f'Time spent: {end - start}')

    start = time.perf_counter()
    test_genvel(engine, engine.input_files['conf'], exe_dir='genvel')
    end = time.perf_counter()
    logger.info(f'Time spent: {end - start}')


if __name__ == '__main__':
    colorama.init(autoreset=True)
    main()
