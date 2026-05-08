#!/usr/bin/env python3
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""pyretisrun - An application for running PyRETIS simulations.

This script is a part of the PyRETIS library and can be used for
running simulations from an input script.

usage: pyretisrun.py [-h] -i INPUT [-V] [-f LOG_FILE] [-l LOG_LEVEL] [-p]

PyRETIS

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Location of PyRETIS input file
  -V, --version         show program's version number and exit
  -f LOG_FILE, --log_file LOG_FILE
                        Specify log file to write
  -l LOG_LEVEL, --log_level LOG_LEVEL
                        Specify log level for log file
  -p, --progress        Display a progress meter instead of text output for
                        the simulation

More information about running PyRETIS can be found at: www.pyretis.org
"""
# pylint: disable=invalid-name
import argparse
import datetime
import logging
import os
import pathlib
import signal
import sys
import traceback
# Other libraries:
import tqdm  # For a progress bar
import colorama  # For coloring text
# PyRETIS library imports:
from pyretis import __version__ as VERSION
from pyretis.info import PROGRAM_NAME, URL, CITE, LOGO
from pyretis.core.pathensemble import generate_ensemble_name
from pyretis.setup import create_simulation
from pyretis.inout.screen import REFERENCE
from pyretis.inout.common import (
    check_python_version,
    create_backup,
)
from pyretis.inout.formats.formatter import (
    get_log_formatter,
    setup_console_logging,
)
from pyretis.inout.settings import (
    parse_settings_file,
    write_settings_file
)


_DATE_FMT = '%d.%m.%Y %H:%M:%S'
logger = setup_console_logging()


def use_tqdm(progress):
    """Return a progress bar if we want one.

    Parameters
    ----------
    progress : boolean
        If True, we should use a progress bar, otherwise not.

    Returns
    -------
    out : object like :py:class:`tqdm.tqdm`
        The progress bar, if requested. Otherwise, just a dummy
        iterator.

    """
    if progress:
        pbar = tqdm.tqdm
    else:
        def empty_tqdm(*args, **kwargs):
            """Return an iterator to replace tqdm."""
            if args:
                return args[0]
            return kwargs.get('iterable', None)
        pbar = empty_tqdm
    return pbar


def hello_world(infile, rundir, logfile):
    """Print out a politically correct greeting for PyRETIS.

    Parameters
    ----------
    infile : string
        String showing the location of the input file.
    rundir : string
        String showing the location we are running in.
    logfile : string
        The output log file

    """
    timestart = datetime.datetime.now().strftime(_DATE_FMT)
    pyversion = sys.version.split()[0]
    logger.banner('\n'.join([LOGO]))

    logger.banner('%s version: %s', PROGRAM_NAME, VERSION)

    logger.banner('Start of execution: %s', timestart)
    logger.banner('Python version: %s', pyversion)

    logger.progress('\nRunning in directory: %s', rundir)
    logger.progress('Input file: %s', infile)
    logger.progress('Log file: %s', logfile)


def bye_bye_world():
    """Print out the goodbye message for PyRETIS."""
    timeend = datetime.datetime.now().strftime(_DATE_FMT)
    msgtxt = f'End of {PROGRAM_NAME} execution: {timeend}'
    logger.progress(msgtxt)
    # display some references:
    references = [f'{PROGRAM_NAME} references:']
    references.append(('-')*len(references[0]))
    for line in CITE.split('\n'):
        if line:
            references.append(line)
    reftxt = '\n'.join(references)
    logger.log(REFERENCE, '\n%s', reftxt)
    urltxt = str(URL)
    logger.log(REFERENCE, urltxt)


def run_md_flux_simulation(sim, sim_settings, progress=False):
    """Run a MD-FLUX simulation.

    Parameters
    ----------
    sim : object like :py:class:`.Simulation`
        This is the simulation to run.
    sim_settings : dict
        The simulation settings.
    progress : boolean, optional
        If True, we will display a progress bar, otherwise, we print
        results to the screen.

    """
    logger.progress('Starting MD-Flux simulation')
    tqd = use_tqdm(progress)
    sim.engine.exe_dir = sim_settings['simulation']['exe_path']
    sim.set_up_output(sim_settings, progress=progress)
    for _ in tqd(sim.run(), initial=sim.cycle['startcycle'],
                 total=sim.cycle['endcycle'], desc='MD-flux'):
        pass
    # Write final restart file:
    sim.write_restart(now=True)
    return True


def run_md_simulation(sim, sim_settings, progress=False):
    """Run a MD simulation.

    Parameters
    ----------
    sim : object like :py:class:`.Simulation`
        This is the simulation to run.
    sim_settings : dict
        The simulation settings.
    progress : boolean, optional
        If True, we will display a progress bar, otherwise, we print
        results to the screen.

    """
    logger.progress('Starting MD simulation')
    tqd = use_tqdm(progress)
    sim.engine.exe_dir = sim_settings['simulation']['exe_path']
    sim.set_up_output(sim_settings, progress=progress)
    for _ in tqd(sim.run(), initial=sim.cycle['startcycle'],
                 total=sim.cycle['endcycle'], desc='MD step'):
        pass
    # Write final restart file:
    sim.write_restart(now=True)
    return True


def explore_simulation(sim, sim_settings, progress=False):
    """Run a RETIS simulation with PyRETIS.

    Parameters
    ----------
    sim : object like :py:class:`.Simulation`
        This is the simulation to run.
    sim_settings : dict
        The simulation settings.
    progress : boolean, optional
        If True, we will display a progress bar, otherwise, we print
        results to the screen.

    """
    sim.set_up_output(
        sim_settings,
        progress=progress
    )

    logtxt = 'Load frames for free energy landscape exploration'
    logger.progress(logtxt)

    # Make sure that the settings are correct. No users don't know better.
    for s_ens in sim_settings.get('ensemble', []):
        s_ens['tis']['freq'] = 0
        s_ens['tis']['allowmaxlength'] = True

    # Here we do the initialisation:
    if not sim.initiate(sim_settings):
        logger.progress('Initiation stopped, will exit now.')
        return False
    sim.write_restart(now=True)

    logtxt = 'Initiation done. Exploring now.'
    logger.progress(logtxt)

    tqd = use_tqdm(progress)

    desc = f'{sim_settings["simulation"]["task"]} Simulation'
    for _ in tqd(sim.run(), initial=sim.cycle['startcycle'],
                 total=sim.cycle['endcycle'], desc=desc):
        pass
    # Write final restart files:
    sim.write_restart(now=True)
    return True


def run_path_simulation(sim, sim_settings, progress=False):
    """Run a RETIS simulation with PyRETIS.

    Parameters
    ----------
    sim : object like :py:class:`.Simulation`
        This is the simulation to run.
    sim_settings : dict
        The simulation settings.
    progress : boolean, optional
        If True, we will display a progress bar, otherwise, we print
        results to the screen.

    """
    sim.set_up_output(
        sim_settings,
        progress=progress
    )

    task = sim_settings['simulation']['task']
    sep = '=' * (len(task) + 26)
    logtxt = f'\n{sep}\n  Initialising {task} simulation.\n{sep}'
    logger.banner(logtxt)

    logger.progress('\nInitialising path ensembles:')

    # Here we do the initialisation:
    if not sim.initiate(sim_settings):
        logger.progress('Initiation stopped, will exit now.')
        return False
    sim.write_restart(now=True)

    logger.progress('\nInitiation done.')

    sep = '=' * (len(task) + 22)
    logtxt = f'\n{sep}\n  Starting {task} simulation.\n{sep}'
    logger.banner(logtxt)

    tqd = use_tqdm(progress)

    desc = f"{sim_settings['simulation']['task']} Simulation. "
    for _ in tqd(sim.run(), initial=sim.cycle['startcycle'],
                 total=sim.cycle['endcycle'], desc=desc):
        pass
    # Write final restart files:
    sim.write_restart(now=True)
    return True


def make_tis_files(_, settings, progress=False):
    """Create TIS simulations input files PyRETIS.

    It just writes out input files for single TIS simulations and
       exit without running a simulation.

    Parameters
    ----------
    settings : list of dicts or Simulation objects
        The settings for the simulations.

    """
    _ = progress
    logtxt = 'Input settings requests: TIS for multiple path ensembles.'
    logger.progress(logtxt)
    logtxt = 'Will create input files for the TIS simulations and exit'
    logger.progress(logtxt)
    i_ens = 0
    for i, ens_settings in enumerate(settings['ensemble']):
        i_ens += 1
        if i == 0 and not settings['simulation']['zero_ensemble']:
            i_ens += 1
        ens_settings['simulation']['zero_ensemble'] = False
        ens_settings['simulation']['task'] = 'tis'
        ensf = generate_ensemble_name(i_ens)
        logtxt = f'Creating input for TIS ensemble: {i_ens} '
        logger.progress(logtxt)
        infile = f'tis-{ensf}.rst'
        logtxt = f'Create file: "{infile}"'
        logger.progress(logtxt)
        exe_dir_file = os.path.join(ens_settings['engine']['exe_path'], infile)
        # Ensure generated TIS input files include a human-friendly
        # orderparameter name without mutating the caller's settings.
        op_section = ens_settings.get('orderparameter', {})
        need_restore = False
        had_name_key = 'name' in op_section
        original_name = op_section.get('name') if had_name_key else None
        if not had_name_key or op_section.get('name') is None:
            ens_settings.setdefault('orderparameter', {})
            ens_settings['orderparameter']['name'] = 'Order Parameter'
            need_restore = True
        try:
            write_settings_file(ens_settings, exe_dir_file, backup=False)
        finally:
            if need_restore:
                if had_name_key:
                    ens_settings['orderparameter']['name'] = original_name
                else:
                    # remove the temporary key we added
                    ens_settings['orderparameter'].pop('name', None)
        logtxt = 'Command for executing:'
        logger.progress(logtxt)
        logtxt = f'pyretisrun -i {infile} -p -f {ensf}.log'
        logger.progress(logtxt)
    return True


def run_generic_simulation(sim, sim_settings, progress=False):
    """Run a generic PyRETIS simulation.

    These are simulations that are just going to complete a given
    number of steps. Other simulation may consist of several
    simulations tied together and these are NOT handled here.

    Parameters
    ----------
    sim : object like :py:class:`.Simulation`
        This is the simulation to run.
    sim_settings : dict
        The simulation settings.
    progress : boolean, optional
        If True, we will display a progress bar, otherwise, we print
        results to the screen.

    """
    logtxt = 'Running simulation'
    logger.progress(logtxt)

    tqd = use_tqdm(progress)
    sim.set_up_output(sim_settings, progress=progress)
    for _ in tqd(sim.run(), desc='Step'):
        pass
    # Write final restart file:
    sim.write_restart(now=True)
    return True


_RUNNERS = {'md-flux': run_md_flux_simulation,
            'md-nve': run_md_simulation,
            'explore': explore_simulation,
            'md': run_md_simulation,
            'make-tis-files': make_tis_files,
            'tis': run_path_simulation,
            'retis': run_path_simulation,
            'pptis': run_path_simulation,
            'repptis': run_path_simulation}


def set_up_simulation(inputfile, runpath):
    """Run all the needed generic set-up.

    Parameters
    ----------
    inputfile : string
        The input file which defines the simulation.
    runpath : string
        The base path we are running the simulation from.

    Returns
    -------
    runner : method
        A method which can be used to execute the simulation.
    sim : object like :py:class:`.Simulation`
        The simulation defined by the input file.
    syst : object like :py:class:`.System`
        The system created.
    sim_settings : dict
        The input settings read from the input file.

    """
    if not os.path.isfile(inputfile):
        raise ValueError(f'Input file "{inputfile}" NOT found!')

    logger.progress('Reading input settings from: %s', inputfile)
    logger.progress('Setting up simulation')

    sim_settings = parse_settings_file(inputfile)
    # NB this is not transmitted to the ensembles
    sim_settings['simulation']['exe_path'] = runpath
    sim_settings['engine']['exe_path'] = runpath
    for ens in sim_settings.get('ensemble', []):
        ens['simulation']['exe_path'] = runpath
        ens['engine']['exe_path'] = runpath

    logtxt = 'Set up and create simulation.'
    logger.info(logtxt)

    sim = create_simulation(sim_settings)

    task = sim_settings['simulation']['task'].lower()
    logger.progress('Setup for simulation "%s" is done.', task)
    runner = _RUNNERS.get(task, run_generic_simulation)
    return runner, sim, sim_settings


def store_simulation_settings(settings, indir, backup):
    """Store the parsed input settings.

    Parameters
    ----------
    settings : dict
        The simulation settings.
    indir : string
        The directory which contains the input script.
    backup : boolean
        If True, an existing settings file will be backed up.

    """
    out_file = os.path.join(indir, 'out.rst')
    logger.info('Full simulation settings written to: %s', out_file)
    write_settings_file(settings, out_file, backup=backup)


def soft_exit_ignore(turn_keyboard_interruption_off=True, exe_dir=None):
    """Manage the KeyboardInterrupt exception.

    Parameters
    ----------
    turn_keyboard_interruption_off : boolean
        If True, instead of regular exiting from the program,
        the file 'EXIT' is created to stop the PyRETIS.
    exe_dir : string, optional
        The path where EXIT file is expected.

    """
    def soft_exit_handler(signum, frame):  # pragma: no cover
        """Handle with a keyboard interruption signal."""
        # pylint: disable=unused-argument
        logger.progress('Attempting soft exit - terminating soon...')
        pathlib.Path(os.path.join(exe_dir, 'EXIT')).touch(exist_ok=True)
    if turn_keyboard_interruption_off:
        return signal.signal(signal.SIGINT, soft_exit_handler)
    return signal.signal(signal.SIGINT, signal.default_int_handler)


def main(infile, indir, exe_dir, progress, log_level):
    """Execute PyRETIS.

    Parameters
    ----------
    infile : string
        The input file to open with settings for PyRETIS.
    indir : string
        The folder containing the settings file.
    exe_dir : string
        The directory we are working from.
    progress : boolean
        Determines if we should use a progress bar or not.
    log_level : integer
        Determines if we should display the error traceback or not.

    """
    simulation = None
    settings = {}
    exit_status = 0

    exit_file = os.path.join(exe_dir, 'EXIT')
    if os.path.isfile(exit_file):
        logger.progress(
            'Exit file found - Remove it before executing PyRETIS.')
        logger.error('*        %s file found         *', exit_file)
        logger.error('Remove the file to execute PyRETIS')
        bye_bye_world()
        return 1

    try:
        run, simulation, settings = set_up_simulation(infile, exe_dir)
        store_simulation_settings(settings, indir, True)
        # Run the simulation:
        soft_exit_ignore(turn_keyboard_interruption_off=True, exe_dir=exe_dir)
        run(simulation, settings, progress=progress)
        soft_exit_ignore(turn_keyboard_interruption_off=False,
                         exe_dir=exe_dir)
    except Exception as error:  # pylint: disable=broad-exception-caught
        exit_status = 1
        logger.error('"%s: %s".', error.__class__.__name__, error.args)
        logger.error('ERROR - execution stopped.')
        logger.error(
            'Please see the LOG for the error message and traceback.'
        )
        # Print the traceback to the log-file, but not to the screen.
        screen = logger.handlers[0]
        lvl = screen.level
        screen.setLevel(logging.CRITICAL + 1)
        logger.error(traceback.format_exc())
        screen.setLevel(lvl)
        if log_level <= logging.DEBUG:
            raise
    finally:
        # Write out the simulation settings as they were parsed and
        # add some additional info:
        if simulation is not None:
            end = getattr(simulation, 'cycle', {'step': None})['step']
            if end is not None:
                settings['simulation']['endcycle'] = end
                logtxt = f'Execution ended at step {end}'
                logger.progress(logtxt)
                logger.success('\nSimulation done\n')
        store_simulation_settings(settings, indir, False)
        bye_bye_world()
    return exit_status


def entry_point():  # pragma: no cover
    """entry_point - The entry point for the pip install of pyretisrun."""
    colorama.init(autoreset=True)
    parser = argparse.ArgumentParser(description=PROGRAM_NAME)
    parser.add_argument('-i', '--input',
                        help=f'Location of {PROGRAM_NAME} input file',
                        required=True)
    parser.add_argument('-V', '--version', action='version',
                        version=f'{PROGRAM_NAME} {VERSION}')
    parser.add_argument('-f', '--log_file',
                        help='Specify log file to write',
                        required=False,
                        default=f'{PROGRAM_NAME.lower()}.log')
    parser.add_argument('-l', '--log_level',
                        help='Specify log level for log file',
                        required=False,
                        default='INFO')
    parser.add_argument('-p', '--progress', action='store_true',
                        help=('Display a progress meter instead of text '
                              'output for the simulation'))
    args_dict = vars(parser.parse_args())

    input_file = args_dict['input']
    # Store directories:
    cwd_dir = os.getcwd()
    input_dir = os.path.dirname(input_file)
    if not os.path.isdir(input_dir):
        input_dir = os.getcwd()

    # Define a file logger:
    create_backup(args_dict['log_file'])
    fileh = logging.FileHandler(args_dict['log_file'], mode='a')
    log_levl = getattr(logging, args_dict['log_level'].upper(),
                       logging.INFO)
    fileh.setLevel(log_levl)
    fileh.setFormatter(get_log_formatter(log_levl))
    logger.addHandler(fileh)
    # Here, we just check the python version. PyRETIS should anyway
    # fail before this for python2.
    check_python_version()

    hello_world(input_file, cwd_dir, args_dict['log_file'])
    sys.exit(main(input_file, input_dir, cwd_dir,
                  args_dict['progress'], log_levl))


if __name__ == '__main__':  # pragma: no cover
    entry_point()
