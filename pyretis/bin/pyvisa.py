#!/usr/bin/env python3
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""pyvisa - An application for analysing PyRETIS simulations.

This script is a part of the PyRETIS library and can be used for
analysing the result from simulations.

Usage::

    pyvisa.py [-h]

Optional arguments::

    -cmp --pyvisa_compressor  use only the compressor tool without GUI.
    -data --pyvisa-data       select the data source.
    -h, --help                show this help message and exit.
    -i INPUT, --input INPUT   location of PyRETIS input files
                              or PyVisA compressed file.
    -oo --only-order          read only the order.txt files.
    -V, --version             show program's version number and exit.

"""
# pylint: disable=invalid-name
import argparse
import datetime
import os
import sys
import colorama  # pylint: disable=import-error
from pyretis import __version__ as VERSION
from pyretis.bin.pyretisanalyse import write_traceback
from pyretis.pyvisa.info import PROGRAM_NAME, CITE, LOGO, URL
from pyretis.inout.common import (
    check_python_version)
from pyretis.inout.formats.formatter import setup_console_logging
from pyretis.pyvisa.orderparam_density import PathDensity, pyvisa_compress
from pyretis.pyvisa.common import recalculate_all
from pyretis.__init__ import HAS_PYVISA_REQ
from pyretis.inout.screen import REFERENCE  # registers custom levels
if HAS_PYVISA_REQ:
    from pyretis.pyvisa.visualize import visualize_main  # pragma: no cover


# Set up for logging:
logger = setup_console_logging()

runpath = os.getcwd()
ERROR_FILE = 'pyvisa_error.txt'


def hello_pyvisa(run_dir, infile):
    """Output a standard greeting for PyVISA.

    Parameters
    ----------
    run_dir : string
        The location where we are executing the analysis.
    infile : string
        String showing the location of the input file.

    """
    pyversion = sys.version.split()[0]
    msgtxt = [LOGO]
    msgtxt += ['                                                    Starting']
    msgtxt += [f'{PROGRAM_NAME} version: {VERSION}']
    msgtxt += [f'Python version: {pyversion}']
    msgtxt += ['', f'Running in directory: {run_dir}']
    msgtxt += [f'Input file: {infile}']
    logger.banner('\n'.join(msgtxt))


def bye_pyvisa():
    """Print out the goodbye message for PyVisA."""
    _DATE_FMT = '%d.%m.%Y %H:%M:%S'
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


def pyvisa_visual(basepath, input_file, pyvisa_dict):  # pragma: no cover
    """Load data to PyVisA.

    Parameters
    ----------
    basepath : string
        The execution folder where the input files are.
    input_file : string
        The input file with settings for the analysis.
    pyvisa_dict : dictionary, optional
        It determines the section of pyvisa to use, it contains:

    """
    # pylint: disable=possibly-used-before-assignment
    if not HAS_PYVISA_REQ:
        msg = ('Requirements not installed. You can still generate the '
               'hdf5 by using the -cmp flag instead')
        raise ImportError(msg)

    if input_file is None:
        visualize_main(basepath)
    # If a compressed file is given as an input
    elif input_file.endswith(('.hdf5', '.zip')):
        visualize_main(basepath, input_file)

    elif input_file.endswith('rst'):
        # If -data is selected with a trajectory
        if pyvisa_dict['pyvisa_data'] is not None:
            visualize_main(basepath, input_file,
                           trajfile=pyvisa_dict['pyvisa_data'])
        else:
            p_data = PathDensity(basepath, iofile=input_file)
            p_data.infos['only_ops'] = pyvisa_dict['only_order']

            # Walk selected ensemble directories
            p_data.walk_dirs(only_ops=pyvisa_dict['only_order'])
            visualize_main(basepath, p_data, rst_file=input_file)

    else:
        msg = f'The format of {input_file} is not supported.'
        raise ImportError(msg)


def main(basepath, input_file, pyvisa_dict=None):
    """Run the analysis.

    Parameters
    ----------
    basepath : string
        The execution folder where the input files are.
    input_file : string
        The input file with settings for the analysis.
    pyvisa_dict : dictionary, optional
        It determines the section of pyvisa to use, it contains:

         * `pyvisa_cmp`, boolean
           If true, only the compressor tool will be executed. A compressed
           file will be produced.
         * `pyvisa_visual`, boolean
           If true, only the visualization tool will be executed.
         * `pyvisa_data`, str
           If given, the file or folder containing the files that
           will be used to feed to PyVisA.
         * `pyvisa_recalculate`, boolean
           If true, use the recalculation tool to compute new op
           and cv values.

    """
    exit_status = 0
    try:
        if input_file is not None and \
                not os.path.isfile(os.path.join(basepath, input_file)):
            raise FileNotFoundError(f'Input file "{input_file}" NOT found!')
        # For recalculation of op and cv's
        if pyvisa_dict.get('pyvisa_recalculate', False):
            recalculate_all(basepath, input_file,
                            data=pyvisa_dict['pyvisa_data'])
        # For compression of the files ONLY
        elif pyvisa_dict.get('pyvisa_compressor', False):
            pyvisa_compress(basepath, input_file, pyvisa_dict)
        else:
            pyvisa_visual(basepath, input_file, pyvisa_dict)

    # pylint: disable=broad-exception-caught
    except Exception as error:  # pragma: no cover
        exit_status = 1
        errtxt = f'{type(error).__name__}: {error.args}'
        logger.error(errtxt)
        logger.error('Execution failed!')
        logger.error('Error traceback is written to: %s', ERROR_FILE)
        write_traceback(os.path.join(basepath, ERROR_FILE))
    finally:
        bye_pyvisa()
    return exit_status


def entry_point():  # pragma: no cover
    """entry_point - The entry point for the pip install of pyretisanalyse."""
    colorama.init(autoreset=True)
    parser = argparse.ArgumentParser(description=PROGRAM_NAME)
    parser.add_argument('-i', '--input',
                        help=(f'Location of {PROGRAM_NAME} input file'),
                        required=False, default=None)
    parser.add_argument('-V', '--version', action='version',
                        version=f'{PROGRAM_NAME} {VERSION}')
    parser.add_argument('-cmp', '--pyvisa_compressor',
                        action='store_true',
                        help='Run PyVisA compressor tool',
                        default=False)
    parser.add_argument('-data', '--pyvisa-data',
                        action='store',
                        help='Select PyVisA data',
                        default=None)
    parser.add_argument('-recalculate', '--pyvisa-recalculate',
                        action='store_true',
                        help='Recalculate op and cv data.',
                        default=False)
    parser.add_argument('-oo', '--only_order', action='store_true',
                        help=('PyVisA: Use only data from order.txt files'),
                        default=False)

    args_dict = vars(parser.parse_args())

    check_python_version()
    infile = args_dict['input']
    basepath = os.getcwd() if infile is None else os.path.dirname(
        os.path.abspath(infile))
    hello_pyvisa(basepath, infile)
    sys.exit(main(basepath, infile, args_dict))


if __name__ == '__main__':  # pragma: no cover
    entry_point()
