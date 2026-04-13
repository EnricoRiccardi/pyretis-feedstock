#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""pyretisanalyse - An application for analysing PyRETIS simulations.

This script is a part of the PyRETIS library and can be used for
analysing the result from simulations.

usage: pyretisanalyse.py [-h] -i INPUT [-V]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Location of PyRETIS input file
  -V, --version         show program's version number and exit
"""
# pylint: disable=invalid-name
import argparse
import logging
import os
import sys
import traceback
import colorama
from pyretis import __version__ as VERSION
from pyretis.info import PROGRAM_NAME, URL, CITE, LOGO
from pyretis.inout.common import create_backup
from pyretis.core.units import CONSTANTS
from pyretis.core.pathensemble import generate_ensemble_name
from pyretis.inout.analysisio.analysisio import run_analysis
from pyretis.inout.common import (
    check_python_version,
    make_dirs,
    name_file,
)
from pyretis.inout.formats.formatter import setup_console_logging
from pyretis.inout.report import generate_report
from pyretis.inout.settings import parse_settings_file
from pyretis.inout.screen import PROGRESS, REFERENCE  # registers custom levels

# Set up for logging:
logger = setup_console_logging()

runpath = os.getcwd()

# Hard-coded patters for report outputs:
REPORTFILES = {
    'md-flux': 'md_flux_report',
    'retis': 'retis_report',
    'make-tis-files': 'tis-multiple_report',
    'tis': 'tis_report',
    'repptis': 'repptis_report',
}

ERROR_FILE = 'error.txt'


def hello_world(infile, run_dir, report_dir):
    """Output a standard greeting for PyRETIS analysis.

    Parameters
    ----------
    infile : string
        String showing the location of the input file.
    run_dir : string
        The location where we are executing the analysis.
    report_dir : string
        String showing the location of where we write the output.

    """
    pyversion = sys.version.split()[0]
    msgtxt = [LOGO]
    msgtxt += ['                                                    Starting']
    msgtxt += ['analysis tool!']
    msgtxt += [f'{PROGRAM_NAME} version: {VERSION}']
    msgtxt += [f'Python version: {pyversion}']
    msgtxt += ['', f'Running in directory: {run_dir}']
    msgtxt += [f'Report directory: {report_dir}']
    msgtxt += [f'Input file: {infile}']
    logger.banner('\n'.join(msgtxt))


def bye_bye_world():
    """Print out the goodbye message for PyRETIS."""
    msgtxt = f'End of {PROGRAM_NAME} analysis execution.'
    logger.progress(msgtxt)
    # display some references:
    references = [f'{PROGRAM_NAME} references:']
    references.append(('-') * len(references[0]))
    for line in CITE.split('\n'):
        if line:
            references.append(line)
    reftxt = '\n'.join(references)
    logger.log(REFERENCE, '\n' + reftxt)
    urltxt = f'{URL}'
    logger.log(REFERENCE, urltxt)


def write_traceback(filename):
    """Write the error traceback to the given file."""
    msg = create_backup(filename)
    if msg:
        logger.warning(msg)
    with open(filename, 'w', encoding='utf-8') as out:
        out.write(traceback.format_exc())


def get_report_name(report_type, ext, prefix=None, path=None):
    """Generate file name for a report.

    Parameters
    ----------
    report_type : string
        Identifier for the report we are writing.
    ext : string
        Extension for the file to write.
    prefix : string, optional
        A prefix to add to the file name. Usually just used
        to mark reports with ensemble number for `report_type` equal
        to 'tis-single'
    path : string
        A directory to use for saving the report to.

    Returns
    -------
    out : string
        The name of the file written.

    """
    name = REPORTFILES[report_type]
    if prefix is not None:
        name = f'{prefix}_{name}'
    return name_file(name, ext, path=path)


def write_file(outname, report_txt):
    """Write a generated report to a given file.

    Parameters
    ----------
    outname : string
        The name of the file to write/create.
    report_txt : string
        This is the generated report as a string.

    Returns
    -------
    out : string
        The name of the file written.

    """
    with open(outname, 'wt', encoding='utf-8') as report_fh:
        report_fh.write(report_txt)
    return outname


def create_reports(settings, analysis_results, report_path):
    """Create some reports to display the output.

    Parameters
    ----------
    settings : dict
        Settings for analysis (and the simulation).
    analysis_results : dict
        Results from the analysis.
    report_path : string
        The path to the directory where the reports should be saved.

    Yields
    ------
    out : string
        The report files created.

    """
    ens_list = settings.get('ensemble', [])
    if settings['simulation']['task'] == 'tis' and len(ens_list) == 1:
        task = 'tis'
        ens_n = settings['ensemble'][0]['tis'].get('ensemble_number', '1')
        pfix = generate_ensemble_name(ens_n)
    else:
        task = settings['simulation']['task']
        pfix = None
    for report_type in settings['analysis']['report']:
        report, ext = generate_report(task, analysis_results,
                                      output=report_type)
        if report is not None:
            reportfile = get_report_name(task, ext, prefix=pfix,
                                         path=report_path)
            write_file(reportfile, report)
            yield reportfile


def main(input_file, run_path, report_dir):
    """Run the analysis.

    Parameters
    ----------
    input_file : string
        The input file with settings for the analysis.
    run_path : string
        The location from which we are running the analysis.
    report_dir : string
        The location where we will write the report.

    """
    try:
        if input_file is None:
            raise FileNotFoundError('Input file required (-i filename).')
        if not os.path.isfile(os.path.join(run_path, input_file)):
            raise FileNotFoundError(f'Input file "{input_file}" NOT found!')
        # Run analysis
        logger.progress('Reading input file "%s"', input_file)
        settings = parse_settings_file(input_file)
        # override exe-path to the one we are executing in now:
        settings['simulation']['exe-path'] = run_path
        units = settings['system']['units']
        # set derived properties:
        settings['system']['beta'] = (settings['system']['temperature'] *
                                      CONSTANTS['kB'][units]) ** -1
        settings['analysis']['report-dir'] = report_dir
        msg_dir = make_dirs(report_dir)
        logger.info(msg_dir)
        task = settings['simulation']['task']
        logger.progress('Simulation task was: "%s"', task)

        results = run_analysis(settings)
        for outfile in create_reports(settings, results, report_dir):
            relfile = os.path.relpath(outfile, start=run_path)
            logger.progress('Report created: %s', relfile)

    except Exception as error:  # Exceptions should subclass BaseException.
        errtxt = f'{type(error).__name__}: {error.args}'
        logger.error(errtxt)
        logger.error('Execution failed!')
        logger.error('Error traceback is written to: %s', ERROR_FILE)
        write_traceback(os.path.join(run_path, ERROR_FILE))
    finally:
        bye_bye_world()


def entry_point():  # pragma: no cover
    """entry_point - The entry point for the pip install of pyretisanalyse."""
    colorama.init(autoreset=True)
    parser = argparse.ArgumentParser(description=PROGRAM_NAME)
    parser.add_argument(
        '-i',
        '--input',
        help=(f'Location of {PROGRAM_NAME} input file'),
        required=False,
        default=None
    )
    parser.add_argument('-V', '--version', action='version',
                        version=f'{PROGRAM_NAME} {VERSION}')

    args_dict = vars(parser.parse_args())

    check_python_version()

    inputfile = args_dict['input']
    reportdir = os.path.join(runpath, 'report')

    hello_world(inputfile, runpath, reportdir)
    main(inputfile, runpath, reportdir)


if __name__ == '__main__':  # pragma: no cover
    entry_point()
