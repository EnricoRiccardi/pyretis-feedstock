# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This file contains common functions for the input/output.

It contains some slave functions that are used in the in/output function
of PyRETIS.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OutputBase (:py:class:`.OutputBase`)
    A base class for handling the output.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

check_python_version (:py:func:`.check_python_version`)
    A method that will give warnings when we use older and possibly
    unsupported Python versions.

create_backup (:py:func:`.create_backup`)
    A function to handle the creation of backups of old files.

make_dirs (:py:func:`.make_dirs`)
    Create directories (for path simulation).

create_empty_ensembles (:py:func:`.create_ensembles`)
    A method to prepare the ensembles inputs in settings

simplify_ensemble_name (:py:func:`.simplify_ensemble_name`)
    Simplify the name of ensembles for creating directories.

generate_file_name (:py:func:`.generate_file_name`)
    Generate file name for an output task, from settings.

"""
import logging
import errno
import os
import re
import sys
from abc import ABCMeta, abstractmethod
from pyretis import MIN_PYTHON_VERSION
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = [
    'check_python_version',
    'create_backup',
    'create_empty_ensembles',
    'make_dirs',
    'OutputBase',
    'simplify_ensemble_name',
    'add_dirname',
    'name_file',
    'generate_file_name',
    'TRJ_FORMATS'
]


# Hard-coded patters for energy analysis output files.
# These are just used to make it simpler to change these default
# names in the future.
ENERFILES = {'energies': 'energies',
             'run_energies': 'runenergies',
             'temperature': 'temperature',
             'run_temp': 'runtemperature',
             'block': '{}block',
             'dist': '{}dist'}
# hard-coded information for the energy terms:
ENERTITLE = {'vpot': 'Potential energy',
             'ekin': 'Kinetic energy',
             'etot': 'Total energy',
             'ham': 'Hamilt. energy',
             'temp': 'Temperature',
             'elec': 'Energy (externally computed)'}
# hard-coded patters for flux analysis output files:
FLUXFILES = {'runflux': 'runflux_{}',
             'block': 'errflux_{}'}
# order files:
ORDERFILES = {'order': 'order',
              'ordervel': 'orderv',
              'run_order': 'runorder',
              'dist': 'orderdist',
              'block': 'ordererror',
              'msd': 'ordermsd'}
# hard-coded patters for path analysis output files:
PATHFILES = {'pcross': '{}_pcross',
             'prun': '{}_prun',
             'prun_sl': '{}_prun_sl',
             'prun_sr': '{}_prun_sr',
             'perror': '{}_perror',
             'perror_sl': '{}_perror_sl',
             'perror_sr': '{}_perror_sr',
             'lengtherror': '{}_lerror',
             'pathlength': '{}_lpath',
             'shoots': '{}_shoots'}
# hard-coded patterns for matched files:
PATH_MATCH = {'total': 'total-probability',
              'progress': 'overall-rrun',
              'error': 'overall-err',
              'match': 'matched-probability'}
# Supported trj formats # TODO (LAMMPS missing)
TRJ_FORMATS = ['.xyz', '.gro', '.trr', '.txt', '.g96']


def create_backup(outputfile):
    """Check if a file exist and create backup if requested.

    This function will check if the given file name exists and if it
    does, it will move that file to a new file name such that the given
    one can be used without overwriting.

    Parameters
    ----------
    outputfile : string
        This is the name of the file we wish to create.

    Returns
    -------
    out : string
        This string is None if no backup is made, otherwise, it will
        just say what file was moved (and to where).

    Note
    ----
    No warning is issued here. This is just in case the `msg` returned
    here will be part of some more elaborate message.

    """
    filename = f'{outputfile}'
    fileid = 0
    msg = None
    while os.path.isfile(filename) or os.path.isdir(filename):
        filename = f'{outputfile}_{fileid:03d}'
        fileid += 1
    if fileid > 0:
        msg = f'Backup existing file "{outputfile}" to "{filename}"'
        os.rename(outputfile, filename)
    return msg


def make_dirs(dirname):
    """Create directories for path simulations.

    This function will create a folder using a specified path.
    If the path already exists and if it's a directory, we will do
    nothing. If the path exists and is a file we will raise an
    `OSError` exception here.

    Parameters
    ----------
    dirname : string
        This is the directory to create.

    Returns
    -------
    out : string
        A string with some info on what this function did. Intended for
        output.

    """
    try:
        os.makedirs(dirname)
        msg = f'Created directory: "{dirname}"'
    except OSError as err:
        if err.errno != errno.EEXIST:  # pragma: no cover
            raise err
        if os.path.isfile(dirname):
            msg = f'"{dirname}" is a file. Will abort!'
            raise OSError(errno.EEXIST, msg) from err
        if os.path.isdir(dirname):
            msg = f'Directory "{dirname}" already exist.'
    return msg


def simplify_ensemble_name(ensemble, fmt='{:03d}'):
    """Simplify path names for file/directory names.

    Here, we are basically translating ensemble names to more friendly
    names for directories and files that is:

    - ``[0^-]`` returns ``000``,
    - ``[0^+]`` returns ``001``,
    - ``[1^+]`` returns ``002``, etc.

    Parameters
    ----------
    ensemble : string
        This is the string to simplify.
    fmt : string. optional
        This is a format to use for the directories.

    """
    match_ensemble = re.search(r'(?<=\[)(\d+)(?=\^)', ensemble)
    if match_ensemble:
        ens = int(match_ensemble.group())
    else:
        match_ensemble = re.search(r'(?<=\[)(\d+)(?=\])', ensemble)
        if match_ensemble:
            ens = int(match_ensemble.group())
        else:
            return ensemble  # Assume that the ensemble is OK as it is.
    match_dir = re.search(r'(?<=\^)(.)(?=\])', ensemble)
    if match_dir:
        dire = match_dir.group()
        if dire != '-':
            ens += 1
    else:
        msg = [f'Could not get direction for ensemble {ensemble}.',
               'We assume "+", note that this might overwrite files']
        logger.warning('\n'.join(msg))
        ens += 1
    return fmt.format(ens)


def add_dirname(filename, dirname):
    """Add a directory as a prefix to a filename, i.e. `dirname/filename`.

    Parameters
    ----------
    filename : string
        The filename.
    dirname : string
        The directory we want to prefix. It can be None, in which
        case we ignore it.

    Returns
    -------
    out : string
        The path to the resulting file.

    """
    if dirname is not None:
        return os.path.join(dirname, filename)
    return filename


def name_file(name, extension, path=None):
    """Return a file name by joining a name and an file extension.

    This function is used to create file names. It will use
    `os.extsep` to create the file names and `os.path.join` to add a
    path name if the `path` is given. The returned file name will be of
    form (example for posix): ``path/name.extension``.

    Parameters
    ----------
    name : string
        This is the name, without extension, for the file.
    extension : string
        The extension to use for the file name.
    path : string, optional
        An optional path to add to the file name.

    Returns
    -------
    out : string
        The resulting file name.

    """
    return add_dirname(os.extsep.join([name, extension]), path)


def generate_file_name(basename, directory, settings):
    """Generate file name for an output task, from settings.

    Parameters
    ----------
    basename : string
        The base file name to use.
    directory : string
        A directory to output to. Can be None to output to the
        current working directory.
    settings : dict
        The input settings

    Returns
    -------
    filename : string
        The file name to use.

    """
    prefix = settings['output'].get('prefix', None)
    if prefix is not None:
        filename = f'{prefix}{basename}'
    else:
        filename = basename
    filename = add_dirname(filename, directory)
    return filename


def check_python_version():
    """Give a warning about old python version(s)."""
    version = sys.version_info[:2]  # Only check major, minor
    if version < MIN_PYTHON_VERSION:
        msgtxt = ('Please upgrade to Python {}.{}.'
                  '\nPython {}.{} is not supported!')
        msgtxt = msgtxt.format(*MIN_PYTHON_VERSION,
                               *version)
        logger.error(msgtxt)
        raise SystemExit(msgtxt)


class OutputBase(metaclass=ABCMeta):
    """A generic class for handling output.

    Attributes
    ----------
    formatter : object like py:class:`.OutputFormatter`
        The object responsible for formatting output.
    target : string
        Determines where the target for the output, for
        instance "screen" or "file".
    first_write : boolean
        Determines if we have written something yet, or
        if this is the first write.

    """

    target = None

    def __init__(self, formatter):
        """Create the object and attach a formatter."""
        self.formatter = formatter
        self.first_write = True

    def output(self, step, data):
        """Use the formatter to write data to the file.

        Parameters
        ----------
        step : int
            The current step number.
        data : list
            The data we are going to output.

        """
        if self.first_write and self.formatter.print_header:
            self.first_write = False
            self.write(self.formatter.header)
        for line in self.formatter.format(step, data):
            self.write(line)

    @abstractmethod
    def write(self, towrite, end='\n'):
        """Write a string to the output defined by this class.

        Parameters
        ----------
        towrite : string
            The string to write.
        end : string, optional
            A "terminator" for the given string.

        Returns
        -------
        status : boolean
            True if we managed to write, False otherwise.

        """
        return

    def formatter_info(self):
        """Return a string with info about the formatter."""
        if self.formatter is not None:
            return self.formatter.__class__
        return None

    def __str__(self):
        """Return basic info."""
        return f'{self.__class__.__name__}\n\t* Formatter: {self.formatter}'


def create_empty_ensembles(settings):
    """Create missing ensembles in the settings.

    Checks the input and allocate it to the right ensemble. In theory
    inouts shall include all these info, but it is not practical.

    Parameters
    ----------
    settings : dict
        The current input settings.

    Returns
    -------
    None, but this method might add data to the input settings.

    """
    ints = settings['simulation']['interfaces']

    # Determine how many ensembles are needed.
    add0 = 0
    if not settings['simulation'].get('flux', False):
        add0 += 1
    if not settings['simulation'].get('zero_ensemble', True):
        add0 += 1

    # if some ensembles have inputs, they need to be kept.
    if 'ensemble' in settings:
        orig_set = settings['ensemble'].copy()
    else:
        orig_set = []

    settings['ensemble'] = []
    add = add0

    # if in the main settings an ensemble_number is defined, then only
    # that ensemble will be considered.
    if 'tis' in settings:
        idx = settings['tis'].get('ensemble_number')
        detect = settings['tis'].get('detect', ints[2])
        if idx is not None:
            settings['ensemble'].append({'interface': ints[1],
                                         'tis': {'ensemble_number': idx,
                                                 'detect': detect}})
            for sav in orig_set:
                settings['ensemble'][0] = {**settings['ensemble'][0], **sav}
            return
    # if one wants to compute the flux, the 000 ensemble is for it
    # todo remove this labelling mismatch, and give to the flux
    # a flux name folder (instead of 000), and leave 000 for the O^+ ens.
    if settings['simulation'].get('flux', False):
        settings['ensemble'].append({'interface': ints[0],
                                     'tis': {'ensemble_number': 0,
                                             'detect': ints[1]}})
        add = add0 + 1
    for i in range(add, len(ints)):
        settings['ensemble'].append({'interface': ints[i - 1],
                                     'tis': {'ensemble_number': i,
                                             'detect': ints[i]}})

    # create the ensembles in setting, keeping eventual inputs.
    # nb. in the settings, specific input for an ensemble can be now given.
    for i_ens, ens in enumerate(settings['ensemble']):
        for sav in orig_set:
            if 'tis' in sav and 'ensemble_number' in sav['tis']:
                if ens['tis']['ensemble_number'] ==\
                        sav['tis']['ensemble_number']:
                    settings['ensemble'][i_ens].update(sav)
            elif ens['interface'] == sav['interface']:
                settings['ensemble'][i_ens].update(sav)

    return
