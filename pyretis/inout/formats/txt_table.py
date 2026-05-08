# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Module defining a table-like output format.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TxtTableFormatter (:py:class:`.TxtTableFormatter`)
    A class for creating a table-like format.

PathTableFormatter (:py:class:`.PathTableFormatter`)
    A class for table-like output from path simulations.

ThermoTableFormatter (:py:class:`.ThermoTableFormatter`)
    A class for thermodynamic (energy) output. Useful for output from
    MD-simulations.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

txt_save_columns (:py:class:`.txt_save_columns`)
    For writing simple column-based output.

"""
import logging
import numpy as np
from pyretis.inout.common import create_backup
from pyretis.inout.formats.formatter import OutputFormatter
from pyretis.core.path import _GENERATED_SHORT
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())

np.set_printoptions(legacy='1.21')

__all__ = [
    'TxtTableFormatter',
    'PathTableFormatter',
    'ThermoTableFormatter',
    'RETISResultFormatter',
    'txt_save_columns',
]


def txt_save_columns(outputfile, header, variables, backup=False):
    """Save variables to a text file using ``numpy.savetxt``.

    Note that the variables are assumed to be numpy.arrays of equal
    shape and that the output file may also be a compressed file in
    gzip format (this is selected by letting the output file name
    end with '.gz').

    Parameters
    ----------
    outputfile : string
        Name of the output file to create.
    header : string
        A string that will be written at the beginning of the file.
    variables : tuple or list of numpy.arrays
        These are the variables that will be saved to the text file.
    backup : boolean, optional
        Determines if we should backup old files or not.

    """
    if backup:
        msg = create_backup(outputfile)
        if msg:
            logger.warning(msg)
    nvar = len(variables)
    mat = np.zeros((len(variables[0]), nvar))
    for i, vari in enumerate(variables):
        try:
            mat[:, i] = vari
        except ValueError:
            msg = 'Could not align variables, skipping (writing zeros)'
            logger.warning(msg)
    np.savetxt(outputfile, mat, header=header)


def _fill_list(the_list, length, fillvalue=None):
    """Fill a list to a specified length.

    Parameters
    ----------
    the_list : list
        The list to fill.
    length : int
        The required length.
    fillvalue : optional
        The value to insert. If None is given the last item in the list
        is repeated.

    """
    if fillvalue is None:
        fillvalue = the_list[-1]
    while len(the_list) < length:
        the_list.append(fillvalue)


class TxtTableFormatter(OutputFormatter):
    """A class for generating table output.

    This class handles formatting of output data to a table-like
    format.

    Attributes
    ----------
    variables : list of strings
        These are the variables we will use in the table.
    fmt : string
        The format to apply to the columns.
    row_fmt : list of strings
        A list of strings used for formatting, used to construct `fmt`.
    title : string
        A title for the table.

    Example
    -------
    For creating a new table, a dictionary is convenient for
    grouping the settings:

    >>> tabl = {
    ...    'title': 'Energy output',
    ...    'var': ['step', 'temp', 'vpot'
    ...            'ekin', 'etot', 'press'],
    ...    'format': {'labels': ['Step', 'Temp', 'Pot',
    ...               'Kin', 'Tot', 'Press'],
    ...               'width': (10, 12),
    ...               'spacing': 2,
    ...               'row_fmt': ['{:> 10d}', '{:> 12.6g}']
    ...              }
    ...    }
    >>> table = TxtTableFormatter(tabl['var'], tabl['title'], **tabl['format'])


    """

    def __init__(self, variables, title, **kwargs):
        """Initialise the TxtTable object.

        Parameters
        ----------
        variables : list of strings
            These are the variables we will use in the table. If the
            header is not specified, then we will create one using
            these variables.
        title : string
            A title for the table.
        kwargs : dict
            Additional settings for the formatter. This may contain:

            * width : list of ints
                The (maximum) width of the columns. If the number of
                items in this list is smaller than the number of
                variables, we simply repeat the last width for the
                number of times we need.
            * labels : list of strings
                Table headers to use for the columns.
            * spacing : integer
                The separation between columns. The default value is 1.
            * row_fmt : list of strings
                The format to apply to the columns. If the number of
                items in this list is smaller than the number of
                variables, we simply repeat the last width for the
                number of times we need.

        """
        spacing = kwargs.get('spacing', 1)
        header = {'spacing': spacing,
                  'labels': kwargs.get('labels', list(variables))}
        width = kwargs.get('width', None)
        if width is None:
            header['width'] = [max(12, len(i)) for i in header['labels']]
        else:
            header['width'] = list(width)

        _fill_list(header['width'], len(header['labels']))

        super().__init__('TxtTableFormatter', header=header)
        self.title = title
        self.variables = variables
        row_fmt = kwargs.get('row_fmt', None)
        if row_fmt is None:
            self.row_fmt = []
            for wid in header['width']:
                if wid - 6 <= 0:
                    self.row_fmt.append(f'{{:> {wid}}}')
                else:
                    self.row_fmt.append(f'{{:> {wid}.{wid-6}g}}')
        else:
            self.row_fmt = row_fmt
        _fill_list(self.row_fmt, len(self.variables))
        str_white = ' ' * spacing
        self.fmt = str_white.join(self.row_fmt)

    def format(self, step, data):
        """Generate output from a dictionary using the requested variables.

        Parameters
        ----------
        step : int
            This is the current step number or a cycle number in a
            simulation.
        data : dict
            This is assumed to a dictionary containing the row items we
            want to format.

        Yields
        ------
        out : string
            A line with the formatted output.

        """
        var = []
        for i in self.variables:
            if i == 'step':
                var.append(step)
            else:
                var.append(data.get(i, float('nan')))
        txt = self.fmt.format(*var)
        yield txt

    def __str__(self):
        """Return a string with some info about the TxtTableFormatter."""
        msg = [f'TxtTableFormatter: "{self.title}"']
        msg += [f'\t* Variables: {self.variables}']
        msg += [f'\t* Format: {self.fmt}']
        return '\n'.join(msg)


class PathTableFormatter(TxtTableFormatter):
    """A special table output class for path ensembles.

    This object will return a table of text with a header and with
    formatted rows for a path ensemble. The table rows will contain
    data from the `PathEnsemble.nstats` variable. This table is just
    meant as output to the screen during a path ensemble simulation.
    """

    def __init__(self):
        """Initialise the PathTableFormatter."""
        title = 'Path Ensemble Statistics'
        var = ['step', 'ACC', 'BWI',
               'BTL', 'FTL', 'BTX', 'FTX']
        table_format = {'labels': ['Cycle', 'Accepted', 'BWI', 'BTL', 'FTL',
                                   'BTX', 'FTX'],
                        'width': (10, 12),
                        'spacing': 2,
                        'row_fmt': ['{:> 10d}', '{:> 12d}']}
        super().__init__(var, title, **table_format)

    def format(self, step, data):
        """Generate the output for the path table.

        Here we overload the :py:meth:`.TxtTableFormatter.format` method
        in order to write path ensemble statistics to (presumably)
        the screen.

        Parameters
        ----------
        step : int
            This is the current step number or a cycle number in a
            simulation.
        data : object like :py:class:`.PathEnsemble`
            This is the path ensemble we are generating output for.

        Yield
        -----
        out : string
            The formatted output.

        """
        row = {}
        for key in self.variables:
            if key == 'step':
                value = step
            else:
                value = data.nstats.get(key, 0)
            row[key] = value
        var = [row.get(i, float('nan')) for i in self.variables]
        yield self.fmt.format(*var)


class ThermoTableFormatter(TxtTableFormatter):
    """A special text table for energy output.

    This object will return a table of text with a header and with
    formatted rows for energy output. Typical use is in MD simulation
    where we want to display energies at different steps in the
    simulations.
    """

    def __init__(self):
        """Initialise the ThermoTableFormatter."""
        title = 'Energy Output'
        var = ['step', 'temp', 'vpot', 'ekin', 'etot', 'press']
        table_format = {'labels': ['Step', 'Temp', 'Pot',
                                   'Kin', 'Tot', 'Press'],
                        'width': (10, 12),
                        'spacing': 2,
                        'row_fmt': ['{:> 10d}', '{:> 12.6g}']}
        super().__init__(var, title, **table_format)


class RETISResultFormatter(TxtTableFormatter):
    """A special table output class for path ensembles in RETIS simulations.

    This object will return a table of text with a header and with
    formatted rows for a path ensemble. The table rows will contain
    data from the `PathEnsemble.nstats` variable. This table is just
    meant as output to the screen during a path ensemble simulation.
    """

    def __init__(self):
        """Initialise the PathTableFormatter."""
        title = 'Path Ensemble Statistics'
        var = ['pathensemble', 'step', 'ACC', 'BWI',
               'BTL', 'FTL', 'BTX', 'FTX']
        table_format = {
            'labels': [
                'Ensemble', 'Cycle', 'Accepted', 'BWI', 'BTL', 'FTL',
                'BTX', 'FTX'
            ],
            'width': (8, 8, 10),
            'spacing': 2,
            'row_fmt': ['{:>10}', '{:> 8d}', '{:> 10d}']
        }
        super().__init__(var, title, **table_format)
        self.print_header = False

    def format(self, step, data):
        """Generate the output for the path table.

        Here we overload the :py:meth:`.TxtTableFormatter.format` method
        in order to write path ensemble statistics to (presumably)
        the screen.

        Parameters
        ----------
        step : int
            This is the current step number or a cycle number in a
            simulation.
        data : object like :py:class:`.PathEnsemble`
            This is the path ensemble we are generating output for.

        Yield
        -----
        out : string
            The formatted output.

        """
        row = {}
        for key in self.variables:
            if key == 'step':
                value = step
            elif key == 'pathensemble':
                value = data.ensemble_name
            else:
                value = data.nstats.get(key, 0)
            row[key] = value
        yield (f'# Results for path ensemble {data.ensemble_name} '
               f'at cycle {step}:')
        path = data.paths[-1]
        move = _GENERATED_SHORT.get(path['generated'][0], 'unknown').lower()
        yield (f'# Generated path with status "{path["status"]}", '
               f'move "{move}" and length {path["length"]}.')
        omax = path['ordermax']
        yield f'# Order parameter max was: {omax[0]} at index {omax[1]}.'
        omin = path['ordermin']
        yield f'# Order parameter min was: {omin[0]} at index {omin[1]}.'
        yield '# Path ensemble statistics:'
        yield self.header
        var = [row.get(i, float('nan')) for i in self.variables]
        yield self.fmt.format(*var)
        yield '\n'
