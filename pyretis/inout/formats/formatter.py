# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Module for defining the generic formatting of output data from PyRETIS.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OutputFormatter (:py:class:`.OutputFormatter`)
    A generic class for formatting output from PyRETIS.

PyretisLogFormatter (:py:class:`.PyretisLogFormatter`)
    A class representing a formatter for the PyRETIS log file.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

apply_format (:py:func:`.apply_format`)
    Apply a format string to a given float value. This method
    can be used for formatting text for tables (i.e. if we want
    a fixed width).

format_number (:py:func:`.format_number`)
    Format a number based on its size.

get_log_formatter (:py:func:`.get_log_formatter`)
    Select a formatter for logging based on a given message level.

"""
import logging
import sys
import colorama
from pyretis.inout.fileio import read_some_lines
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


# hard-coded formats to use for Log files:
LOG_FMT = '[%(levelname)s]: %(message)s'
LOG_DEBUG_FMT = ('[%(levelname)s] [%(name)s, %(funcName)s() at'
                 ' line %(lineno)d]: %(message)s')


__all__ = ['OutputFormatter', 'PyretisLogFormatter',
           'ColoredLogFormatter', 'apply_format',
           'format_number', 'get_log_formatter', 'setup_console_logging',
           'LOG_FMT', 'LOG_DEBUG_FMT']


def _make_header(labels, width, spacing=1):
    """Format a table header with the given labels.

    Parameters
    ----------
    labels : list of strings
        The strings to use for the table header.
    width : list of ints
        The widths to use for the table.
    spacing : int
        The spacing between the columns in the table

    Returns
    -------
    out : string
        A header for the table.

    """
    heading = []
    for i, col in enumerate(labels):
        try:
            wid = width[i]
        except IndexError:
            wid = width[-1]
        if i == 0:
            fmt = f'# {{:>{wid-2}s}}'
        else:
            fmt = f'{{:>{wid}s}}'
        heading.append(fmt.format(col))
    str_white = ' ' * spacing
    return str_white.join(heading)


def apply_format(value, fmt):
    """Apply a format string to a given float value.

    Here we check the formatting of a float. We are *forcing* a
    *maximum length* on the resulting string. This is to avoid problems
    like: '{:7.2f}'.format(12345.7) which returns '12345.70' with a
    length 8 > 7. The intended use of this function is to avoid such
    problems when we are formatting numbers for tables. Here it is done
    by switching to an exponential notation. But note however that this
    will have implications for how many decimal places we can show.

    Parameters
    ----------
    value : float
        The float to format.
    fmt : string
        The format to use.

    Note
    ----
    This function converts numbers to have a fixed length. In some
    cases this may reduce the number of significant digits. Remember
    to also output your numbers without this format in case a specific
    number of significant digits is important!

    """
    maxlen = fmt.split(':')[1].split('.')[0]
    align = ''
    if not maxlen[0].isalnum():
        align = maxlen[0]
        maxlen = maxlen[1:]
    maxlen = int(maxlen)
    str_fmt = fmt.format(value)
    if len(str_fmt) > maxlen:  # switch to exponential:
        if value < 0:
            deci = maxlen - 7
        else:
            deci = maxlen - 6
        new_fmt = f'{{:{align}{maxlen}.{deci}e}}'
        return new_fmt.format(value)
    return str_fmt


def format_number(number, minf, maxf, fmtf='{0:<16.9f}', fmte='{0:<16.9e}'):
    """Format a number based on its size.

    Parameters
    ----------
    number : float
        The number to format.
    minf : float
        If the number is smaller than `minf` then apply the
        format with scientific notation.
    maxf : float
        If the number is greater than `maxf` then apply the
        format with scientific notation.
    fmtf : string, optional
        Format to use for floats.
    fmte : string, optional
        Format to use for scientific notation.

    Returns
    -------
    out : string
        The formatted number.

    """
    if minf <= number <= maxf:
        return fmtf.format(number)
    return fmte.format(number)


class OutputFormatter:
    """A generic class for formatting output from PyRETIS.

    Attributes
    ----------
    name : string
        A string which identifies the formatter.
    header : string
        A header (or table heading) with information about the
        output data.
    print_header : boolean
        Determines if we are to print the header or not. This is useful
        for classes making use of the formatter to determine if they
        should write out the header or not.

    """

    _FMT = '{}'

    def __init__(self, name, header=None):
        """Initialise the formatter.

        Parameters
        ----------
        name : string
            A string which identifies the output type of this formatter.
        header : dict, optional
            The header for the output data

        """
        self.name = name
        self._header = None
        self.print_header = True
        if header is not None:
            if 'width' in header and 'labels' in header:
                self._header = _make_header(header['labels'],
                                            header['width'],
                                            spacing=header.get('spacing', 1))
            else:
                self._header = header.get('text', None)
        else:
            self.print_header = False

    @property
    def header(self):
        """Define the header as a property."""
        return self._header

    @header.setter
    def header(self, value):
        """Set the header."""
        self._header = value

    def format(self, step, data):
        """Use the formatter to generate output.

        Parameters
        ----------
        step : integer
            This is assumed to be the current step number for
            generating the output.
        data : list, dict or similar
            This is the data we are to format. Here we assume that
            this is something we can iterate over.

        """
        out = [f'{step}']
        for i in data:
            out.append(self._FMT.format(i))
        yield ' '.join(out)

    @staticmethod
    def parse(line):
        """Parse formatted data.

        This method is intended to be the "inverse" of the :py:meth:`.format`
        method. In this particular case, we assume that we read floats from
        columns in a file. One input line corresponds to a "row" of data.

        Parameters
        ----------
        line : string
            The string we will parse.

        Returns
        -------
        out : list of floats
            The parsed input data.

        """
        return [int(col) if i == 0 else
                float(col) for i, col in enumerate(line.split())]

    def load(self, filename):
        """Read generic data from a file.

        Since this class defines how the data is formatted it is also
        convenient to have methods for reading the data defined here.
        This method will read entire blocks of data from a file into
        memory. This will be slow for large files and this method
        could be converted to also yield the individual "rows" of
        the blocks, rather than the full blocks themselves.

        Parameters
        ----------
        filename : string
            The path/file name of the file we want to open.

        Yields
        ------
        data : list of tuples of int
            This is the data contained in the file. The columns are the
            step number, interface number and direction.

        See Also
        --------
        :py:func:`.read_some_lines`.

        """
        for blocks in read_some_lines(filename, self.parse):
            data_dict = {'comment': blocks['comment'],
                         'data': blocks['data']}
            yield data_dict

    def __str__(self):
        """Return basic info about the formatter."""
        return self.name


def get_log_formatter(level):
    """Select a log format based on a given level.

    Here, it is just used to get a slightly more verbose format for
    the debug level.

    Parameters
    ----------
    level : integer
        This integer defines the log level.

    Returns
    -------
    out : object like :py:class:`logging.Formatter`
        An object that can be used as a formatter for a logger.

    """
    if level <= logging.DEBUG:
        return PyretisLogFormatter(LOG_DEBUG_FMT)
    return PyretisLogFormatter(LOG_FMT)


class PyretisLogFormatter(logging.Formatter):  # pragma: no cover
    """Hard-coded formatter for the PyRETIS log file.

    This formatter will just adjust multi-line messages to have some
    indentation.
    """

    def format(self, record):
        """Apply the PyRETIS log format."""
        out = logging.Formatter.format(self, record)
        if '\n' in out:
            heading, _ = out.split(record.message)
            if len(heading) < 12:
                out = out.replace('\n', '\n' + ' ' * len(heading))
            else:
                out = out.replace('\n', '\n' + ' ' * 4)
        return out


class ColoredLogFormatter(logging.Formatter):  # pragma: no cover
    """Colored formatter for console output.

    Maps logging levels to ANSI colors via colorama for colored
    console output via the standard ``logger``.
    """

    COLORS = {
        logging.DEBUG: colorama.Fore.CYAN,
        logging.INFO: colorama.Fore.BLUE,
        25: colorama.Fore.GREEN,   # PROGRESS level
        26: colorama.Fore.CYAN,    # BANNER level (logo, version)
        27: colorama.Fore.WHITE,   # REFERENCE level (citations, URLs)
        logging.WARNING: colorama.Fore.YELLOW,
        logging.ERROR: colorama.Fore.RED,
        logging.CRITICAL: colorama.Fore.RED + colorama.Style.BRIGHT,
    }

    # Levels that show only the message (no [LEVEL]: prefix) on console
    _CLEAN_LEVELS = {25, 26, 27}  # PROGRESS, BANNER, REFERENCE

    def format(self, record):
        """Apply color based on the log level."""
        color = self.COLORS.get(record.levelno, '')
        if record.levelno in self._CLEAN_LEVELS:
            # Show just the message, no prefix, for progress/banner/reference
            out = record.getMessage()
        else:
            out = logging.Formatter.format(self, record)
        if color:
            return color + out + colorama.Style.RESET_ALL
        return out


def setup_console_logging(level: int = None):
    """Set up the root logger with a colored console handler.

    This configures the root logger at DEBUG level (so all messages
    are passed to handlers) and attaches a colored StreamHandler to
    stdout that displays only PROGRESS (25) and above on the console.
    INFO and DEBUG messages are suppressed on the console and are
    intended for file handlers added separately.

    Parameters
    ----------
    level : int, optional
        The minimum log level for the console handler.
        Defaults to PROGRESS (25) if not specified.

    Returns
    -------
    root_logger : logging.Logger
        The configured root logger.

    """
    # Import here to avoid circular imports at module level
    from pyretis.inout.screen import PROGRESS  # noqa: PLC0415
    if level is None:
        level = PROGRESS
    root_logger = logging.getLogger('')
    root_logger.setLevel(logging.DEBUG)
    # Guard against duplicate handlers when multiple modules call this:
    for handler in root_logger.handlers:
        if (isinstance(handler, logging.StreamHandler) and
                getattr(handler, 'stream', None) is sys.stdout):
            return root_logger
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(level)
    console.setFormatter(ColoredLogFormatter(LOG_FMT))
    root_logger.addHandler(console)
    return root_logger
