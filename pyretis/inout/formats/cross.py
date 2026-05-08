# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Module for formatting crossing data from PyRETIS.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CrossFormatter (:py:class:`.CrossFormatter`)
    A class for formatting crossing data from flux simulations.

CrossFile (:py:class:`.CrossFile`)
    A class for handling PyRETIS crossing files.

"""
import logging
from pyretis.inout.formats.formatter import OutputFormatter
from pyretis.inout.fileio import FileIO
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = [
    'CrossFormatter',
    'CrossFile',
]


class CrossFormatter(OutputFormatter):
    """A class for formatting crossing data from flux simulations.

    This class handles formatting of crossing data. The format for
    the crossing file is three columns:

    1) The first column is the step number (an integer).

    2) The second column is the interface number (an integer). These are
       numbered from 1 (_NOT_ from 0).

    3) The direction we are moving in - `+` for the positive direction
       or `-` for the negative direction. Internally this is converted
       to an integer (`+1` or `-1`).

    """

    # Format for crossing files:
    CROSS_FMT = '{:>10d} {:>4d} {:>3s}'

    def __init__(self):
        """Initialise the crossing formatter."""
        header = {'labels': ['Step', 'Int', 'Dir'], 'width': [10, 4, 3]}
        super().__init__('CrossFormatter', header=header)

    def format(self, step, data):
        """Generate output data to be written to a file or screen.

        This will format the crossing data in a space delimited form.


        Parameters
        ----------
        step : int
            This is the current step number. It is only used here for
            debugging and can possibly be removed. However, it's useful
            to have here since this gives a common interface for all
            formatters.
        data : list of tuples
            The tuples are crossing with interfaces (if any) on the form
            `(timestep, interface, direction)` where the direction
            is '-' or '+'.

        Yields
        ------
        out : string
            The line(s) to be written.

        See Also
        --------
        :py:meth:`.check_crossing` in :py:mod:`pyretis.core.path`
        calculates the tuple `cross` which is used in this routine.

        Note
        ----
        We add 1 to the interface number here. This is for
        compatibility with the old FORTRAN code where the interfaces
        are numbered 1, 2, ... rather than 0, 1, ... .

        """
        msgtxt = f'Generating crossing data at step: {step}'
        logger.debug(msgtxt)
        for cro in data:
            if cro:
                yield self.CROSS_FMT.format(cro[0], cro[1] + 1, cro[2])

    @staticmethod
    def parse(line):
        """Parse crossing data.

        The format is described in :py:meth:`.format`, this method will
        parse this format for a _single_ row of data.

        Parameters
        ----------
        line : string
            A line to parse.

        Returns
        -------
        out[0] : integer
            The step number.
        out[1] : integer
            The interface number.
        out[2] : integer
            The direction, left (-1) or right (1).

        Note
        ----
        A '1' will be subtracted from the interface in analysis. This is
        just for backward compatibility with the old FORTRAN code.

        """
        linessplit = line.strip().split()
        step, inter = int(linessplit[0]), int(linessplit[1])
        direction = -1 if linessplit[2] == '-' else 1
        return step, inter, direction


class CrossFile(FileIO):
    """A class for handling PyRETIS crossing files."""

    def __init__(self, filename, file_mode, backup=True):
        """Create the cross file with correct formatter."""
        super().__init__(filename, file_mode, CrossFormatter(), backup=backup)
