# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Module for formatting path ensemble data from PyRETIS.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PathEnsembleFormatter (:py:class:`.PathEnsembleFormatter`)
    A class for formatting path ensemble data.

PathEnsembleFile (:py:class:`.PathEnsembleFile`)
    A class for handling PyRETIS path ensemble files.

"""
import logging
from pyretis.core.pathensemble import PathEnsemble
from pyretis.inout.formats.formatter import OutputFormatter
from pyretis.inout.fileio import FileIO
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = ['PathEnsembleFormatter', 'PathEnsembleFile']


class PathEnsembleFormatter(OutputFormatter):
    """A class for formatting path ensemble data.

    This class will effectively define the data which we store
    for each path ensemble. The data is stored in columns with the
    format defined below and contains:

    1) The current cycle number: `'Step'`.

    2) The number of accepted paths: `'No.-acc'`.

    3) The number of shooting moves attempted: `'No.shoot'`.

    4) Starting location (with respect to interfaces): `'l'`.
       This can be `L` if a path starts on the left side, `R` if it
       starts on the right side and `*` if it did not reach the
       interface.

    5) Marker for crossing he middle interface: (`'m'`). This is
       either `M` (when the middle interface is crossed by the path)
       or `*` (if the middle interface is no crossed).

    6) End point for the path: `'r'`. The possible values here are
       similar to the values for the starting location given above.

    7) The length of the generated path: `'Length'`.

    8) The status of the path: `'Acc'`. This is one of the possible
       path statuses defined in :py:mod:`pyretis.core.path`.

    9) The type of move executed for generating the path: `'Mc'`.
       This is one of the moves defined in :py:mod:`pyretis.core.path`.

    10) The smallest order parameter in the path: `'Min-O'`.

    11) The largest order parameter in the path: `'Max-O'`.

    12) The index in the path for the smallest order parameter: `'Idx-Min'`.

    13) The index in the path for the largest order parameter: `'Idx-Max'`.

    14) The order parameter for the shooting point, immediately after
        the shooting move: `'O-shoot'`.

    15) The index in the source path for the shooting point: `'Idx-sh'`.

    16) The index in the new path for the shooting point: `'Idx-shN'`.

    17) The statistical weight of the path: `'Weight'`.

    """

    PATH_FMT = (
        '{0:>10d} {1:>10d} {2:>10d} {3:1s} {4:1s} {5:1s} {6:>7d} '
        '{7:3s} {8:2s} {9:>16.9e} {10:>16.9e} {11:>7d} {12:>7d} '
        '{13:>16.9e} {14:>7d} {15:7d} {16:>16.9e}'
    )
    HEADER = {
        'labels': ['Step', 'No.-acc', 'No.-shoot',
                   'l', 'm', 'r', 'Length', 'Acc', 'Mc',
                   'Min-O', 'Max-O', 'Idx-Min', 'Idx-Max',
                   'O-shoot', 'Idx-sh', 'Idx-shN', 'Weight'],
        'width': [10, 10, 10, 1, 1, 1, 7, 3, 2, 16, 16, 7, 7, 16, 7, 7, 7],
    }

    def __init__(self):
        """Initialise the formatter for path ensemble data."""
        super().__init__('PathEnsembleFormatter', header=self.HEADER)
        self.print_header = True

    def format(self, step, data):
        """Format path ensemble data.

        Here we generate output based on the last path from the
        given path ensemble.

        Parameters
        ----------
        step : int
            The current cycle/step number.
        data : object like :py:class:`.PathEnsemble`
            The path ensemble we will format here.

        Yields
        ------
        out : string
            The formatted line with the data from the last path.

        """
        path_ensemble = data
        path_dict = path_ensemble.paths[-1]
        interfaces = ['*' if i is None else i for i in path_dict['interface']]

        yield self.PATH_FMT.format(
            step,
            path_ensemble.nstats['ACC'],
            path_ensemble.nstats['nshoot'],
            interfaces[0],
            interfaces[1],
            interfaces[2],
            path_dict['length'],
            path_dict['status'],
            path_dict['generated'][0],
            path_dict['ordermin'][0],
            path_dict['ordermax'][0],
            path_dict['ordermin'][1],
            path_dict['ordermax'][1],
            path_dict['generated'][1],
            path_dict['generated'][2],
            path_dict['generated'][3],
            path_dict['weight'],
        )

    @staticmethod
    def parse(line):
        """Parse a line to a simplified representation of a path.

        Parameters
        ----------
        line : string
            The line of text to parse.

        Returns
        -------
        out : dict
            The dict with the simplified path information.

        """
        linec = line.strip()
        if linec.startswith('#'):
            # This is probably the comment
            return None
        # stip trailing comments
        linec = linec.split('#')[0]
        data = [i.strip() for i in linec.split()]
        if len(data) < 16:
            logger.warning(
                'Incorrect number of columns in path data, skipping line.'
            )
            return None
        if len(data) == 16:
            path_info = {'cycle': int(data[0]),
                         'generated': [str(data[8]), float(data[13]),
                                       int(data[14]), int(data[15])],
                         'interface': (str(data[3]), str(data[4]),
                                       str(data[5])),
                         'length': int(data[6]),
                         'ordermax': (float(data[10]), int(data[12])),
                         'ordermin': (float(data[9]), int(data[11])),
                         'status': str(data[7]),
                         'weight': 1.}

        else:
            path_info = {
                'cycle': int(data[0]),
                'generated': [str(data[8]), float(data[13]),
                              int(data[14]), int(data[15])],
                'interface': (str(data[3]), str(data[4]), str(data[5])),
                'length': int(data[6]),
                'ordermax': (float(data[10]), int(data[12])),
                'ordermin': (float(data[9]), int(data[11])),
                'status': str(data[7]),
                'weight': float(data[16]),
            }

        return path_info

    def load(self, filename):
        """Yield the different paths stored in the file.

        The lines are read on-the-fly, converted and yielded one-by-one.

        Parameters
        ----------
        filename : string
            The path/filename to open.

        Yields
        ------
        out : dict
            The information for the current path.

        """
        try:
            with open(filename, 'r', encoding="utf8") as fileh:
                for line in fileh:
                    path_data = self.parse(line)
                    if path_data is not None:
                        path_data['filename'] = filename
                        yield path_data
        except IOError as error:
            logger.warning('I/O error (%d): %s', error.errno, error.strerror)
        except Exception as error:  # pragma: no cover
            logger.warning('Error: %s', error)
            raise


class PathEnsembleFile(PathEnsemble, FileIO):
    """A class for handling PyRETIS path ensemble files.

    This class inherits from the :py:class:`.PathEnsemble` and
    this makes it possible to run the analysis directly on this
    file.

    """

    def __init__(self, filename, file_mode, ensemble_settings=None,
                 backup=True):
        """Set up the file-like object.

        Parameters
        ----------
        filename : string
            The file to open.
        file_mode : string
            Determines the mode for opening the file.
        ensemble_settings : dict, optional
            Ensemble specific settings.
        backup : boolean, optional
            Determines how we handle existing files when the mode
            is set to writing.

        """
        default_settings = {
            'ensemble_number': 0,
            'interfaces': [0.0, 1.0, 2.0],
            'detect': None
        }
        settings = {}
        for key, val in default_settings.items():
            try:
                settings[key] = ensemble_settings[key]
            except (TypeError, KeyError):
                settings[key] = val
                logger.warning(
                    'No "%s" ensemble setting given for "%s". Using defaults',
                    key,
                    self.__class__
                )
        PathEnsemble.__init__(self, settings['ensemble_number'],
                              settings['interfaces'])
        FileIO.__init__(self, filename, file_mode, PathEnsembleFormatter(),
                        backup=backup)

    def get_paths(self):
        """Yield paths from the file."""
        return self.load()
