# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Module for formatting snapshot data from PyRETIS.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SnapshotFormatter (:py:class:`.SnapshotFormatter`)
    Generic class for formatting system snapshots (coordinates). Note
    it is intended to format objects like :py:class:`.System`. The
    format is very similar to the XYZ format.

SnapshotFile (:py:class:`.SnapshotFile`)
    A class for a collection of snapshots.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

adjust_coordinate (:py:func:`.adjust_coordinate`)
    Helper method to add dimensions when formatting data in 1D or 2D
    to a format that requires 3D data.

read_txt_snapshots (:py:func:`.read_txt_snapshots`)
    For reading PyRETIS snapshots from a file.

"""
import logging
import numpy as np
from pyretis.core.common import initiate_instance
from pyretis.inout.formats.formatter import OutputFormatter
from pyretis.inout.fileio import FileIO
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())

np.set_printoptions(legacy='1.21')

__all__ = [
    'adjust_coordinate',
    'read_txt_snapshots',
    'SnapshotFormatter',
    'SnapshotFile',
]


class SnapshotFormatter(OutputFormatter):
    """Generic class for formatting system snapshots.

    Attributes
    ----------
    write_vel : boolean
        If True, we will also format velocities
    fmt : string
        Format to use for position output.
    fmt_vel : string
        Format to use for position and velocity output.

    """

    data_keys = ('atomname', 'x', 'y', 'z', 'vx', 'vy', 'vz')
    _FMT_FULL = '{} {} {} {}'
    _FMT_FULL_VEL = '{} {} {} {} {} {} {}'
    _FMT = '{:5s} {:15.9f} {:15.9f} {:15.9f}'
    _FMT_VEL = '{:5s} {:15.9f} {:15.9f} {:15.9f} {:15.9f} {:15.9f} {:15.9f}'

    def __init__(self, write_vel=True, fmt=None):
        """Initialise the formatter.

        Parameters
        ----------
        write_vel : boolean, optional
            If True, the formatter will attempt to output velocities.
        fmt : string, optional
            Selects the format to use.

        """
        super().__init__('SnapshotFormatter', header=None)
        self.print_header = False
        self.write_vel = write_vel
        if fmt == 'full':
            self.fmt = self._FMT_FULL
            self.fmt_vel = self._FMT_FULL_VEL
        else:
            self.fmt = self._FMT
            self.fmt_vel = self._FMT_VEL

    def format(self, step, data):
        """Generate the snapshot output.

        Parameters
        ----------
        step : integer
            The current step number for generating the output.
        data : object like :py:class:`.System`
            The system we are generating output for.

        """
        yield from self.format_snapshot(step, data)

    def _format_without_vel(self, particles):
        """Format positions of particles for output.

        Parameters
        ----------
        particles : object like :py:class:`.Particles`
            The particles for which we will format information.

        Yields
        ------
        out : string
            The formatted output, to be written.

        """
        pos = adjust_coordinate(particles.pos)
        for namei, posi in zip(particles.name, pos):
            yield self.fmt.format(namei, *posi)

    def _format_with_vel(self, particles):
        """Format positions of particles for output.

        Parameters
        ----------
        particles : object like :py:class:`.Particles`
            The particles for which we will format information.

        Yields
        ------
        out : string
            The formatted output, to be written.

        """
        pos = adjust_coordinate(particles.pos)
        vel = adjust_coordinate(particles.vel)
        for namei, posi, veli in zip(particles.name, pos, vel):
            yield self.fmt_vel.format(namei, posi[0], posi[1], posi[2],
                                      veli[0], veli[1], veli[2])

    @staticmethod
    def parse(line):
        """Not implemented - line parser for snapshots.

        For snapshots, we use a specialized reader that will read entire
        system snapshots. To avoid confusion, we just give a warning here.

        """
        logger.warning('The line parser is not implemented for the'
                       'snapshot reader.')

    def format_snapshot(self, step, system):
        """Format the given snapshot.

        Parameters
        ----------
        step : int
            The current simulation step.
        system : object like :py:class:`.System`
            The system object with the positions to format.

        Returns
        -------
        out : list of strings
            The formatted snapshot.

        """
        npart = system.particles.npart
        buff = [
            f'{npart}',
            f'Snapshot, step: {step} box: {system.box.print_length()}',
        ]
        if self.write_vel:
            formatter = self._format_with_vel
        else:
            formatter = self._format_without_vel
        for lines in formatter(system.particles):
            buff += [lines]
        return buff

    def load(self, filename):
        """Read snapshots from a given file.

        Parameters
        ----------
        filename : string
            The path/filename to open.

        Yields
        ------
        out : dict
            This dict contains the snapshot.

        """
        yield from read_txt_snapshots(filename, data_keys=self.data_keys)


def read_txt_snapshots(filename, data_keys=None):
    """Read snapshots from a text file.

    Parameters
    ----------
    filename : string
        The file to read from.
    data_keys : tuple of strings, optional
        This tuple determines the data we are to read. It can
        be of type ``('atomname', 'x', 'y', 'z', ...)``.

    Yields
    ------
    out : dict
        A dictionary with the snapshot.

    """
    lines_to_read = 0
    snapshot = None
    if data_keys is None:
        data_keys = ('atomname', 'x', 'y', 'z', 'vx', 'vy', 'vz')
    read_header = False
    with open(filename, 'r', encoding="utf8") as fileh:
        for lines in fileh:
            if read_header:
                snapshot = {'header': lines.strip()}
                snapshot['box'] = get_box_from_header(snapshot['header'])
                read_header = False
                continue
            if lines_to_read == 0:  # new snapshot
                if snapshot is not None:
                    yield snapshot
                try:
                    lines_to_read = int(lines.strip())
                except ValueError:
                    logger.error('Error in the input file %s', filename)
                    raise
                read_header = True
                snapshot = None
            else:
                lines_to_read -= 1
                data = lines.strip().split()
                for i, (val, key) in enumerate(zip(data, data_keys)):
                    if i == 0:
                        value = val.strip()
                    else:
                        value = float(val)
                    try:
                        snapshot[key].append(value)
                    except KeyError:
                        snapshot[key] = [value]
    if snapshot is not None:
        yield snapshot


def get_box_from_header(header):
    """Get box lengths from a text header.

    Parameters
    ----------
    header : string
        Header from which we will extract the box.

    Returns
    -------
    out : numpy.array or None
        The box lengths.

    """
    low = header.lower()
    if low.find('box:') != -1:
        txt = low.split('box:')[1].strip()
        return np.array([float(i) for i in txt.split()])
    return None


def adjust_coordinate(coord):
    """Adjust the dimensionality of coordinates.

    This is a helper method for trajectory writers.

    A lot of the different formats expects us to have 3 dimensional
    data. This method just adds dummy dimensions equal to zero.

    Parameters
    ----------
    coord : numpy.array
        The real coordinates.

    Returns
    -------
    out : numpy.array
        The "zero-padded" coordinates.

    """
    if len(coord.shape) == 1:
        npart, dim = len(coord), 1
    else:
        npart, dim = coord.shape
    if dim == 3:
        # correct dimensionality, just stop here:
        return coord
    adjusted = np.zeros((npart, 3))
    try:
        for i in range(dim):
            adjusted[:, i] = coord[:, i]
    except IndexError:
        if dim == 1:
            adjusted[:, 0] = coord
    return adjusted


class SnapshotFile(FileIO):
    """A class for a collection of snapshots."""

    def __init__(self, filename, file_mode, backup=True, format_settings=None):
        """Create the snapshot file with the possible settings."""
        if format_settings is not None:
            fmt = initiate_instance(SnapshotFormatter, format_settings)
        else:
            fmt = SnapshotFormatter()
        super().__init__(filename, file_mode, fmt, backup=backup)
