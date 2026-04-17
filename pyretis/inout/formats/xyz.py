# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Module for handling the output/input of trajectory data.

This module defines some classes for writing out and reading snapshots
and trajectories in a XYZ-like format. By XYZ-like we here mean
that we also include velocities as part of the file.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

format_xyz_data (:py:func:`.format_xyz_data`)
    A method for formatting position/velocity data into a
    XYZ-like format. This can be used by external engines to
    convert to a standard format.

read_xyz_file (:py:func:`.read_xyz_file`)
    A method for reading snapshots from a XYZ file.

reverse_xyz_file (:py:func:`.reverse_xyz_file`)
    Method to read an XYZ-file in reverse, helps :py:func:`.xyz_merge`
    merge trajectories.

txt_to_xyz (:py:func:`.txt_to_xyz`)
    Method for converting/extracting paths from the internal trajectory
    format and writing them as a .xyz file.

write_xyz_file (:py:func:`.write_xyz_file`)
    Just a convenience method for writing to a new file.

write_xyz_trajectory (:py:func:`.write_xyz_trajectory`)
    A helper method to write xyz trajectories. This will
    also attempt to write box information to the XYZ-header.

xyz_merge (:py:func:`.xyz_merge`)
    A method to merge a forward and a backward XYZ trajectory. This
    method is mainly used when creating a whole trajectory for
    visualisation.

"""
import io
import logging
import numpy as np
from pyretis.inout.formats.snapshot import (
    adjust_coordinate,
    read_txt_snapshots,
)
from pyretis.inout.formats.path import PathIntFormatter
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())

np.set_printoptions(legacy='1.21')

# define formats for the trajectory output:
_XYZ_FMT = '{0:5s} {1:8.3f} {2:8.3f} {3:8.3f}'
_XYZ_BIG_FMT = '{:5s}' + 3*' {:15.9f}'
_XYZ_BIG_VEL_FMT = _XYZ_BIG_FMT + 3*' {:15.9f}'


__all__ = [
    'format_xyz_data',
    'read_xyz_file',
    'reverse_xyz_file',
    'write_xyz_file',
    'write_xyz_trajectory',
    'xyz_merge',
    'convert_snapshot'
]


def read_xyz_file(filename):
    """Read files in XYZ format.

    This method will read a XYZ file and yield the different snapshots
    found in the file.

    Parameters
    ----------
    filename : string
        The file to open.

    Yields
    ------
    out : dict
        This dict contains the snapshot.

    Examples
    --------
    >>> from pyretis.inout.formats.xyz import read_xyz_file
    >>> for snapshot in read_xyz_file('traj.xyz'):
    ...     print(snapshot['x'][0])

    Note
    ----
    The positions will **NOT** be converted to a specified set of units.

    """
    xyz_keys = ('atomname', 'x', 'y', 'z', 'vx', 'vy', 'vz')
    yield from read_txt_snapshots(filename, data_keys=xyz_keys)


def convert_snapshot(snapshot):
    """Convert a XYZ snapshot to numpy arrays.

    Parameters
    ----------
    snapshot : dict
        The dict containing a snapshot read from a XYZ-file.

    Returns
    -------
    box : numpy.array, 1D
        The box dimensions if we manage to read it.
    xyz : numpy.array
        The positions.
    vel : numpy.array
        The velocities.
    names : list of strings
        The atom names found in the file.

    """
    names = snapshot['atomname']
    box = snapshot.get('box', None)
    natom = len(names)
    xyz = np.zeros((natom, 3))
    vel = np.zeros_like(xyz)
    for i, dim in enumerate(('x', 'y', 'z')):
        xyz[:, i] = snapshot[dim]
        key = f'v{dim}'
        if key in snapshot:
            vel[:, i] = snapshot[key]
    return box, xyz, vel, names


def format_xyz_data(pos, vel=None, names=None, header=None, fmt=None):
    """Format XYZ data for outputting.

    Parameters
    ----------
    pos : numpy.array
       The positions to write.
    vel : numpy.array, optional
       The velocities to write.
    names : list, optional
        The atom names.
    header : string, optional
        Header to use for writing the XYZ-file.
    fmt : string, optional
        A format to use for the writing

    Yields
    ------
    out : string
        The formatted lines.

    """
    npart = len(pos)
    pos = adjust_coordinate(pos)

    if fmt is None:
        fmt = _XYZ_BIG_FMT if vel is None else _XYZ_BIG_VEL_FMT

    if vel is not None:
        vel = adjust_coordinate(vel)
    yield str(npart)

    if header is None:
        yield 'PyRETIS XYZ writer'
    else:
        yield str(header)

    if names is None:
        logger.warning('No atom name given. Using "X"')

    if vel is None:
        for i in range(npart):
            namei = 'X' if names is None else names[i]
            yield fmt.format(namei, pos[i, 0], pos[i, 1], pos[i, 2])
    else:
        for i in range(npart):
            namei = 'X' if names is None else names[i]
            yield fmt.format(namei, pos[i, 0], pos[i, 1], pos[i, 2],
                             vel[i, 0], vel[i, 1], vel[i, 2])


def write_xyz_file(filename, pos, vel=None, names=None, header=None):
    """Create a new XYZ file with the given file name.

    Parameters
    ----------
    filename : string
        The file to create.
    pos : numpy.array
       The positions to write.
    vel : numpy.array, optional
       The velocities to write.
    names : list, optional
        The atom names.
    header : string, optional
        Header to use for writing the XYZ-file.

    Note
    ----
    We will here just overwrite if the file already exists.

    Examples
    --------
    >>> import numpy as np
    >>> from pyretis.inout.formats.xyz import write_xyz_file
    >>> xyz = 10 * np.random.rand(10, 3)
    >>> write_xyz_file('conf.xyz', xyz)
    >>> vel = 10 * np.random.rand(10, 3)
    >>> write_xyz_file('confv.xyz', xyz, vel)

    """
    with open(filename, 'w', encoding='utf-8') as output_file:
        for line in format_xyz_data(pos, vel=vel, names=names,
                                    header=header):
            output_file.write(f'{line}\n')


def write_xyz_trajectory(filename, pos, vel, names, box, step=None,
                         append=True):
    """Write XYZ snapshot to a trajectory.

    This is intended as a lightweight alternative for just
    dumping snapshots to a trajectory file.

    Parameters
    ----------
    filename : string
        The file name to dump to.
    pos : numpy.array
        The positions we are to write.
    vel : numpy.array
        The velocities we are to write.
    names : list of strings
        Atom names to write.
    box : numpy.array
        The box dimensions/vectors
    step : integer, optional
        If the ``step`` is given, then the step number is
        written to the header.
    append : boolean, optional
        Determines if we append or overwrite an existing file.

    Note
    ----
    We will here append to the file.

    """
    npart = len(pos)

    filemode = 'a' if append else 'w'
    with open(filename, filemode, encoding='utf-8') as output_file:
        output_file.write(f'{npart}\n')
        header = ['#']
        if step is not None:
            header.append(f'Step: {step}')
        if box is not None:
            header.append(f'Box: {" ".join([f"{i:9.4f}" for i in box])}')
        header.append('\n')
        header_str = ' '.join(header)
        output_file.write(header_str)
        for i in range(npart):
            line = _XYZ_BIG_VEL_FMT.format(names[i], pos[i, 0], pos[i, 1],
                                           pos[i, 2], vel[i, 0], vel[i, 1],
                                           vel[i, 2])
            output_file.write(f'{line}\n')


def xyz_merge(backward, forward, merged):
    """Merge forward and backward trajectories into one.

    Parameters
    ----------
    backward : string
        File path to the backward trajectory.
    forward : string
        File path to the forward trajectory.
    merged : string
        File path to the merged (output) trajectory. This will
        overwrite existing files with the same file path!

    Note
    ----
    The velocities in the backward trajectory will not be reversed,
    the purpose of this method is mainly to make whole trajectories
    for visualisation purposes.

    """
    reverse_xyz_file(backward, merged)
    with open(forward, 'r', encoding='utf-8') as infile, \
            open(merged, 'a+', encoding='utf-8') as output:
        for lines in infile:
            output.write(lines)


def _reverse_xyz_buffer(buff, output):
    """Reverse the order in a XYZ buffer and extract frames.

    Parameters
    ----------
    buff : string
        A buffer read from a .xyz file.
    output : object like ``file object``
        Where we write the reversed frames.

    Returns
    -------
    out : string
        The left-over data which did not make a full frame.

    """
    frame = []
    finish_next = False
    for lines in reversed(buff.split('\n')):
        if not lines:
            continue
        frame.append(lines)
        if finish_next:
            output.write('\n'.join(frame[::-1]))
            output.write('\n')
            output.flush()
            frame = []
            finish_next = False
            continue
        if lines.find('# Step') != -1:
            # need just one more line to complete the frame
            finish_next = True
    if not frame:
        return None
    return '\n'.join(frame[::-1])


def reverse_xyz_file(filename, outputfile):
    """Reverse the order for frames in a given XYZ-file.

    Parameters
    ----------
    filename : string
        The input .xyz file to open.
    outputfile : string
        The .xyz file to write.

    Note
    ----
    This method will *NOT* reverse velocities.

    """
    buff_size = io.DEFAULT_BUFFER_SIZE
    left_over = None
    with open(filename, 'r', encoding='utf-8') as fileh, \
            open(outputfile, 'w', encoding='utf-8') as outfh:
        fileh.seek(0, 2)  # Go to the end
        current_pos = fileh.tell()
        done = False
        while current_pos >= 0 and not done:
            next_pos = current_pos - buff_size
            if next_pos > 0:
                fileh.seek(next_pos)
                buff = fileh.read(buff_size)
                if left_over is not None:
                    buff += left_over
                left_over = _reverse_xyz_buffer(buff, outfh)
                current_pos = fileh.tell() - buff_size
            else:
                # data < buff_size is left, just read it all:
                fileh.seek(0)
                buff = fileh.read(current_pos)
                if left_over is not None:
                    buff += left_over
                left_over = _reverse_xyz_buffer(buff, outfh)
                done = True


def txt_to_xyz(inputfile, outputfile, atoms, selection=None, nzero=6):
    """Convert a .txt trajectory to a .xyz trajectory.

    The input trajectory is a trajectory in the internal
    trajectory format.

    Parameters
    ----------
    inputfile : string
        The input trajectory file.
    outputfile : string
        A template name for output trajectories.
    atoms : list of strings
        Atom names to write to the XYZ-file.
    selection : string, optional
        The selection can be used to select trajectories for output
        based on the status.
    nzero : integer, optional
        The number of zeros we use to pad trajectory names.

    Returns
    -------
    out : list of strings
        The files we created.

    """
    all_output = []
    out = ''.join([outputfile, '-{:0', f'{nzero}d', '}-{}.xyz'])
    for traj in PathIntFormatter().load(inputfile):
        split = traj['comment'][0].split()
        cycle = int(split[2][:-1])
        status = split[4]
        if selection is not None and status.lower() != selection.lower():
            continue
        output = out.format(cycle, status)
        all_output.append(output)
        with open(output, 'w', encoding='utf-8') as outh:
            for j, snapshot in enumerate(traj['data']):
                for lines in format_xyz_data(snapshot['pos'],
                                             vel=snapshot['vel'],
                                             names=atoms,
                                             header=f'Snapshot: {j}'):
                    outh.write(f'{lines}\n')
    return all_output
