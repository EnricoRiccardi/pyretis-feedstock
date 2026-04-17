# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Method to re-calculate order parameters for external trajectories.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

recalculate_order (:py:func:`.recalculate_order`)
    Generic method for recalculating order parameters.

recalculate_from_trj (:py:func:`.recalculate_from_trj`)
    Recalculate order parameters using a GROMACS .trr file.

recalculate_from_xyz (:py:func:`.recalculate_from_xyz`)
    Recalculate order parameters using a .xyz file.

recalculate_from_frame (:py:func:`.recalculate_from_frame`)
    Recalculate order parameters using a .gro or .g96 file.

"""
import collections
import logging
import os
import numpy as np
from pyretis.core import System, ParticlesExt
from pyretis.core.box import box_matrix_to_list
from pyretis.inout.formats.gromacs import (
    read_trr_file,
    read_gromos96_file,
    read_gromacs_gro_file
)
from pyretis.inout.formats.xyz import read_xyz_file, convert_snapshot
from pyretis.inout.formats.path import PathExtFile

logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())

__all__ = [
    'recalculate_from_trj',
    'recalculate_from_xyz',
    'recalculate_from_frame',
    'recalculate_order',
]


def recalculate_from_trj(order_parameter, trr_file, options):
    """Re-calculate order parameters from a .trr file.

    Parameters
    ----------
    order_parameter : object like :py:class:`.OrderParameter`
        The order parameter to use.
    trr_file : string
        The path to the trr file we should read.
    options: dict
        It contains:

        * `reverse`: boolean, optional
          If True, we reverse the velocities.
        * `maxidx`: integer, optional
          This is the maximum frame we will read. Can be used in case
          the .trr file contains extra frames not needed by us.
        * `minidx`: integer, optional
          This is the first frame we will read. Can be used in case we
          want to skip some frames from the .trr file.
        * `idx`: integer, optional
          This allows the selection of a single frame to recompute.

    Yields
    ------
    out : list of lists of floats
        The order parameters, calculated per frame.

    """
    system = System(box=None)  # Add dummy system.
    msg = 'Re-calculate from %s: Step %s, time %s'
    minidx, maxidx = options.get('minidx'), options.get('maxidx')
    if options.get('idx', False):
        maxidx = options['idx']
        minidx = options['idx']
    for i, (header, data) in enumerate(read_trr_file(trr_file)):
        if maxidx is not None and i > maxidx:
            break
        if minidx is not None and i < minidx:
            continue
        logger.info(msg, os.path.basename(trr_file),
                    header['step'], header['time'])
        if system.particles is None:
            system.particles = ParticlesExt(dim=data['x'].shape[1])
        system.particles.pos = data['x']
        if 'v' in data:
            if options.get('reverse', False):
                system.particles.vel = -1.0 * data['v']
            else:
                system.particles.vel = data['v']
        else:
            logger.warning('No velocities found in .trr file! Set to 0.')
            data['v'] = np.zeros_like(data['x'])
            system.particles.vel = data['v']
        length = box_matrix_to_list(data['box'])
        system.update_box(length)
        system.particles.config = (trr_file, i)
        order = order_parameter.calculate(system)
        yield order


def recalculate_from_xyz(order_parameter, traj_file, options):
    """Re-calculate order parameters from a .xyz file.

    Parameters
    ----------
    order_parameter : object like :py:class:`.OrderParameter`
        The order parameter to use.
    traj_file : string
        The path to the trajectory file we should read.
    options: dict
        It contains:

        * `reverse`: boolean, optional
          If True, we reverse the velocities.
        * `maxidx`: integer, optional
          This is the maximum frame we will read. Can be used in case
          the .trr file contains extra frames not needed by us.
        * `minidx`: integer, optional
          This is the first frame we will read. Can be used in case we
          want to skip some frames from the .trr file.
        * `box`: list of floats
          It contains the box vector lenght. It is required in the case
          that .xyz do not normally contains the simulation box dimension.

    Yields
    ------
    out : list of lists of floats
        The order parameters as a list.

    """
    system = System(box=None)
    msg = 'Re-calculate from %s: Step %s'
    minidx, maxidx = options.get('minidx'), options.get('maxidx')
    reverse = options.get('reverse')

    for i, snapshot in enumerate(read_xyz_file(traj_file)):
        if maxidx is not None and i > maxidx:
            break
        if minidx is not None and i < minidx:
            continue
        logger.info(msg, os.path.basename(traj_file), i)
        box, xyz, vel, _ = convert_snapshot(snapshot)
        if box is None:
            box = options.get('box')
        if reverse:
            vel *= -1
        if system.particles is None:
            system.particles = ParticlesExt(dim=tuple(xyz.shape)[1])
        system.particles.config = (traj_file, i)
        system.particles.pos = xyz
        system.particles.vel = vel
        system.update_box(box)
        order = order_parameter.calculate(system)
        yield order


def recalculate_from_frame(order_parameter, frame_file, options):
    """Re-calculate order parameters from a .g96/.gro file.

    Here we assume that there is *ONE* frame in the ``frame_file``.

    Parameters
    ----------
    order_parameter : object like :py:class:`.OrderParameter`
        The order parameter to use.
    frame_file : string
        The path to the frame file we should read.
    options: dict
        It contains:

        * `ext`: string
          File extension for the ``frame_file``.
        * `reverse`: boolean, optional
          If True, we reverse the velocities.

    Returns
    -------
    out : list of lists of floats
        The order parameters for the current frame.

    """
    system = System(box=None)
    msg = f'Re-calculate from {os.path.basename(frame_file)}:'
    logger.info(msg)
    if options['ext'] == '.g96':
        _, xyz, vel, box = read_gromos96_file(frame_file)
    elif options['ext'] == '.gro':
        _, xyz, vel, box = read_gromacs_gro_file(frame_file)
    else:
        raise ValueError(f"Unknown format {options['ext']}")
    if options.get('reverse'):
        vel *= -1
    if system.particles is None:
        system.particles = ParticlesExt(dim=xyz.shape[1])
    system.particles.config = (frame_file, 0)
    system.particles.pos = xyz
    system.particles.vel = vel
    system.update_box(box)
    return [order_parameter.calculate(system)]


def recalculate_order(order_parameter, traj_file, options):
    """Re-calculate order parameters.

    Parameters
    ----------
    order_parameter : object like :py:class:`.OrderParameter`
        The order parameter to use.
    traj_file : string
        Path to the trajectory file to recalculate for.
    options: dict
        It contains:

        * `reverse`: boolean, optional
          If True, we reverse the velocities.
        * `maxidx`: integer, optional
          This is the maximum frame we will read. Can be used in case
          the .trr file contains extra frames not needed by us.
        * `minidx`: integer, optional
          This is the first frame we will read. Can be used in case we
          want to skip some frames from the .trr file.

    """
    _, ext = os.path.splitext(traj_file)

    helpers = {
        '.trr': recalculate_from_trj,
        '.xtc': recalculate_from_trj,
        '.xyz': recalculate_from_xyz,
        '.g96': recalculate_from_frame,
        '.gro': recalculate_from_frame
    }

    options['ext'] = ext
    all_order = helpers[ext](order_parameter, traj_file, options)

    if options.get('reverse'):
        return reversed(list(all_order))
    return all_order


def get_traj_files(traj_file_name):
    """Read an external file and get trajectory information.

    Parameters
    ----------
    traj_file_name : string
        The file to read trajectories from.

    Returns
    -------
    files : dict
        A dictionary with the files constituting the trajectory.

    """
    files = collections.OrderedDict()
    with PathExtFile(traj_file_name, 'r') as trajfile:
        # Just get the first trajectory:
        traj = next(trajfile.load())
        for snapshot in traj['data']:
            filename = snapshot[1]
            reverse = int(snapshot[-1]) == -1
            if filename not in files:
                files[filename] = {'minidx': None, 'maxidx': None,
                                   'reverse': reverse}
            idx = int(snapshot[2])
            minidx = files[filename]['minidx']
            if minidx is None or idx < minidx:
                files[filename]['minidx'] = idx
            maxidx = files[filename]['maxidx']
            if maxidx is None or idx > maxidx:
                files[filename]['maxidx'] = idx
    return files
