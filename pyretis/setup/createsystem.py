# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This module handles the set-up of initial positions and a box.

The initial positions can either be generated on a lattice, or it can
be read from a file.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set_up_box (:py:func:`.set_up_box`)
    Create a simulation box from simulation settings.

create_initial_positions (:py:func:`.create_initial_positions`)
    Get initial positions based on settings. This will either be
    read from a file or generated on a lattice.

create_system (:py:func:`.create_system`)
    Set up a system from given settings. This method will probably
    also need to set/get the initial positions and velocities for
    the particles and set up the simulation box.

create_velocities (:py:func:`.create_velocities`)
    Create velocities from settings for a system with particles.

initial_positions_file (:py:func:`.initial_positions_file`)
    Get initial positions from a file.

initial_positions_lattice (:py:func:`.initial_positions_lattice`)
    Get initial positions by generating a lattice.

"""
import logging
import os
import numpy as np
from pyretis.tools import generate_lattice
from pyretis.core.box import create_box
from pyretis.core.system import System
from pyretis.core.particles import Particles, get_particle_type
from pyretis.core.units import CONVERT, create_conversion_factors
from pyretis.inout.settings import look_for_input_files
from pyretis.inout.formats.snapshot import read_txt_snapshots
from pyretis.inout.formats.xyz import read_xyz_file
from pyretis.inout.formats.gromacs import read_gromacs_file, read_gromos96_file
from pyretis.engines.openmm import make_pyretis_system as openmm_system
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = ['create_initial_positions', 'create_system', 'create_velocities',
           'set_up_box', 'initial_positions_file', 'initial_positions_lattice']


PERIODIC_TABLE = {'H': 1.007975, 'He': 4.002602, 'Li': 6.9675,
                  'Be': 9.0121831, 'B': 10.8135, 'C': 12.0106,
                  'N': 14.006855, 'O': 15.9994, 'F': 18.998403163,
                  'Ne': 20.1797, 'Na': 22.98976928, 'Mg': 24.3055,
                  'Al': 26.9815385, 'Si': 28.085, 'P': 30.973761998,
                  'S': 32.0675, 'Cl': 35.4515, 'Ar': 39.948,
                  'K': 39.0983, 'Ca': 40.078, 'Sc': 44.955908,
                  'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961,
                  'Mn': 54.938044, 'Fe': 55.845, 'Co': 58.933194,
                  'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.38,
                  'Ga': 69.723, 'Ge': 72.63, 'As': 74.921595,
                  'Se': 78.971, 'Br': 79.904, 'Kr': 83.798,
                  'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90584,
                  'Zr': 91.224, 'Nb': 92.90637, 'Mo': 95.95,
                  'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42,
                  'Ag': 107.8682, 'Cd': 112.414, 'In': 114.818,
                  'Sn': 118.71, 'Sb': 121.76, 'Te': 127.6,
                  'I': 126.90447, 'Xe': 131.293, 'Cs': 132.90545196,
                  'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116,
                  'Pr': 140.90766, 'Nd': 144.242, 'Sm': 150.36,
                  'Eu': 151.964, 'Gd': 157.25, 'Tb': 158.92535,
                  'Dy': 162.5, 'Ho': 164.93033, 'Er': 167.259,
                  'Tm': 168.93422, 'Yb': 173.045, 'Lu': 174.9668,
                  'Hf': 178.49, 'Ta': 180.94788, 'W': 183.84,
                  'Re': 186.207, 'Os': 190.23, 'Ir': 192.217,
                  'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.592,
                  'Tl': 204.3835, 'Pb': 207.2, 'Bi': 208.9804,
                  'Th': 232.0377, 'Pa': 231.03588, 'U': 238.02891}


# Variable that defines some information for reading input files:
# TODO: Move this to the place where the file readers are defined?
READFILE = {'xyz': {'reader': read_xyz_file,
                    'units': {'length': 'A', 'velocity': 'A/fs'}},
            'txt': {'reader': read_txt_snapshots,
                    'units': None},
            'gro': {'reader': read_gromacs_file,
                    'units': {'length': 'nm', 'velocity': 'nm/ps'}},
            'g96': {'reader': read_gromos96_file,
                    'units': {'length': 'nm', 'velocity': 'nm/ps'}}}


def list_get(input_list, index):
    """Get an item from a list and handle out-of bounds errors.

    This method is intended to be used when we are picking items from
    a list and possibly we want a number of items which is larger than
    the number of items in the list. Here, we then just return the last
    element.

    Parameters
    ----------
    input_list : list
        The list to pick from.
    index : integer
        The index to pick.

    """
    try:
        return input_list[index]
    except IndexError:
        return input_list[-1]


def _assign_mass_from_file(filename, unit):
    """Assign masses and imasses to particles from a configuration file.

    This function is useful when velocities will be internally modified
    while using an external engine.

    Parameters
    ----------
    filename : string
        The configuration file.
    unit : string
        The system of units.

    """
    snapshot, _ = _get_snapshot_from_file({'input_file': filename}, unit)
    masses = []
    for i, atomname in enumerate(snapshot['atomname']):
        masses.append([guess_particle_mass(i+1, atomname, unit)*1.])
    return masses, np.reciprocal(masses)


def guess_particle_mass(particle_no, particle_type, unit):
    """Guess a particle mass from it's type.

    Parameters
    ----------
    particle_no : integer
        Just used to identify the particle number.
    particle_type : string
        Used to identify the particle.
    unit : string
        The system of units. This is used in case we try to get the
        mass from the periodic table where the units are in `g/mol`.

    """
    logger.info(('Mass not specified for particle no. %i\n'
                 'Will guess from particle type "%s"'), particle_no,
                particle_type)
    mass = PERIODIC_TABLE.get(particle_type, None)
    if mass is None:
        particle_mass = 1.0
        logger.info(('-> Could not find mass. '
                     'Assuming %f (internal units)'), particle_mass)
    else:
        if ('g/mol', unit) not in CONVERT['mass']:
            create_conversion_factors(unit)
        particle_mass = CONVERT['mass']['g/mol', unit] * mass
        logger.info(('-> Using a mass of %f g/mol '
                     '(%f in internal units)'), mass, particle_mass)
    return particle_mass


def initial_positions_lattice(settings):
    """Generate initial positions based on given settings.

    We assume here the input values are given with the correct units
    as dictated by ``settings['system']['units']``.

    Parameters
    ----------
    settings : dict
        The input settings for the simulation.

    Returns
    -------
    particles : object like :py:class:`.Particles`
        The particles we created.
    size : list of floats
        A size for the region we created. This can be used to create
        a box.

    """
    pos_settings = settings['particles']['position']
    ptype = settings['particles'].get('ptype', [0])
    pname = settings['particles'].get('name', ['Ar'])

    lattice, size = generate_lattice(
        pos_settings['generate'].lower(),
        pos_settings.get('repeat'),
        lcon=pos_settings.get('lcon', None),
        density=pos_settings.get('density', None)
    )

    ndim = settings['system'].get('dimensions', None)
    if ndim is None:
        ndim = len(size)
    else:
        if ndim != len(size):
            msgtxt = ('Inconsistent dimenaions in settings and generated '
                      f'lattice.\nSettings gives {ndim}D while the generated '
                      f'lattice is {len(size)}D.')
            logger.error(msgtxt)
            raise ValueError(msgtxt)

    particles = Particles(dim=ndim)

    for i, pos in enumerate(lattice):
        particle_type = list_get(ptype, i)
        particle_name = list_get(pname, i)
        # Infer the mass from the input masses, or try to get it
        # from the periodic table:
        p_mass = settings['particles'].get('mass', {}).get(particle_name)
        if p_mass is None:
            p_mass = guess_particle_mass(i + 1, particle_name,
                                         settings['system']['units'])
        particles.add_particle(
            pos,
            np.zeros_like(pos),
            np.zeros_like(pos),
            mass=p_mass,
            name=particle_name,
            ptype=particle_type
        )
    logger.info('Generated %i particles on lattice "%s".',
                particles.npart, pos_settings['generate'].lower())
    logger.info('Lattice is %iD.', ndim)

    size = np.array(size)
    box = {'low': size[:, 0], 'high': size[:, 1]}
    return particles, box


def _get_snapshot_from_file(pos_settings, units):
    """Get a configuration snapshot from a file.

    This snapshot will be used to set up the initial configuration.

    Parameters
    ----------
    pos_settings : dict
        A dict with information on what we should read.
    units : string
        The internal units.

    Returns
    -------
    snapshot : dict
        The snapshot we found in the file. It will at least have the
        keys with the positions ('x', 'y', 'z') and atom name
        'atomnames'. It may have information about velocities
        ('vx', 'vy', 'vz') and the box ('box').
    convert : dict
        Dictionary with conversion factors to internal units.

    """
    filename = pos_settings.get('input_file', None)
    if filename is None:
        msg = ('Requested reading (initial) configuration from file, '
               'but no input file given!')
        logger.error(msg)
        raise ValueError(msg)
    fmt = pos_settings.get('format', os.path.splitext(filename)[1][1:])
    snaps = []
    convert = None
    if fmt not in READFILE:
        msg = (f'Input configuration "{filename}" has unknown '
               f'format "{fmt}".')
        logger.error(msg)
        logger.error('Supported formats are: %s.', list(READFILE))
        raise ValueError(msg)

    reader = READFILE[fmt]['reader']
    read_units = READFILE[fmt]['units']
    if read_units is None:
        convert = {'length': 1.0, 'velocity': 1.0}
    else:
        convert = {
            'length': CONVERT['length'][read_units['length'], units],
            'velocity': CONVERT['velocity'][read_units['velocity'], units]
        }
    logger.info(
        'Reading initial configuration from "%s" (format: "%s").',
        filename,
        fmt,
    )
    snaps = [snap for snap in reader(filename)]

    snapshot = None
    if len(snaps) == 1:
        snapshot = snaps[0]
    elif len(snaps) > 1:
        msg = ('Found several frames ({}) in input file.'
               ' Will use the last one!').format(len(snaps))
        logger.warning(msg)
        snapshot = snaps[-1]
    else:
        msg = f"Could not find any configurations in input file: {filename}"
        logger.error(msg)
        raise ValueError(msg)
    return snapshot, convert


def initial_positions_file(settings):
    """Get initial positions from an input file.

    Parameters
    ----------
    settings : dict
        The input settings for the simulation.

    Returns
    -------
    particles : object like :py:class:`.Particles`
        The particles we created.
    size : list of floats
        A size for the region we created. This can be used to create
        a box.
    vel_read : boolean
        True if we read velocities from the input file.

    """
    ndim = settings['system'].get('dimensions', 3)
    pos_settings = settings['particles']['position']
    ptype = settings['particles'].get('ptype', None)
    pname = settings['particles'].get('name', None)
    pmass = settings['particles'].get('mass', {})
    ptypes = {}  # To automatically set particle types based on name.
    snapshot, convert = _get_snapshot_from_file(pos_settings,
                                                settings['system']['units'])
    vel_read = False
    particles = Particles(dim=ndim)
    for i, atomname in enumerate(snapshot['atomname']):
        pos = []
        vel = []
        for key in ['x', 'y', 'z'][:ndim]:
            pos.append(snapshot[key][i])
            vel_key = f'v{key}'
            if vel_key in snapshot:
                vel.append(snapshot[vel_key][i])
        pos = np.array(pos) * convert['length']
        if len(vel) != ndim:
            vel = np.zeros_like(pos)
        else:
            vel = np.array(vel) * convert['velocity']
            vel_read = True
        # Get particle type from the atom names or from input list:
        if ptype is None:
            if atomname not in ptypes:
                ptypes[atomname] = len(ptypes)
            particle_type = ptypes[atomname]
        else:
            particle_type = list_get(ptype, i)
        if pname is None:
            particle_name = atomname
        else:
            particle_name = list_get(pname, i)
        # Infer the mass from the input masses, or try to get it
        # from the periodic table:
        try:
            particle_mass = pmass[particle_name]
        except KeyError:
            particle_mass = guess_particle_mass(i + 1, particle_name,
                                                settings['system']['units'])
        particles.add_particle(pos,
                               vel, np.zeros_like(pos),
                               mass=particle_mass, name=particle_name,
                               ptype=particle_type)
    try:
        box = {'cell': [i * convert['length'] for i in snapshot['box']]}
        if ndim < 3:
            box['cell'] = box['cell'][:ndim]
    except (KeyError, IndexError, TypeError) as err:
        logger.debug('No box read from file: %s.', err)
        box = None
    logger.info('Read %d particle(s) from "%s".', particles.npart,
                pos_settings['input_file'])
    if vel_read:
        logger.info('Read velocities from file: "%s".',
                    pos_settings['input_file'])
    return particles, box, vel_read


def create_initial_positions(settings):
    """Set up the initial positions from the given settings.

    The settings can specify the initial positions as a file or
    to be generated on a lattice by PyRETIS.

    Parameters
    ----------
    settings : dict
        Settings for creating the initial positions.

    Returns
    -------
    out[0] : object like :py:class:`.Particles`
        The particles we created.
    out[1] : list
        The size associated with the particles. Can be used to create a
        box.
    out[2] : boolean
        True if we have read/created velocities different from just
        zeros. This is only True if we have read from a file with
        velocities.

    """
    logger.debug('Settings used for initial positions: %s',
                 settings['particles']['position'])
    particles = None
    if 'generate' in settings['particles']['position']:
        particles, box = initial_positions_lattice(settings)
        return particles, box, False
    if 'input_file' in settings['particles']['position']:
        # First check if we need to add a path to the file:
        filename = settings['particles']['position']['input_file']
        if (not os.path.isfile(filename) and
                'exe_path' in settings['simulation']):
            filename = os.path.join(settings['simulation']['exe_path'],
                                    filename)
            settings['particles']['position']['input_file'] = filename
        particles, box, vel = initial_positions_file(settings)
        return particles, box, vel
    msg = 'Unknown settings for initial positions: {}'
    msgtxt = msg.format(settings['particles']['position'])
    logger.error(msgtxt)
    raise ValueError(msgtxt)


def set_up_box(settings, boxs, dim=3):
    """Set up a box from given settings.

    Parameters
    ----------
    settings : dict
        The dict with the simulation settings.
    boxs : dict or None
        If no box settings are given, we can still create a box,
        inferred from the positions of the particles. This dict
        contains the settings to do so.
    dim : integer, optional
        Number of dimensions for the box. This is used only as a last
        resort when no information about the box is given.

    Returns
    -------
    box : object like :py:class:`.BoxBase` or None
        The box if we managed to create it, otherwise None.

    """
    msg = 'Box created {}:\n{}'
    box = None
    if settings.get('box', None) is not None:
        box = create_box(**settings['box'])
        msgtxt = msg.format('from settings', box)
        logger.info(msgtxt)
        debugtxt = f"Settings used:\n{settings['box']}"
        logger.debug(debugtxt)
    else:
        if boxs is not None:
            box = create_box(**boxs)
            msgtxt = msg.format('from initial positions', box)
            logger.info(msgtxt)
            msgwarn = 'The box was assumed periodic in all directions.'
            logger.warning(msgwarn)
        else:
            if dim > 0:
                box = create_box(periodic=[False]*dim)
                msgtxt = msg.format('without specifications', box)
                logger.info(msgtxt)
                msgwarn = 'The box was assumed nonperiodic in all directions.'
                logger.warning(msgwarn)
    return box


def create_velocities(system, settings, vel):
    """Create velocities from settings for a system.

    Parameters
    ----------
    system : object like :py:class:`.System`
        The system to create velocities for. It's needed since
        we need to know the degrees of freedom.
    settings : dict
        Settings to use for creating the velocities.
    vel : boolean
        If True, we already read velocities. They will now be
        overwritten. We just make some warnings about this.

    Returns
    -------
    out : boolean
        True if we actually generated velocities.

    """
    vel_settings = settings['particles'].get('velocity', {})
    if vel:
        logger.info(
            'Velocities read from input configuration (temperature: %6.2g).',
            system.calculate_temperature(),
        )
    if 'generate' in vel_settings:
        if vel:
            logger.warning(
                'Will generate and overwrite velocities already set.'
            )
        gen_settings = {'distribution': vel_settings['generate']}
        for key in ('seed', 'momentum', 'temperature', 'rgen'):
            try:
                gen_settings[key] = vel_settings[key]
            except KeyError:
                pass
        system.generate_velocities(**gen_settings)
        logger.info(
            'Generated new velocities with average temperature: %6.2g',
            system.calculate_temperature(),
        )
        logger.debug('Settings used for generating velocities: %s',
                     gen_settings)
        return True
    if 'scale' in vel_settings:
        target = vel_settings['scale']
        # Just set the velocities to some temperature for now.
        # The scaling is done later by calling `system.extra_setup()`.
        gen_settings = {
            'distribution': 'maxwell',
            'momentum': system.particles.npart != 1,
            'temperature': settings['system'].get('temperature', 1.0),
        }
        system.generate_velocities(**gen_settings)
        msg = f'Scaling velocities to total energy {target}'
        logger.debug(msg)
        system.post_setup.append(('rescale_velocities', (target,)))
        return True
    return False


def create_system(settings):
    """Create a system from input settings.

    Parameters
    ----------
    settings : dict
        The dict with the simulation settings.

    Returns
    -------
    system : object like :py:class:`.System`
        The system object we create here.

    """
    if 'restart' in settings:
        system = System()
        system.load_restart_info(settings['restart']['system'])
        return system

    vel = None
    if ('engine' not in settings or
            settings['engine'].get('type', 'internal') == 'internal'):
        particles, box, vel = create_initial_positions(settings)
        box = set_up_box(settings, box, dim=particles.dim)
    elif settings['engine'].get('type', None).lower() == 'openmm':
        return openmm_system(settings)
    else:
        # Engine is not None and not internal => external.
        klass = get_particle_type('external')
        particles = klass(dim=3)
        if 'input_path' in settings['engine']:
            required_file = {}
            # Get the engine input files
            if 'cp2k' in settings['engine']['class'].lower():
                ext = settings['engine'].get('cp2k_format', 'xyz')
                default_conf = f'initial.{ext}'
                required_file = {'conf': settings['engine'].get(
                    'conf', default_conf)}
            elif 'gromacs' in settings['engine']['class'].lower():
                required_file = {
                    'conf':
                        f"conf.{settings['engine'].get('gmx_format', 'gro')}"}

            input_file = look_for_input_files(settings['engine']['input_path'],
                                              required_file)
            particles.set_pos((input_file['conf'], None))
            # gromacs tests fail without this if-statement.
            # _get_snapshot_from_file() has issues with reading .g96 format.
            if settings['engine']['class'].lower() == 'cp2k':
                particles.mass, particles.imass = _assign_mass_from_file(
                    input_file['conf'],
                    settings['system']['units'])
        else:
            particles.set_pos((None, None))
        particles.set_vel(False)
        box = None

    system = System(
        temperature=settings['system']['temperature'],
        units=settings['system']['units'],
        box=box
        )
    system.particles = particles

    # figure out what to do with velocities:
    if 'position' in settings['particles'] and \
            'input_file' not in settings['particles']['position']:
        vel_gen = create_velocities(system, settings, vel)
        if not (vel_gen or vel):
            logger.warning('Velocities not created/read: Set to zero!')

    return system
