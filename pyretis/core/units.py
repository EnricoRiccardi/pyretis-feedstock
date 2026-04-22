# *- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
r"""This module defines natural constants and unit conversions.

This module defines some natural constants and conversions between units
which can be used by the PyRETIS program.
The :ref:`natural constants <natural-constants>` are mainly used for
conversions but it is also used to define the Boltzmann constant for
internal use in PyRETIS. The :ref:`unit conversions <unit-conversions>`
are mainly useful for input/output.

All numerical values are from the National Institute of Standards and
Technology and can be accessed through a web interface
http://physics.nist.gov/constants or in plain text. [7]_

Internally, all computations are carried out in units which are defined
by a length scale, an energy scale and a mass scale. The time scale
is defined by these choices.

Charges are typically given (in the input) in units of the electron
charge. The internal unit for charge is not yet implemented, but a
choice is here to include the
factor :math:`\frac{1}{\sqrt{4\pi\varepsilon_0}}`. An internal
calculation of :math:`q_1 q_2` will then include coulombs constant in
the correct units.

The different sets of unit systems are described below in
the section on :ref:`unit systems <unit-conversions-systems>`.


.. _natural-constants:

Natural constants
~~~~~~~~~~~~~~~~~
The keys for ``CONSTANTS`` defines the natural constant and its units,
for instance ``CONSTANTS['kB']['J/K']`` is the Boltzmann constants in
units of Joule per Kelvin. The currently defined natural constants are:

- ``kB`` : The Boltzmann constant. [1]_

- ``NA`` : The Avogadro constant. [2]_

- ``e`` : The elementary charge. [3]_

- ``c0`` : The velocity of light in vacuum. [4]_

- ``mu0``: Vacuum permeability. [5]_

- ``e0``: Vacuum permittivity (or permittivity of free space or electric
  constant). [6]_


.. _unit-conversions:

Unit conversions
~~~~~~~~~~~~~~~~
For defining the different unit conversions a simple set of base
conversions are defined. These represent some common units that are
convenient for input and output. For each dimension [12]_ we define some
units and the conversion between these. The base units are:

Charge:
  * ``e``: Electron charge.
  * ``C``: Coulomb.

Energy:
  * ``kcal``: Kilocalorie.
  * ``kcal/mol``: Kilocalorie per mol.
  * ``J``: Joule.
  * ``J/mol``: Joule per mol.
  * ``kJ/mol``: Kilojoule per mol.
  * ``eV``: Electronvolt.
  * ``hartree``: Hartree (atomic unit of energy).

Force:
  * ``N``: Newton.
  * ``pN``: Piconewton.
  * ``dyn``: Dyne.

Length:
  * ``A``: Ångström.
  * ``nm``: Nanometre.
  * ``bohr``: Bohr radius.
  * ``m``: Meter.

Mass:
  * ``g/mol``: Grams per mol, numerically equal to the atomic mass unit.
  * ``g``: Gram.
  * ``kg``: Kilogram.

Pressure:
  * ``Pa``: Pascal.
  * ``bar``: Bar.
  * ``atm``: Atmosphere.

Temperature:
  * ``K``: Kelvin.

Time:
  * ``s``: Second.
  * ``ps``: Picosecond.
  * ``fs``: Femtosecond
  * ``ns``: Nanosecond.
  * ``us``: Microsecond.
  * ``ms``: Millisecond.

Velocity:
  * ``m/s``: Meter per second.
  * ``nm/ps``: Nanometer per picosecond.
  * ``A/fs``: Ångström per femtosecond.
  * ``A/ps``: Ångström per picosecond.


.. _unit-conversions-systems:

Unit conversions and internal systems of units
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following system of units are defined for PyRETIS:

- ``lj``: A Lennard-Jones type of units.

- ``real``: A system of units similar to the LAMMPS unit real. [8]_

- ``metal``: A system of units similar to the LAMMPS unit metal. [8]_

- ``au``: Atomic units. [9]_

- ``electron``: A system of units similar to the LAMMPS unit
  electron. [8]_

- ``si``: A system of units similar to the LAMMPS unit si. [8]_

- ``gromacs``: A system of units similar to the units used
  by GROMACS. [10]_

- ``reduced``: A reduced system of units.

- ``cp2k``: A system of units consistent with CP2K.


The defining units for the Lennard-Jones units (``lj``) are typically
based on the Lennard-Jones parameters for one of the components, e.g.
:math:`\varepsilon`, :math:`\sigma` and the atomic mass
of argon (119.8 kB, 3.405 Å, 39.948 g/mol). [11]_ The defining
units for the other systems are given in the table below:


.. table:: Defining units for energy systems

  +-------------+--------------+-------------+--------------------+
  | System name | Energy unit  | Length unit | Mass unit          |
  +=============+==============+=============+====================+
  | real        | 1 kcal/mol   | 1 Å         |  1 g/mol           |
  +-------------+--------------+-------------+--------------------+
  | metal       | 1 eV         | 1 Å         |  1 g/mol           |
  +-------------+--------------+-------------+--------------------+
  | au          | 1 hartree    | 1 bohr      |  9.10938291e-31 kg |
  +-------------+--------------+-------------+--------------------+
  | electron    | 1 hartree    | 1 bohr      | 1 g/mol            |
  +-------------+--------------+-------------+--------------------+
  | si          | 1 J          | 1 m         | 1 kg               |
  +-------------+--------------+-------------+--------------------+
  | gromacs     | 1 kJ/mol     | 1 nm        | 1 g/mol            |
  +-------------+--------------+-------------+--------------------+
  | cp2k        | 1 hartree    | 1 Å         | 1 g/mol            |
  +-------------+--------------+-------------+--------------------+


These units are also used for the input and define the time unit.
Further, all system of units expect an input temperature in Kelvin
(``K``) and all systems, with the exception of ``si``, expect a
charge in units of electron charges. The ``si`` system uses here
Coulomb as its unit of charge. The time units for the different
energy systems are given in the table below.


.. table:: Time units for different systems

  +-------------+----------------------+
  | System name | Time unit            |
  +=============+======================+
  | real        | 48.8882129084 fs     |
  +-------------+----------------------+
  | metal       | 10.1805056505 fs     |
  +-------------+----------------------+
  | au          |  0.0241888423521 fs  |
  +-------------+----------------------+
  | electron    | 1.03274987345 fs     |
  +-------------+----------------------+
  | si          | 1.0 s                |
  +-------------+----------------------+
  | gromacs     | 1.0 ps               |
  +-------------+----------------------+
  | cp2k        | 0.0457102889683 fs   |
  +-------------+----------------------+


The interpretation here is that if you are for instance using the system
``real`` and would like to have a time step equal to 0.5 fs, then the
input time step to PyRETIS should be ``0.5 fs / 48.8882129084 fs``.
NOTE: When using external engines, PyRETIS will not do any assumptions
on the input time/length etc and simply assume that the input values
are given in correct units for the external engine. In this case, the
only time PyRETIS will make use of the unit system is when the
Boltzmann constant is used together with energies. That is, the
Boltzmann constant must be in units consistent with the energy output
from the external engine.


References and footnotes
------------------------

.. [1] https://en.wikipedia.org/wiki/Boltzmann_constant

.. [2] https://en.wikipedia.org/wiki/Avogadro_constant

.. [3] https://en.wikipedia.org/wiki/Elementary_charge

.. [4] https://en.wikipedia.org/wiki/Speed_of_light

.. [5] https://en.wikipedia.org/wiki/Vacuum_permeability

.. [6] https://en.wikipedia.org/wiki/Vacuum_permittivity

.. [7] National Institute of Standards and Technology,
   http://physics.nist.gov/cuu/Constants/Table/allascii.txt

.. [8] The LAMMPS manual, http://lammps.sandia.gov/doc/units.html

.. [9] https://en.wikipedia.org/wiki/Atomic_units

.. [10] The GROMACS manual, tables 2.1 and 2.2 on page. 8,
   http://manual.gromacs.org/documentation/5.1.1/manual-5.1.1.pdf

.. [11] Rowley et al., J. Comput. Phys., vol. 17, pp. 401-414, 1975,
   doi: http://dx.doi.org/10.1016/0021-9991

.. [12] Note that 'dimension' here is, strictly speaking, not a true
        dimension, for instance, we define conversions for the dimension
        `velocity` which in reality is composed of the dimensions
        `length` and `time`.


Examples
--------
>>> from pyretis.core.units import CONVERT  # doctest: +ELLIPSIS
>>> print(CONVERT['length'])
{('A', 'nm'): 0.1, ('A', 'bohr'): 1.8897261254578281, ('A', 'm'): 1e-10}
>>> from pyretis.core.units import create_conversion_factors
>>> create_conversion_factors('lj', length=(3.405, 'A'), energy=(119.8, 'kB'),
...                           mass=(39.948, 'g/mol'), charge='e')
>>> print(CONVERT['length']['bohr', 'nm'])
0.052917721067
>>> print(CONVERT['length']['lj', 'nm'])
0.3405
>>> print(CONVERT['length']['bohr', 'lj'])
0.155411809301
>>> create_conversion_factors('cgs', length=(0.01, 'm'),
...                           energy=(1.0e-7, 'J'),
...                           mass=(1.0, 'g'), charge='e')
>>> print(round(CONVERT['force']['cgs', 'dyn'], 2))
1.0
>>> print(round(CONVERT['time']['cgs', 's'], 2))
1.0

"""
import logging
import math
from collections import deque
import numpy as np
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())

np.set_printoptions(legacy='1.21')

__all__ = ['generate_conversion_factors', 'generate_inverse',
           'bfs_convert', 'convert_bases', 'print_table',
           'write_conversions', 'read_conversions', 'units_from_settings',
           'create_conversion_factors']


CAL = 4184.  # Define 1 kcal = `CAL` J.
CONSTANTS = {}
"""A dict containing the natural constants. Natural constants can be
accessed with `CONSTANTS['kB']['eV/K']` etc."""
CONSTANTS['kB'] = {'eV/K': 8.6173303e-5, 'J/K': 1.38064852e-23,
                   'kJ/mol/K': 8.3144598e-3}
CONSTANTS['kB']['J/mol/K'] = CONSTANTS['kB']['kJ/mol/K'] * 1000.
CONSTANTS['kB']['kJ/K'] = CONSTANTS['kB']['J/K'] / 1000.0
CONSTANTS['kB']['kcal/K'] = CONSTANTS['kB']['J/K'] / CAL
CONSTANTS['kB']['kcal/mol/K'] = CONSTANTS['kB']['J/mol/K'] / CAL
CONSTANTS['NA'] = {'1/mol': 6.022140857e23}
CONSTANTS['c0'] = {'m/s': 299792458.}
CONSTANTS['mu0'] = {'H/m': 4.0 * np.pi * 1.0e-7}
CONSTANTS['e'] = {'C':  1.6021766208e-19}
CONSTANTS['e0'] = {'F/m': 1.0 / (CONSTANTS['mu0']['H/m'] *
                                 CONSTANTS['c0']['m/s']**2)}
# Set value of kB in the system of units we have defined.
# These values will be used in the simulations.
CONSTANTS['kB']['lj'] = 1.0
CONSTANTS['kB']['reduced'] = 1.0
CONSTANTS['kB']['si'] = CONSTANTS['kB']['J/K']
CONSTANTS['kB']['real'] = CONSTANTS['kB']['kcal/mol/K']
CONSTANTS['kB']['metal'] = CONSTANTS['kB']['eV/K']
CONSTANTS['kB']['au'] = 1.0
CONSTANTS['kB']['electron'] = 3.16681534e-6
CONSTANTS['kB']['gromacs'] = CONSTANTS['kB']['kJ/mol/K']
CONSTANTS['kB']['openmm'] = CONSTANTS['kB']['kJ/mol/K']
CONSTANTS['kB']['cp2k'] = 3.16681534e-6  # hartree


DIMENSIONS = {'length', 'mass', 'time', 'energy', 'velocity', 'charge',
              'temperature', 'pressure', 'force'}
"""A dictionary with the known dimensions. Note that not all of these
are true dimensions, for instance, we are using velocity as a dimension
here. This is just because it is convenient to use this to get
conversion factors for velocities."""
# For each dimension we want conversion factors and units

CONVERT = {key: {} for key in DIMENSIONS}
"""A dictionary with conversion factors. It is used to convert between
different units, e.g. `CONVERT['length']['bohr', 'nm]` will convert
from the length unit `bohr` to the length unit `nm`."""

UNITS = {key: {} for key in DIMENSIONS}
"""A dictionary of sets. Each set defines the known base unit for a
dimension, i.e. `UNITS['length']` is the set with known base units for
the length: `UNITS['length'] = {'A', 'nm', 'bohr', 'm'}`"""

# Define a few units and some base conversions:
UNITS['length'] = {'A', 'nm', 'bohr', 'm'}
CONVERT['length']['A', 'nm'] = 0.1
# http://physics.nist.gov/cuu/Constants/Table/allascii.txt
CONVERT['length']['A', 'bohr'] = 1.0 / 0.52917721067
CONVERT['length']['A', 'm'] = 1.0e-10

UNITS['mass'] = {'g/mol', 'g', 'kg'}
CONVERT['mass']['g', 'kg'] = 1.e-3
CONVERT['mass']['kg', 'g'] = 1.e3
CONVERT['mass']['g/mol', 'g'] = 1.0 / CONSTANTS['NA']['1/mol']
CONVERT['mass']['g', 'g/mol'] = CONSTANTS['NA']['1/mol']
CONVERT['mass']['g/mol', 'kg'] = (CONVERT['mass']['g', 'kg'] /
                                  CONSTANTS['NA']['1/mol'])
CONVERT['mass']['kg', 'g/mol'] = (CONSTANTS['NA']['1/mol'] /
                                  CONVERT['mass']['g', 'kg'])

UNITS['time'] = {'s', 'ps', 'fs', 'ns', 'us', 'ms'}
CONVERT['time']['s', 'ps'] = 1.0e12
CONVERT['time']['s', 'fs'] = 1.0e15
CONVERT['time']['s', 'ns'] = 1.0e9
CONVERT['time']['s', 'us'] = 1.0e6
CONVERT['time']['s', 'ms'] = 1.0e3

UNITS['energy'] = {'kcal', 'kcal/mol', 'J', 'J/mol', 'kJ/mol', 'eV',
                   'hartree'}
CONVERT['energy']['kcal', 'kcal/mol'] = CONSTANTS['NA']['1/mol']
CONVERT['energy']['kcal', 'J'] = CAL
CONVERT['energy']['kcal', 'J/mol'] = (CONSTANTS['NA']['1/mol'] *
                                      CONVERT['energy']['kcal', 'J'])
CONVERT['energy']['kcal', 'kJ/mol'] = (CONVERT['energy']['kcal', 'J/mol'] *
                                       1.0e-3)
CONVERT['energy']['kcal', 'eV'] = (CONSTANTS['kB']['eV/K'] /
                                   CONSTANTS['kB']['kcal/K'])
CONVERT['energy']['kcal', 'hartree'] = (CONVERT['energy']['kcal', 'eV'] *
                                        (1.0 / 27.21138602))
# 27.21138602 is hartree to eV from:
# http://physics.nist.gov/cuu/Constants/Table/allascii.txt

UNITS['velocity'] = {'m/s', 'nm/ps', 'A/fs', 'A/ps'}
CONVERT['velocity']['m/s', 'nm/ps'] = 1.0e9 / 1.0e12
CONVERT['velocity']['m/s', 'A/fs'] = 1.0e10 / 1.0e15
CONVERT['velocity']['m/s', 'A/ps'] = 1.0e10 / 1.0e12

UNITS['charge'] = {'e', 'C'}
CONVERT['charge']['e', 'C'] = CONSTANTS['e']['C']
CONVERT['charge']['C', 'e'] = 1.0 / CONSTANTS['e']['C']

UNITS['pressure'] = {'Pa', 'bar', 'atm'}
CONVERT['pressure']['Pa', 'bar'] = 1.0e-5
CONVERT['pressure']['Pa', 'atm'] = 1.0 / 101325.

UNITS['temperature'] = {'K'}

UNITS['force'] = {'N', 'pN', 'dyn'}
CONVERT['force']['N', 'pN'] = 1.0e12
CONVERT['force']['N', 'dyn'] = 1.0e5

# For completeness, add self-convert:
for i in DIMENSIONS:
    for j in UNITS[i]:
        CONVERT[i][j, j] = 1.0

# Definitions for systems of units:
UNIT_SYSTEMS = {'lj': {}, 'real': {}, 'metal': {}, 'au': {},
                'electron': {}, 'si': {}, 'gromacs': {}, 'openmm': {},
                'reduced': {}}
"""A dictionary containing basic information about the different
unit systems. E.g. `UNIT_SYSTEMS['lj']['length']` contains the length
unit for the `'lj'` unit system."""
UNIT_SYSTEMS['lj'] = {'length': (3.405, 'A'),
                      'energy': (119.8, 'kB'),
                      'mass': (39.948, 'g/mol'),
                      'charge': 'e'}
UNIT_SYSTEMS['reduced'] = {'length': (1, 'A'),
                           'energy': (1, 'kB'),
                           'mass': (1.0, 'g/mol'),
                           'charge': 'e'}
UNIT_SYSTEMS['real'] = {'length': (1.0, 'A'),
                        'energy': (1.0, 'kcal/mol'),
                        'mass': (1.0, 'g/mol'),
                        'charge': 'e'}
UNIT_SYSTEMS['metal'] = {'length': (1.0, 'A'),
                         'energy': (1.0, 'eV'),
                         'mass': (1.0, 'g/mol'),
                         'charge': 'e'}
UNIT_SYSTEMS['au'] = {'length': (1.0, 'bohr'),
                      'energy': (1.0, 'hartree'),
                      'mass': (9.10938356e-31, 'kg'),
                      'charge': 'e'}
UNIT_SYSTEMS['electron'] = {'length': (1.0, 'bohr'),
                            'energy': (1.0, 'hartree'),
                            'mass': (1.0, 'g/mol'),
                            'charge': 'e'}
UNIT_SYSTEMS['si'] = {'length': (1.0, 'm'),
                      'energy': (1.0, 'J'),
                      'mass': (1.0, 'kg'),
                      'charge': 'e'}
UNIT_SYSTEMS['gromacs'] = {'length': (1.0, 'nm'),
                           'energy': (1.0, 'kJ/mol'),
                           'mass': (1.0, 'g/mol'),
                           'charge': 'e'}
UNIT_SYSTEMS['openmm'] = {'length': (1.0, 'nm'),
                          'energy': (1.0, 'kJ/mol'),
                          'mass': (1.0, 'g/mol'),
                          'charge': 'e'}
UNIT_SYSTEMS['cp2k'] = {'length': (1.0, 'A'),
                        'energy': (1.0, 'hartree'),
                        'mass': (9.10938356e-31, 'kg'),
                        'charge': 'e'}


def _add_conversion_and_inverse(conv_dict, value, unit1, unit2):
    """Add a specific conversion and it's inverse.

    This function is mainly here to ensure that we don't forget to add
    the inverse conversions.

    Parameters
    ----------
    conv_dict : dict
        This is where we store the conversion.
    value : float
        The conversion factor to add
    unit1 : string
        The unit we are converting from.
    unit2 : string
        The unit we are converting to.

    Returns
    -------
    None, but updates the given `conv_dict`.

    """
    conv_dict[unit1, unit2] = value
    conv_dict[unit2, unit1] = 1.0 / conv_dict[unit1, unit2]


def _generate_conversion_for_dim(conv_dict, dim, unit):
    """Generate conversion factors for the specified dimension.

    It will generate all conversions from a specified unit for a given
    dimension considering all other units defined in `UNITS`.

    Parameters
    ----------
    conv_dict : dict
        A dictionary with conversions which we wish to update.
    dim : string
        The dimension to consider
    unit : string
        The unit we create conversions for.

    Returns
    -------
    None, but updates the given `conv_dict`

    """
    convertdim = conv_dict[dim]
    for unit_to in UNITS[dim]:
        if unit == unit_to:  # just skip
            continue
        value = bfs_convert(convertdim, unit, unit_to)[1]
        if value is not None:
            _add_conversion_and_inverse(convertdim, value, unit, unit_to)
        else:
            logger.warning('Could not convert %s -> %s for dimension %s',
                           unit, unit_to, dim)


def generate_conversion_factors(unit, distance, energy, mass, charge='e'):
    u"""Create conversions for a system of units from fundamental units.

    This will create a system of units from the three fundamental units
    distance, energy and mass.

    Parameters
    ----------
    unit : string
        This is a label for the unit
    distance : tuple
        This is the distance unit. The form is assumed to be
        `(value, unit)` where the unit is one of the known distance
        units, 'nm', 'A', 'm'.
    energy : tuple
        This is the energy unit. The form is assumed to be
        `(value, unit)` where unit is one of the known energy
        units, 'J', 'kcal', 'kcal/mol', 'kb'.
    mass : tuple
        This is the mass unit. The form is assumed to be `(value, unit)`
        where unit is one of the known mass units, 'g/mol', 'kg', 'g'.
    charge : string, optional
        This selects the base charge. It can be 'C' or 'e' for Coulomb
        or the electron charge. This will determine how we treat
        Coulomb's constant.

    """
    CONVERT['length'][unit, distance[1]] = distance[0]
    CONVERT['mass'][unit, mass[1]] = mass[0]
    if energy[1] == 'kB':  # in case the energy is given in units of kB.
        CONVERT['energy'][unit, 'J'] = energy[0] * CONSTANTS['kB']['J/K']
    else:
        CONVERT['energy'][unit, energy[1]] = energy[0]
        # let us also check if we can define kB now:
        if unit not in CONSTANTS['kB']:
            try:
                kboltz = CONSTANTS['kB'][f'{energy[1]}/K'] / energy[0]
                CONSTANTS['kB'][unit] = kboltz
            except KeyError as err:
                msg = f'For "{unit}" you need to define kB'
                raise ValueError(msg) from err
    # First, set up simple conversions:
    for dim in ('length', 'energy', 'mass'):
        _generate_conversion_for_dim(CONVERT, dim, unit)
    # We can now set up time conversions (since it's using length, mass and
    # energy:
    try:
        value = (CONVERT['length'][unit, 'm']**2 *
                 CONVERT['mass'][unit, 'kg'] /
                 CONVERT['energy'][unit, 'J'])**0.5
    except KeyError as err:
        keys = [(unit, 'm'), (unit, 'kg'), (unit, 'J')]
        dims = ['length', 'mass', 'energy']
        msg = ""
        for key, dim in zip(keys, dims):
            if key not in CONVERT[dim]:
                msg += 'Missing {key} in CONVERT["{dim}"]. '
        raise ValueError(msg) from err
    _add_conversion_and_inverse(CONVERT['time'], value, unit, 's')
    # And velocity (since it's using time and length):
    value = CONVERT['length'][unit, 'm'] / CONVERT['time'][unit, 's']
    _add_conversion_and_inverse(CONVERT['velocity'], value, unit, 'm/s')
    # And pressure (since it's using energy and length):
    value = CONVERT['energy'][unit, 'J'] / CONVERT['length'][unit, 'm']**3
    _add_conversion_and_inverse(CONVERT['pressure'], value, unit, 'Pa')
    # And force (since it's using energy and length):
    value = CONVERT['energy'][unit, 'J'] / CONVERT['length'][unit, 'm']
    _add_conversion_and_inverse(CONVERT['force'], value, unit, 'N')
    # Generate the rest of the conversions:
    for dim in ('time', 'velocity', 'pressure', 'force'):
        _generate_conversion_for_dim(CONVERT, dim, unit)
    # Now, figure out the Temperature conversion:
    kboltz = CONSTANTS['kB']['J/K'] * CONVERT['energy']['J', unit]
    # kboltz in now in units of 'unit'/K, temperature conversion is:
    value = CONSTANTS['kB'][unit] / kboltz
    _add_conversion_and_inverse(CONVERT['temperature'], value, unit, 'K')
    # convert permittivity:
    if charge == 'C':
        CONSTANTS['e0'][unit] = CONSTANTS['e0']['F/m']
    else:
        CONSTANTS['e0'][unit] = (CONSTANTS['e0']['F/m'] *
                                 CONVERT['charge']['C', 'e']**2 /
                                 (CONVERT['force']['N', unit] *
                                  CONVERT['length']['m', unit]**2))
    value = np.sqrt(4.0 * np.pi * CONSTANTS['e0'][unit])
    _add_conversion_and_inverse(CONVERT['charge'], value, unit, charge)
    # convert [charge] * V/A to force, in case it's needed in the future:
    # qE = CONVERT['energy']['J', unit] / CONVERT['charge']['C', 'e']
    _generate_conversion_for_dim(CONVERT, 'charge', unit)


def generate_inverse(conversions):
    """Generate all simple inverse conversions.

    A simple inverse conversion is something we can obtain by doing
    a ``1 / unit`` type of conversion.

    Parameters
    ----------
    conversions : dictionary
        The unit conversions as `convert[quantity]`.
        Note that this variable will be updated by this function.

    Returns
    -------
    out : None
        Will not return anything, but will update the given parameter
        `conversions`.

    """
    newconvert = {}
    for unit in conversions:
        unit_from, unit_to = unit
        newunit = (unit_to, unit_from)
        if newunit not in conversions:
            newconvert[newunit] = 1.0 / conversions[unit]
    for newunit, value in newconvert.items():
        conversions[newunit] = value


def bfs_convert(conversions, unit_from, unit_to):
    """Generate unit conversion between the provided units.

    The unit conversion can be obtained given that a "path" between
    these units exist. This path is obtained by a Breadth-first search.

    Parameters
    ----------
    conversions : dictionary
        The unit conversion as `convert[quantity]`.
    unit_from : string
        Starting unit.
    unit_to : string
        Target unit.

    Returns
    -------
    out[0] : tuple
        A tuple containing the two units: `(unit_from, unit_to)`.
    out[1] : float
        The conversion factor.
    out[2] : list of tuples
        The 'path' between the units, i.e. the traversal from
        `unit_from` to `unit_to`. `out[2][i]` gives the
        `(unit_from, unit_to)` tuple for step `i` in the conversion.

    """
    if unit_from == unit_to:
        return (unit_from, unit_to), 1.0, None
    if (unit_from, unit_to) in conversions:
        return (unit_from, unit_to), conversions[unit_from, unit_to], None
    que = deque([unit_from])
    visited = [unit_from]
    parents = {unit_from: None}
    while que:
        node = que.popleft()
        if node == unit_to:
            break
        for unit in conversions:
            unit1, unit2 = unit
            if unit1 != node:
                continue
            if unit2 not in visited:
                visited.append(unit2)
                que.append(unit2)
                parents[unit2] = node
    path = []
    if unit_to not in parents:
        logger.warning('Could not determine conversion %s -> %s', unit_from,
                       unit_to)
        return (unit_from, unit_to), None, None
    node = unit_to
    while parents[node]:
        new = [None, node]
        node = parents[node]
        new[0] = node
        path.append(tuple(new))
    factor = 1
    for unit in path[::-1]:
        factor *= conversions[unit]
    return (unit_from, unit_to), factor, path[::-1]


def convert_bases(dimension):
    """Create all conversions between base units.

    This function will generate all conversions between base units
    defined in a `UNITS[dimension]` dictionary. It assumes that one of
    the bases have been used to defined conversions to all other bases.

    Parameters
    ----------
    dimension : string
        The dimension to convert for.

    """
    convert = CONVERT[dimension]
    # start by generating inverses
    generate_inverse(convert)
    for key1 in UNITS[dimension]:
        for key2 in UNITS[dimension]:
            if key1 == key2:
                continue
            value = bfs_convert(convert, key1, key2)[1]
            if value is not None:
                unit1 = (key1, key2)
                unit2 = (key2, key1)
                convert[unit1] = value
                convert[unit2] = 1.0 / convert[unit1]
            else:
                logger.warning(('Could not convert base %s -> %s for '
                                'dimension %s'), key1, key2, dimension)


def generate_system_conversions(system1, system2):
    """Generate conversions between two different systems.

    Parameters
    ----------
    system1 : string
        The system we convert from.
    system2 : string
        The system we convert to.

    Returns
    -------
    out : None
        Returns nothing, but updates `CONVERT`.

    """
    for dim in CONVERT:
        convert = CONVERT[dim]
        value = bfs_convert(convert, system1, system2)[1]
        if value is not None:
            convert[system1, system2] = value
            convert[system2, system1] = 1.0 / convert[system1, system2]
        else:
            logger.warning('Could not convert %s -> %s for dimension %s',
                           system1, system2, dim)


def print_table(unit, system=False):  # pragma: no cover
    """Print out tables with conversion factors.

    This is a table in rst format which is useful for displaying
    conversions.

    Parameters
    ----------
    unit : string
        The unit to print out conversions for.
    system : boolean, optional
        Determines if information for system conversions should be printed.

    Returns
    -------
    out: None
        Does not return anything, but will print to the screen.

    """
    row_fmt = '  | {:10s} | {:16.8e} | {:16.8e} |'
    row_head = '  | {:10s} | {:16s} | {:16s} |'
    row_line = ''.join(['+-', '-' * 10, '-+-', '-' * 16, '-+-',
                        '-' * 16, '-+'])
    row_line = '  {row_line}'
    if system:
        title = f'Conversions for {unit} to other systems'
        print(f'\n.. _conversions-{unit}-system:')
    else:
        title = f'Conversions for {unit}'
        print(f'\n.. _conversions-{unit}:')
    print(f'\n{title}')
    print('-' * len(title))
    for dim in sorted(CONVERT):
        header = f'.. table:: Conversion factors for: {dim.capitalize()}'
        if system:
            print(f'\n.. _{dim}-{unit}-system:')
        else:
            print(f'\n.. _{dim}-{unit}:')
        print(f'\n{header}\n')
        row = row_head.format('Unit', f'{unit} -> unit', f'unit -> {unit}')
        print(row_line)
        print(row)
        print(row_line.replace('-', '='))
        for unt in sorted(CONVERT[dim]):
            un1, un2 = unt
            if system and (un1 in UNITS[dim] or un2 in UNITS[dim]):
                continue
            if un1 == unit:
                row = row_fmt.format(un2, CONVERT[dim][unt],
                                     CONVERT[dim][un2, un1])
                print(row)
                print(row_line)
    print('\n')
    print(f'Value of kB: {CONSTANTS["kB"][unit]}')
    print('\n')


def write_conversions(filename='units.txt'):
    """Print out the information currently stored in CONVERT.

    This function is intended for creating a big list of conversion
    factors that can be included in this script for defining unit
    conversions.

    Parameters
    ----------
    filename : string, optional
        The file to write units to.

    Returns
    -------
    out : None
        Will not return anything, but writes the given file.

    """
    with open(filename, 'w', encoding="utf-8") as fileh:
        for dim in sorted(CONVERT):
            convert = CONVERT[dim]
            for unit in sorted(convert):
                out = f'{dim} {unit[0]} {unit[1]} {convert[unit]}\n'
                fileh.write(out)
                fileh.flush()


def read_conversions(filename='units.txt', select_units=None):
    """Load conversion factors from a file.

    This will load unit conversions from a file.

    Parameters
    ----------
    filename : string, optional
        The file to load units from.
    select_units : string, optional
        If `select_units` is different from None, it can be used to
        pick out only conversions for a specific system of units,
        e.g. 'real' or 'gromacs', ... etc.

    Returns
    -------
    out : dict
        A dictionary with the conversions.

    """
    convert = {}
    with open(filename, 'r', encoding="utf-8") as fileh:
        for lines in fileh:
            try:
                dim, unit1, unit2, conv = lines.strip().split()
                conv = float(conv)
                if dim not in CONVERT:
                    raise ValueError
            except ValueError:
                logger.warning('Skipping line "%s" in "%s"', lines.strip(),
                               filename)
                continue
            if dim not in convert:
                convert[dim] = {}
            if not select_units:
                convert[dim][unit1, unit2] = conv
            else:
                if select_units in (unit1, unit2):
                    convert[dim][unit1, unit2] = conv
    return convert


def _check_input_unit(unit, dim, input_unit):
    """Check input units for :py:func:`.create_conversion_factors`.

    Parameters
    ----------
    unit : string
        Name for the unit system we are dealing with.
    dim : string
        The dimension we are looking at, typically 'length', 'mass' or
        'energy'.
    input_unit : tuple
        This is the input unit on the form (value, string) where the
        value is the numerical value and the string the unit,
        e.g. `(1.0, nm)`.

    Returns
    -------
    out : tuple
        The `input_unit` if it passes the tests, otherwise an exception
        will be raised. If the `input_unit` is `None` the default values
        from `UNIT_SYSTEMS` will be returned if they have been defined.

    Raises
    ------
    ValueError
        If the unit in `input_unit` is unknown or malformed.

    """
    if input_unit is not None:
        try:
            value, unit_dim = input_unit
            if unit_dim not in UNITS[dim] and unit_dim != 'kB':
                msg = f'Unknown {dim} unit: {unit_dim}'
                raise LookupError(msg)
            return value, unit_dim
        except TypeError as err:
            msg = f'Could not understand {dim} unit: {input_unit}'
            raise TypeError(msg) from err
    else:  # Try do get values from default:
        try:
            value, unit_dim = UNIT_SYSTEMS[unit][dim]
            return value, unit_dim
        except KeyError as err:
            msg = f'Could not determine {dim} unit for {unit}'
            raise ValueError(msg) from err


def create_conversion_factors(unit, length=None, energy=None, mass=None,
                              charge=None):
    """Set up conversion factors for a system of units.

    Parameters
    ----------
    unit : string
        A name for the system of units
    length : tuple, optional
        This is the length unit given as (float, string) where the
        float is the numerical value and the string the unit,
        e.g. `(1.0, nm)`.
    energy : tuple, optional
        This is the energy unit given as (float, string) where the
        float is the numerical value and the string the unit,
        e.g. `(1.0, eV)`.
    mass : tuple, optional
        This is the mass unit given as (float, string) where the
        float is the numerical value and the string the unit,
        e.g. `(1.0, g/mol)`.
    charge : string, optional
        This is the unit of charge given as a string, e.g. 'e' or 'C'.

    Returns
    -------
    None but will update `CONVERT` so that the conversion factors are
    available.

    """
    # First just set up conversions for the base units:
    for key in DIMENSIONS:
        convert_bases(key)
    # Check inputs and generate factors:
    length = _check_input_unit(unit, 'length', length)
    energy = _check_input_unit(unit, 'energy', energy)
    mass = _check_input_unit(unit, 'mass', mass)
    if charge is None:
        try:
            charge = UNIT_SYSTEMS[unit]['charge']
        except KeyError as err:
            msg = f'Undefined charge unit for {unit}'
            raise ValueError(msg) from err
    else:
        if charge not in UNITS['charge']:
            msg = f'Unknown charge unit "{charge}" requested.'
            raise ValueError(msg)
    generate_conversion_factors(unit, length, energy, mass, charge=charge)


def _format_unit_value(value):
    """Format a unit definition value for logging."""
    if isinstance(value, tuple):
        return f'{value[0]:g} {value[1]}'
    return str(value)


def _pick_display_unit(dimension, unit):
    """Pick a readable base unit for displaying conversions."""
    best = None
    best_score = None
    for base_unit in UNITS[dimension]:
        value = CONVERT[dimension].get((unit, base_unit))
        if value is None or value <= 0.0:
            continue
        score = abs(math.log10(value))
        if best_score is None or score < best_score:
            best_score = score
            best = (value, base_unit)
    return best


def _default_units_from_engine(settings):
    """Return the default unit system defined by the selected engine.

    Parameters
    ----------
    settings : dict
        The simulation settings used to resolve the engine and its
        default unit system.

    Returns
    -------
    out[0] : string
        The default unit system.
    out[1] : string
        The resolved engine class name.

    Raises
    ------
    KeyError
        If the selected engine does not define default units.
    ValueError
        If the resolved engine returns a unit definition of the wrong
        type.

    """
    # pylint: disable=import-outside-toplevel
    from pyretis.engines import get_default_units, resolve_engine_class

    engine_class = resolve_engine_class(settings)
    engine_name = settings.get('engine', {}).get('class')
    if engine_class is not None:
        engine_name = engine_class.__name__
    unit = get_default_units(settings)
    if unit is None:
        raise KeyError('units')
    if not isinstance(unit, str):
        msg = (
            f'Expected default units from engine "{engine_name}" to be a '
            f'string, but got {type(unit).__name__}.'
        )
        raise ValueError(msg)
    return unit.lower().strip(), engine_name


def _log_units_details(settings, resolved):
    """Write detailed unit information to the log file.

    Parameters
    ----------
    settings : dict
        The simulation settings.
    resolved : dict
        Information about the resolved unit system.

    Raises
    ------
    ValueError
        If the unit source marker is unknown.

    """
    unit = resolved['unit']
    source = resolved['source']
    if source == 'system':
        logger.info('Unit system "%s" taken from system.units.', unit)
    elif source == 'unit-system':
        logger.info('Unit system "%s" taken from unit-system.name.', unit)
    elif source == 'engine':
        logger.info('Unit system "%s" selected as default for engine "%s".',
                    unit, resolved.get('engine'))
    else:
        raise ValueError(f'Unknown unit source "{source}"')

    if ('unit-system' in settings and
            settings['unit-system'].get('name', '').lower() == unit):
        definitions = {
            key: settings['unit-system'][key]
            for key in ('length', 'energy', 'mass', 'charge')
            if key in settings['unit-system']
        }
    else:
        definitions = UNIT_SYSTEMS.get(unit, {})

    time_unit = _pick_display_unit('time', unit)
    if time_unit is None:
        time_info = 'unknown'
    else:
        time_info = f'{time_unit[0]:g} {time_unit[1]}'

    logger.info(
        'Units used in calculations: length=%s, energy=%s, mass=%s, '
        'charge=%s, time=%s, temperature=K.',
        _format_unit_value(definitions.get('length', 'unknown')),
        _format_unit_value(definitions.get('energy', 'unknown')),
        _format_unit_value(definitions.get('mass', 'unknown')),
        _format_unit_value(definitions.get('charge', 'unknown')),
        time_info,
    )


def units_from_settings(settings):
    """Set up units from given input settings.

    Parameters
    ----------
    settings : dict
        A dict defining the units. The unit system is taken from
        ``system.units`` when explicitly given. Otherwise, the selected
        engine must define default units.

    Returns
    -------
    msg : string
        Just a string with some information about the units
        created. This can be used for printing out some info to
        the user.

    Raises
    ------
    KeyError
        If the ``system`` section is missing, or if no unit system can
        be determined.
    ValueError
        If the provided unit-system settings are inconsistent.

    """
    try:
        system_settings = settings['system']
    except KeyError as err:
        raise KeyError('system') from err

    resolved = {'source': 'system'}
    if 'units' in system_settings and system_settings['units'] is not None:
        unit = system_settings['units'].lower().strip()
    elif 'unit-system' in settings and settings['unit-system'].get('name'):
        unit = settings['unit-system']['name'].lower().strip()
        system_settings['units'] = unit
        resolved['source'] = 'unit-system'
    else:
        unit, engine_name = _default_units_from_engine(settings)
        unit = unit.lower().strip()
        system_settings['units'] = unit
        resolved = {'source': 'engine', 'engine': engine_name}

    resolved['unit'] = unit
    if 'unit-system' in settings:
        try:
            unit2 = settings['unit-system']['name'].lower()
        except KeyError as err:
            msg = 'Could not find "name" setting for section "unit-system"!'
            logger.error(msg)
            raise ValueError(msg) from err
        if not unit2 == unit:
            msg = f'Inconsistent unit settings "{unit}" != "{unit2}".'
            logger.error(msg)
            raise ValueError(msg)
        setts = {}
        for key in ('length', 'energy', 'mass', 'charge'):
            try:
                setts[key] = settings['unit-system'][key]
            except KeyError as err:
                msg = f'Could not find "{key}" for section "unit-system"!'
                logger.error(msg)
                raise ValueError(msg) from err
        logger.debug('Creating unit system: "%s".', unit)
        create_conversion_factors(unit, **setts)
        msg = f'Created unit system: "{unit}".'
    else:
        logger.debug('Creating units: "%s".', unit)
        create_conversion_factors(unit)
        msg = f'Created units: "{unit}".'
        if resolved['source'] == 'engine':
            msg = (
                f'Created units: "{unit}" '
                f'(default for engine "{resolved.get("engine")}").'
            )
    _log_units_details(settings, resolved)
    return msg


def _examples():  # pragma: no cover
    """Just some examples of usage."""
    # Create new system and print it out:
    new_system = {'pyretis': {'length': (10, 'A'),
                              'energy': (1000, 'J/mol'),
                              'mass': (1.0, 'g/mol'),
                              'charge': 'e'}}
    for uni, sub_dict in new_system.items():
        create_conversion_factors(
            uni,
            length=sub_dict['length'],
            energy=sub_dict['energy'],
            mass=sub_dict[uni]['mass'],
            charge=sub_dict['charge'])
        print_table(uni)
    # Also add some base conversions for systems:
    create_conversion_factors('lj')
    create_conversion_factors('gromacs')
    create_conversion_factors('real')
    create_conversion_factors('cp2k')
    # Show how we can convert between systems:
    for key in ('energy', 'time'):
        for sys1 in ('lj', 'gromacs', 'real', 'cp2k'):
            for sys2 in ('lj', 'gromacs', 'real', 'cp2k'):
                if sys1 == sys2:
                    continue
                print(f'\nTo convert "{key}" '
                      f'between systems "{sys1}" & "{sys2}"')
                _, value, path = bfs_convert(CONVERT[key], sys1, sys2)
                txt_path = [f'{nodes[0]} -> {nodes[1]}' for nodes in path]
                txt = ' -> '.join(txt_path)
                print(f'Conversion value: {value}')
                print(f'Conversion path: {txt}')
    # To generate conversions between different systems:
    for sys1 in UNIT_SYSTEMS:
        for sys2 in UNIT_SYSTEMS:
            if sys1 != sys2:
                generate_system_conversions(sys1, sys2)
    write_conversions(filename='units-example.txt')


if __name__ == '__main__':  # pragma: no cover
    _examples()
