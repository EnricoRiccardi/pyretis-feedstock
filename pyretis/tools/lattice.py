# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Some methods for generating initial lattice structures.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

generate_lattice (:py:func:`.generate_lattice`)
    Generate points on a simple lattice.

Examples
--------
>>> from pyretis.tools.lattice import generate_lattice
>>> xyz, size = generate_lattice('diamond', [1, 1, 1], lcon=1)

"""
import logging
import itertools
import numpy as np
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name

__all__ = ['generate_lattice']


UNIT_CELL = {'sc': np.array([[0.0, 0.0, 0.0]]),
             '1d': np.array([[0.0]]),
             'sq': np.array([[0.0, 0.0]]),
             'sq2': np.array([[0.0, 0.0], [0.5, 0.5]]),
             'bcc': np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
             'fcc': np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0],
                              [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]]),
             'hcp': np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0],
                              [0.5, 5.0/6.0, 0.5], [0.0, 1.0/3.0, 0.5]]),
             'diamond': np.array([[0.0, 0.0, 0.0], [0.0, 0.5, 0.5],
                                  [0.5, 0.0, 0.5], [0.5, 0.5, 0.0],
                                  [0.25, 0.25, 0.25], [0.25, 0.75, 0.75],
                                  [0.75, 0.25, 0.75], [0.75, 0.75, 0.25]])}


def generate_lattice(lattice, repeat=None, lcon=None, density=None):
    """Generate points on a simple lattice.

    The lattice is one of the defined keys in the global variable
    `UNIT_CELL`. This lattice will be repeated a number of times.
    The lattice spacing can be given explicitly, or it can be given
    implicitly by the number density.

    Parameters
    ----------
    lattice : string
        Select the kind of lattice. The following options are currently
        defined in `UNIT_CELL`:

         * `1d` : 1D lattice
         * `sc` : Simple cubic lattice.
         * `sq` : Square lattice (2D) with one atom in the unit cell.
         * `sq2` : Square lattice with two atoms in the unit cell.
         * `bcc` : Body-centred cubic lattice.
         * `fcc` : Face-centred cubic lattice.
         * `hcp` : Hexagonal close-packed lattice.
         * `diamond` : Diamond structure.

    repeat : list of integers, optional.
        How many time the cell is replicated.
    lcon : float, optional
        The lattice constant.
    density : float, optional
        The desired density. If this is given, `lcon` is calculated.
        Note that density will be interpreted as given in internal
        units.

    Returns
    -------
    positions : numpy.array
        The lattice positions.
    size : list of floats
        The corresponding size(s), can be used to define a simulation
        box.

    """
    try:
        unit_cell = UNIT_CELL[lattice.lower()]
    except KeyError as err:
        msg = [f'Unknown lattice "{lattice}" requested!']
        msg += [f'Input lattice should be a string in {UNIT_CELL.keys()}']
        raise ValueError('\n'.join(msg)) from err
    except AttributeError as err:
        msgtxt = f'Input lattice should be a string in {UNIT_CELL.keys()}'
        raise ValueError(msgtxt) from err
    ndim = len(unit_cell[0])
    npart = len(unit_cell)
    if density is not None:
        lcon = (npart / density)**(1.0 / float(ndim))
    if repeat is None and density is None and lcon is None:
        logger.debug('The construction of the simulation '
                     'box is attempted with insufficient '
                     'information. Please check')
        repeat = [1]*ndim
        lcon = 1
    elif lcon is None:
        msgtxt = 'Could not determine lattice constant!'
        raise ValueError(msgtxt)
    elif len(repeat) < ndim:
        msgtxt = 'To few "repeat" values given: Expected {} but got {}.'
        raise ValueError(msgtxt.format(ndim, len(repeat)))
    positions = []
    for i in itertools.product(*[range(nri) for nri in repeat[:ndim]]):
        pos = lcon * (np.array(i) + unit_cell)
        positions.extend(pos)
    size = [[0.0, i * lcon] for i in repeat[:ndim]]
    positions = np.array(positions)
    return positions, size
