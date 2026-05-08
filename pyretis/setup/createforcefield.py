# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This module handles the creation of force fields from simulation settings.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

create_potentials (:py:func:`.create_potentials`)
    Method for creating potentials from a dictionary of settings.
    Note that this method will make use of :py:func:`.create_potential`.

create_force_field (:py:func:`.create_force_field`)
    Method to create a force field from input settings.
"""
import logging
from pyretis.setup.common import create_potential
from pyretis.forcefield import ForceField
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = ['create_force_field', 'create_potentials']


def create_potentials(settings):
    """Create potential functions from given simulations settings.

    This method will basically loop over the given potential settings
    and just run :py:func:`.create_potential` for each setting.

    Parameters
    ----------
    settings : dict
        This dictionary contains the settings for the simulation.

    Returns
    -------
    out[0] : list
        A list of potential functions.
    out[1] : list
        A list of parameters for the potential functions.

    """
    potentials = settings.get('potential', [])
    try:
        ndim = settings['system']['dimensions']
    except KeyError:
        ndim = None
    out_pot, out_par = [], []
    for i, pot_settings in enumerate(potentials):
        potential_function = create_potential(settings, pot_settings)
        if potential_function is None:
            msg = 'The following potential settings were ignored!\n{}'
            msgtxt = msg.format(pot_settings)
            logger.warning(msgtxt)
        pdim = getattr(potential_function, 'dim', None)
        if pdim is not None and ndim is not None:
            if ndim != pdim:
                msg = ('Inconsistent dimensions in potential!'
                       '\nSettings gives: {}D, potential {} is {}D')
                msgtxt = msg.format(ndim, i, pdim)
                logger.error(msgtxt)
                raise ValueError(msgtxt)
        out_pot.append(potential_function)
        out_par.append(pot_settings.get('parameter', None))
    return out_pot, out_par


def create_force_field(settings):
    """Create a force field from input settings.

    This method will create the required potential functions with the
    specified parameters from `settings`.

    Parameters
    ----------
    settings : dict
        This dictionary contains the settings for a single potential.

    Returns
    -------
    out : object like :py:class:`.ForceField`
        This object represents the force field.

    """
    try:
        desc = settings['forcefield']['description']
    except KeyError:
        desc = 'Generic force field'
    potentials, pot_param = create_potentials(settings)
    ffield = ForceField(desc, potential=potentials, params=pot_param)
    msg = ['Created force field:', f'{ffield}']
    msgtxt = '\n'.join(msg)
    logger.info(msgtxt)
    return ffield
