# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Define a factory for potentials."""
import logging
from pyretis.core.common import generic_factory
from pyretis.forcefield.potentials import (PairLennardJonesCut,
                                           PairLennardJonesCutnp,
                                           DoubleWellWCA,
                                           DoubleWell,
                                           RectangularWell)
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = ['potential_factory']


def potential_factory(settings):
    """Create a potential according to the given settings.

    This function is included as a convenient way of setting up and
    selecting a potential function.

    Parameters
    ----------
    settings : dict
        This defines how we set up and select the potential.

    Returns
    -------
    out[0] : object like :py:class:`.PotentialFunction`
        This object represents the potential.

    """
    potential_map = {'doublewell': {'cls': DoubleWell},
                     'rectangularwell': {'cls': RectangularWell},
                     'pairlennardjonescut': {'cls': PairLennardJonesCut},
                     'pairlennardjonescutnp': {'cls': PairLennardJonesCutnp},
                     'doublewellwca': {'cls': DoubleWellWCA}}
    return generic_factory(settings, potential_map, name='potential')
