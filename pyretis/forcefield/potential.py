# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Define the class for a generic potential function.

This module defines the generic class for potential functions.
This class is sub-classed in all potential functions.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PotentialFunction (:py:class:`.PotentialFunction`)
    A class for representing generic potential functions.
"""
import logging
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = ['PotentialFunction']


class PotentialFunction:
    """Base class for a generic potential function.

    Generic class for potential functions.

    Attributes
    ----------
    desc : string
        Short description of the potential.
    dim : int
        Represents the spatial dimensionality of the potential.
    params : dict
        The parameters for the potential. This dict defines,
        on initiation, the parameters the potential will handle
        and store.

    """

    def __init__(self, dim=1, desc=''):
        """Initialise the potential.

        Parameters
        ----------
        dim : int, optional
            Represents the dimensionality.
        desc : string, optional
            Description of the potential function. Used to print out
            information about the potential.

        """
        self.dim = dim
        self.desc = desc
        self.params = {}

    def set_parameters(self, parameters):
        """Update all parameters. Input is assumed to be a dict."""
        for key in parameters:
            if key in self.params:
                self.params[key] = parameters[key]
            else:
                msg = 'Could not find "{}" in parameters. Ignoring!'
                msg = msg.format(key)
                logger.warning(msg)
        self.check_parameters()

    def check_parameters(self):
        """Check the consistency of the parameters.

        Returns
        -------
        out : boolean
            True if the check(s) pass.

        """
        if not self.params:
            logger.warning('No parameters are set for the potential')
            return False
        return True

    def __str__(self):
        """Return the string description of the potential."""
        msg = [f'Potential: {self.desc}']
        strinfo = '{}: {}'
        for key in sorted(self.params):
            msg.append(strinfo.format(key, self.params[key]))
        return '\n'.join(msg)
