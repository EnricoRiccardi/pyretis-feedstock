# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This file contains a WCA double well potential.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DoubleWellWCA (:py:class:`.DoubleWellWCA`)
    This class defines a double well WCA potential.
"""
import logging
import numpy as np
from pyretis.forcefield.potential import PotentialFunction
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = ['DoubleWellWCA']


class DoubleWellWCA(PotentialFunction):
    r"""A double well potential.

    This class defines a double well WCA potential. The potential energy
    (:math:`V_\text{pot}`) for a pair of particles separated by a
    distance :math:`r` is given by,

    .. math::

       V_\text{pot} = h (1 - (r - r_0 - w)^2/w^2)^2,

    where :math:`h` gives the 'height' of the potential, :math:`r_0` the
    minimum and :math:`w` the width. These parameters are stored in the
    attributes `height`, `rzero` and `width` respectively.

    Attributes
    ----------
    params : dict
        Contains the parameters. These are:

        * `height`: A float describing the "height" of the potential.

        * `height4`: A float equal to ``4.0 * height``.
          (This variable is just included for convenience).

        * `rzero`: A float defining the two minimums. One is located at
          ``rzero``, the other at ``rzero+2*width``.

        * `types`: A set defining what kind of particle pairs to
          consider for this interaction. If `types` is not set (i.e.
          equal to None), it will be assumed to apply to **ALL**
          particles.

        * `width`: A float describing the "width" of the potential.

        * `width2`: A float equal to ``width*width`` (for convenience).

    """

    def __init__(self, dim=3, desc='A WCA double well potential'):
        """Initialise the Double Well WCA potential.

        Parameters
        ----------
        dim : int, optional
            The dimensionality of the potential.
        desc : string, optional
            Description of the force field.

        """
        super().__init__(dim=dim, desc=desc)
        self.params = {'height': 0.0,
                       'height4': 0.0,
                       'rwidth': 0.0,
                       'rzero': 0.0,
                       'types': [],
                       'width': 0.0,
                       'width2': 0.0}

    def set_parameters(self, parameters):
        """Add new potential parameters to the potential.

        Parameters
        ----------
        parameters : dict
            The new parameters, they are assume to be dicts on the form
            ``{'types': set([(0,0)]), 'rzero': 1.0, 'width': 0.25,
            'height': 6.0}``

        """
        for key in parameters:
            if key in self.params:
                self.params[key] = parameters[key]
            else:
                msg = f'Unknown parameter {key} - ignored!'
                logger.warning(msg)
        if self.params['types'] is not None:
            self.params['types'] = set(self.params['types'])
        self.params['width2'] = self.params['width']**2
        self.params['rwidth'] = self.params['rzero'] + self.params['width']
        self.params['height4'] = 4.0 * self.params['height']

    def _activate(self, itype, jtype):
        """Determine if we should calculate a interaction or not.

        Parameters
        ----------
        itype : string
            Particle type for particle i.
        jtype : string
            Particle type for particle j.

        """
        if self.params['types'] is None:
            return True
        pair1, pair2 = (itype, jtype), (jtype, itype)
        return (pair1 in self.params['types'] or
                pair2 in self.params['types'])

    def min_max(self):
        """Return the minima & maximum of the `DoubleWellWCA` potential.

        The minima are located at ``rzero`` & ``rzero + 2*width``.
        The maximum is located at ``rzero + width``.

        Returns
        -------
        out[0] : float
            Minimum number one, located at: ``rzero``.
        out[1] : float
            Minimum number two, located at: ``rzero + 2*width``.
        out[2] : float
            Maximum, located at: ``rzero + width``.

        """
        rzero = self.params['rzero']
        width = self.params['width']
        return rzero, rzero + 2.0 * width, rzero + width

    def potential(self, system):
        """Calculate the potential energy.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system we evaluate the potential in.

        Returns
        -------
        v_pot : float
            The potential energy.

        """
        particles = system.particles
        box = system.box
        v_pot = 0.0
        rwidth = self.params['rwidth']
        width2 = self.params['width2']
        height = self.params['height']
        for pair in particles.pairs():
            i, j, itype, jtype = pair
            if self._activate(itype, jtype):
                delta = box.pbc_dist_coordinate(particles.pos[i] -
                                                particles.pos[j])
                delr = np.sqrt(np.dot(delta, delta))
                v_pot += (height * (1.0 - (((delr - rwidth)**2) / width2))**2)
        return v_pot

    def force(self, system):
        """Calculate the force.

        We also calculate the virial here, since the force is evaluated.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system we evaluate the potential in.

        Returns
        -------
        forces : numpy.array
            The force as a numpy.array of the same shape as the
            positions in `particles.pos`.
        virial : numpy.array
            The virial, as a symmetric matrix with dimensions (dim, dim)
            where dim is given by the box.

        """
        particles = system.particles
        box = system.box
        forces = np.zeros(particles.pos.shape)
        virial = np.zeros((box.dim, box.dim))
        for pair in particles.pairs():
            i, j, itype, jtype = pair
            if self._activate(itype, jtype):
                delta = box.pbc_dist_coordinate(particles.pos[i] -
                                                particles.pos[j])
                delr = np.sqrt(np.dot(delta, delta))
                diff = delr - self.params['rwidth']
                forceij = (
                    self.params['height4'] *
                    (1.0 - diff**2 / self.params['width2']) *
                    (diff / self.params['width2'])
                )
                forceij = forceij * delta / delr
                forces[j] -= forceij
                forces[i] += forceij
                virial += np.outer(forceij, delta)
        return forces, virial

    def potential_and_force(self, system):
        """Calculate the force & potential.

        We also calculate the virial here, since the force is evaluated.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system we evaluate the potential in.

        Returns
        -------
        out[0] : float
            The potential energy as a float.
        out[1] : numpy.array
            The force as a numpy.array of the same shape as the
            positions in `particles.pos`.
        out[2] : numpy.array
            The virial, as a symmetric matrix with dimensions (dim, dim)
            where dim is given by the box.

        """
        particles = system.particles
        forces = np.zeros(particles.pos.shape)
        virial = np.zeros((system.box.dim, system.box.dim))
        v_pot = 0.0
        for pair in particles.pairs():
            i, j, itype, jtype = pair
            if self._activate(itype, jtype):
                delta = system.box.pbc_dist_coordinate(particles.pos[i] -
                                                       particles.pos[j])
                delr = np.sqrt(np.dot(delta, delta))
                diff = delr - self.params['rwidth']
                v_pot += (self.params['height'] *
                          (1.0 - diff**2 / self.params['width2'])**2)
                forceij = (self.params['height4'] *
                           (1.0 - diff**2 / self.params['width2']) *
                           (diff / self.params['width2']))
                forceij = forceij * delta / delr
                forces[i] += forceij
                forces[j] -= forceij
                virial += np.outer(forceij, delta)
        return v_pot, forces, virial
