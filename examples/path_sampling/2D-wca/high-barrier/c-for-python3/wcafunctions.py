# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Wrapper for the C extensions for the WCA potential example."""
import sys
import os
import logging
import numpy as np
from pyretis.forcefield import PotentialFunction
from pyretis.orderparameter import OrderParameter
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())
# Just to handle imports of the library:
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
try:
    import wcaforces
except ImportError:
    MSG = ('Could not import external C library for forces.'
           '\nPlease compile with: "python setup.py build_ext --inplace"')
    logger.critical(MSG)
    raise ImportError(MSG)
try:
    import wcalambda
except ImportError:
    MSG = ('Could not import external C library for order parameter.'
           '\nPlease compile with "python setup.py build_ext --inplace"!')
    logger.critical(MSG)
    raise ImportError(MSG)


__all__ = ['WCAPotential', 'WCAOrderParameter']


class WCAPotential(PotentialFunction):
    r"""class WCAPotential(PotentialFunction).

    This class implements a WCA double well + WCA fluid potential.

    Attributes
    ----------
    params : dict
        The parameters for the potential. This dict is assumed to
        contain parameters for the potential:

        * `lj1` : float
          Lennard-Jones parameters used for calculation of the force.
          Calculated as: ``48.0 * epsilon * sigma**12``
        * `lj2` : float
          Lennard-Jones parameters used for calculation of the force.
          Calculated as: ``24.0 * epsilon * sigma**6``
        * `lj3` : float
          Lennard-Jones parameters used for calculation of the potential.
          Calculated as: ``4.0 * epsilon * sigma**12``
        * `lj4` : float
          Lennard-Jones parameters used for calculation of the potential.
          Calculated as: ``4.0 * epsilon * sigma**6``
        * `offset` : float
          Potential values for shifting the potential if requested.
          This is the potential evaluated at the cutoff.
        * `rcut2` : numpy.array
          The squared cut-off for each interaction type.

    """

    def __init__(self, dim=2, shift=True,
                 desc='WCA pair potential (C)'):
        """Set up the WCA potential.

        Parameters
        ----------
        dim : int
            The dimensionality to use.
        shift : boolean
            Determines if the potential should be shifted or not

        """
        super().__init__(dim=dim, desc=desc)
        self.shift = shift
        self.params = {'height': 0.0,
                       'height4': 0.0,
                       'rwidth': 0.0,
                       'rzero': 0.0,
                       'idxi': 0,
                       'idxj': 0,
                       'width': 0.0,
                       'width2': 0.0,
                       'epsilon': 0.0,
                       'sigma': 0.0,
                       'rcut': 0.0,
                       'lj1': 0.0, 'lj2': 0.0, 'lj3': 0.0, 'lj4': 0.0,
                       'rcut2': 0.0, 'offset': 0.0}

    def set_parameters(self, parameters):
        """Update all parameters.

        Here, we generate pair interactions, since that is what this
        potential actually is using.

        Parameters
        ----------
        parameters : dict
            The input base parameters.

        """
        for key in parameters:
            if key in self.params:
                self.params[key] = parameters[key]
            else:
                msg = f'Unknown parameter {key} - ignored!'
                logging.warning(msg)
        self.params['width2'] = self.params['width']**2
        self.params['rwidth'] = self.params['rzero'] + self.params['width']
        self.params['height4'] = 4.0 * self.params['height']
        sig_ij = parameters['sigma']
        eps_ij = parameters['epsilon']
        rcut = parameters['rcut']
        self.params['lj1'] = 48.0 * eps_ij * sig_ij**12
        self.params['lj2'] = 24.0 * eps_ij * sig_ij**6
        self.params['lj3'] = 4.0 * eps_ij * sig_ij**12
        self.params['lj4'] = 4.0 * eps_ij * sig_ij**6
        self.params['rcut2'] = rcut**2
        vcut = 0.0
        if self.shift:
            try:
                vcut = 4.0 * eps_ij * ((sig_ij / rcut)**12 -
                                       (sig_ij / rcut)**6)
            except ZeroDivisionError:
                vcut = 0.0
        self.params['offset'] = vcut

    def potential(self, system):
        """Calculate the potential energy.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system we are operating on.

        Returns
        -------
        out : float
            The potential energy.

        """
        particles = system.particles
        box = system.box
        v_pot = wcaforces.potential2D(particles.pos,
                                      box.length,
                                      box.ilength,
                                      self.params['rwidth'],
                                      self.params['width2'],
                                      self.params['height'],
                                      self.params['idxi'],
                                      self.params['idxj'],
                                      self.params['lj3'],
                                      self.params['lj4'],
                                      self.params['offset'],
                                      self.params['rcut2'],
                                      particles.npart)
        return v_pot

    def force(self, system):
        """Calculate the forces.

        We also calculate the virial here, since the force
        is evaluated.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system we are operating on.

        Returns
        -------
        out[0] : numpy.array
            The forces.
        out[1] : numpy.array
            The virial.

        """
        particles = system.particles
        box = system.box
        forces = np.zeros_like(particles.pos)
        virial = np.zeros((box.dim, box.dim))
        wcaforces.force2D(particles.pos,
                          box.length,
                          box.ilength,
                          forces,
                          virial,
                          self.params['rwidth'],
                          self.params['width2'],
                          self.params['height4'],
                          self.params['idxi'],
                          self.params['idxj'],
                          self.params['lj1'],
                          self.params['lj2'],
                          self.params['rcut2'],
                          particles.npart, box.dim)
        return forces, virial

    def potential_and_force(self, system):
        """Calculate potential and forces.

        Since the force is evaluated, the virial is also calculated.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system we are operating on.

        Note
        ----
        Currently, the virial is only calculated for all the particles.
        It is not calculated as a virial per atom. The virial
        per atom might be useful to obtain a local pressure or stress,
        however, this needs some consideration. Perhaps it's best to
        fully implement this as a method of planes or something similar.

        Returns
        -------
        out[0] : float
            The potential energy as a float.
        out[1] : numpy.array
            The force as a numpy.array of the same shape as the
            positions in `particles.pos`.
        out[2] : numpy.array
            The virial, as a symmetric matrix with dimensions
            (dim, dim) where dim is given by the box/system dimensions.

        """
        particles = system.particles
        box = system.box
        forces = np.zeros_like(particles.pos)
        virial = np.zeros((box.dim, box.dim))
        vpot = wcaforces.force_and_pot2D(particles.pos,
                                         box.length,
                                         box.ilength,
                                         forces,
                                         virial,
                                         self.params['rwidth'],
                                         self.params['width2'],
                                         self.params['height'],
                                         self.params['height4'],
                                         self.params['idxi'],
                                         self.params['idxj'],
                                         self.params['lj1'],
                                         self.params['lj2'],
                                         self.params['lj3'],
                                         self.params['lj4'],
                                         self.params['offset'],
                                         self.params['rcut2'],
                                         particles.npart, box.dim)
        return vpot, forces, virial


class WCAOrderParameter(OrderParameter):
    """WCAOrderParameter(OrderParamete).

    This class represents the order parameter for the WCA example.
    This is the order parameter for the high barrier case and it's
    simply the bond length.

    Attributes
    ----------
    index : tuple of ints
        The index for the particles to use.

    """

    def __init__(self, index):
        """Set up the order parameter.

        Parameters
        ----------
        index : tuple of ints
            The index for the particles to use.

        """
        super().__init__(description='WCA order parameter')
        self.index = index

    def calculate(self, system):
        """Calculate the order parameter and return it.

        Parameters
        ----------
        system : object like :py:class:`.System`
            This object is used for the actual calculation, typically
            only `system.particles.pos` and/or `system.particles.vel`
            will be used. In some cases system.forcefield can also be
            used to include specific energies for the order parameter.

        Returns
        -------
        out : float
            The order parameter.

        """
        particles = system.particles
        lamb = wcalambda.orderp(particles.pos,
                                system.box.length,
                                system.box.ilength,
                                self.index[1], self.index[0])
        lambv = wcalambda.orderv(particles.pos,
                                 particles.vel,
                                 system.box.length,
                                 system.box.ilength,
                                 self.index[1], self.index[0])
        return [lamb, lambv]


class WCAOrderParameterp(OrderParameter):
    """WCAOrderParameterp(OrderParameter).

    This class represents the order parameter for the WCA example.
    This is the order parameter for the high barrier case and it's
    simply the bond length. This is just a pure python implementation.

    Attributes
    ----------
    index : tuple of ints
        The index for the particles to use.

    """

    def __init__(self, index):
        """Set up the order parameter.

        Parameters
        ----------
        index : tuple of ints
            The index for the particles to use.

        """
        super().__init__(description='WCA order parameter')
        self.index = index

    def calculate(self, system):
        """Calculate the order parameter and return it.

        Parameters
        ----------
        system : object like :py:class:`.System`
            This object is used for the actual calculation, typically
            only `system.particles.pos` and/or `system.particles.vel`
            will be used. In some cases system.forcefield can also be
            used to include specific energies for the order parameter.

        Returns
        -------
        out : float
            The order parameter.

        """
        particles = system.particles
        delta = system.box.pbc_dist_coordinate(particles.pos[self.index[1]] -
                                               particles.pos[self.index[0]])
        lamb = np.sqrt(np.dot(delta, delta))
        delta_v = particles.vel[self.index[1]] - particles.vel[self.index[0]]
        lamb_v = np.dot(delta, delta_v) / lamb
        return [lamb, lamb_v]
