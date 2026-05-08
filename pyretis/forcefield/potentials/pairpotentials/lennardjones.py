# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Module defining Lennard-Jones pair potentials.

This module defines the Lennard-Jones potential for PyRETIS.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PairLennardJonesCut (:py:class:`.PairLennardJonesCut`)
    A class representing a Lennard-Jones 6-12 potential implemented
    in pure python.

PairLennardJonesCutnp (:py:class:`.PairLennardJonesCutnp`)
    A class representing a Lennard-Jones 6-12 potential implemented
    using numpy.
"""
import logging
import numpy as np
from pyretis.forcefield.potential import PotentialFunction
from .pairpotential import generate_pair_interactions


logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = ['PairLennardJonesCut', 'PairLennardJonesCutnp']


class PairLennardJonesCut(PotentialFunction):
    r"""Lennard-Jones 6-12 potential in pure Python.

    This class implements as simple Lennard-Jones 6-12 potential which
    employs a simple cut-off and can be shifted. The potential energy
    (:math:`V_\text{pot}`) is defined in the usual way for an
    interacting pair of particles a distance :math:`r` apart,

    .. math::

       V_\text{pot} = 4 \varepsilon \left( x^{12} - x^{6} \right),

    where :math:`x = \sigma/r` and :math:`\varepsilon`
    and :math:`\sigma` are the potential parameters. The parameters are
    stored as attributes of the potential, we store one set for each
    kind of pair interaction. Parameters can be generated with a
    specific mixing rule by the force field.

    This implementation is in pure python (yes we are double looping!)
    and it is slow. It should not be used for production, please
    consider the numpy aware `PairLennardJonesCutnp` instead.

    Attributes
    ----------
    params : dict
        The parameters for the potential. This dict is assumed to
        contain parameters for pairs, i.e. for interactions.
    _lj : dict of dict
        Lennard-Jones parameters used for calculation of the force.
        Keys are the pairs (particle types) that may interact.

        [1] Calculated as: ``48.0 * epsilon * sigma**12``
        [2] Calculated as: ``24.0 * epsilon * sigma**6``
        [3] Calculated as: ``4.0 * epsilon * sigma**12``
        [4] Calculated as: ``4.0 * epsilon * sigma**6``
    _offset : dict
        Potential values for shifting the potential if requested.
        This is the potential evaluated at the cut-off.
    _rcut2 : dict
        The squared cut-off for each interaction type.
        Keys are the pairs (particle types) that may interact.

    """

    def __init__(self, dim=3, shift=True, mixing='geometric',
                 desc='Lennard-Jones pair potential'):
        """Initialise the Lennard-Jones potential.

        Parameters
        ----------
        dim : int, optional
            The dimensionality to use.
        shift : boolean, optional
            Determines if the potential should be shifted or not.
        mixing : string, optional
            Determines how we should mix potential parameters.
        desc : string, optional
            Description of the potential.

        """
        super().__init__(dim=dim, desc=desc)
        self.shift = shift
        self._lj = {1: {}, 2: {}, 3: {}, 4: {}}
        self._rcut2 = {}
        self._offset = {}
        self.params = {}
        self.mixing = mixing

    def set_parameters(self, parameters):
        """Update all parameters.

        Here, we generate pair interactions, since that is what this
        potential actually is using.

        Parameters
        ----------
        parameters : dict
            The input pair parameters.

        """
        self.params = {}
        pair_param = generate_pair_interactions(parameters, self.mixing)
        for pair, item in pair_param.items():
            eps_ij = item['epsilon']
            sig_ij = item['sigma']
            rcut = item['rcut']
            self._lj[1][pair] = 48.0 * eps_ij * sig_ij**12
            self._lj[2][pair] = 24.0 * eps_ij * sig_ij**6
            self._lj[3][pair] = 4.0 * eps_ij * sig_ij**12
            self._lj[4][pair] = 4.0 * eps_ij * sig_ij**6
            self._rcut2[pair] = rcut**2
            vcut = 0.0
            if self.shift:
                try:
                    vcut = 4.0 * eps_ij * ((sig_ij / rcut)**12 -
                                           (sig_ij / rcut)**6)
                except ZeroDivisionError:
                    vcut = 0.0
            self._offset[pair] = vcut
            self.params[pair] = item

    def __str__(self):
        """Generate a string with the potential parameters.

        It will generate a string with both pair and atom parameters.

        Returns
        -------
        out : string
            Table with the parameters of all interactions.

        """
        strparam = [self.desc]
        strparam += ['Potential parameters, Lennard-Jones:']
        useshift = 'yes' if self.shift else 'no'
        strparam.append(f'Shift potential: {useshift}')
        atmformat = '{0:12s} {1:>9s} {2:>9s} {3:>9s}'
        atmformat2 = '{0:12s} {1:>9.4f} {2:>9.4f} {3:>9.4f}'
        strparam.append('Pair parameters:')
        strparam.append(atmformat.format('Atom/pair', 'epsilon', 'sigma',
                                         'cut-off'))
        for pair in sorted(self.params):
            eps_ij = self.params[pair]['epsilon']
            sig_ij = self.params[pair]['sigma']
            rcut = np.sqrt(self._rcut2[pair])
            stri = f'{pair[0]}-{pair[1]}'
            strparam.append(atmformat2.format(stri, eps_ij, sig_ij, rcut))
        return '\n'.join(strparam)

    def potential(self, system):
        """Calculate the potential energy for the Lennard-Jones interaction.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system for which we calculate the potential.

        Returns
        -------
        The potential energy as a float.

        """
        particles = system.particles
        box = system.box
        v_pot = 0.0
        for pair in particles.pairs():
            i, j, itype, jtype = pair
            delta = box.pbc_dist_coordinate(particles.pos[i] -
                                            particles.pos[j])
            rsq = np.dot(delta, delta)
            if rsq < self._rcut2[itype, jtype]:
                r2inv = 1.0/rsq
                r6inv = r2inv**3
                v_pot += (r6inv * (self._lj[3][itype, jtype] * r6inv -
                                   self._lj[4][itype, jtype]) -
                          self._offset[itype, jtype])
        return v_pot

    def force(self, system):
        """Calculate the force for the Lennard-Jones interaction.

        We also calculate the virial here, since the force
        is evaluated.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system for which we calculate the force.

        Returns
        -------
        out[0] : numpy.array
            The force as a numpy.array.
        out[1] : numpy.array
            The virial as a numpy.array.

        """
        particles = system.particles
        forces = np.zeros(particles.pos.shape)
        virial = np.zeros((system.box.dim, system.box.dim))
        for pair in particles.pairs():
            i, j, itype, jtype = pair
            delta = system.box.pbc_dist_coordinate(particles.pos[i] -
                                                   particles.pos[j])
            if np.dot(delta, delta) < self._rcut2[itype, jtype]:
                r2inv = 1.0 / np.dot(delta, delta)
                r6inv = r2inv**3
                forcelj = r2inv * r6inv * (self._lj[1][itype, jtype] * r6inv -
                                           self._lj[2][itype, jtype])
                forceij = forcelj * delta
                forces[i] += forceij
                forces[j] -= forceij
                virial += np.outer(forceij, delta)
        return forces, virial

    def potential_and_force(self, system):
        """Calculate potential and force for the Lennard-Jones interaction.

        Since the force is evaluated, the virial is also calculated.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system for which we calculate the potential and force.

        Note
        ----
        Currently, the virial is only calculated for all the particles.
        It is not calculated as per atom virial. The virial
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
        v_pot = 0.0
        forces = np.zeros(system.particles.pos.shape)
        virial = np.zeros((system.box.dim, system.box.dim))
        for pair in system.particles.pairs():
            i, j, itype, jtype = pair
            delta = system.box.pbc_dist_coordinate(system.particles.pos[i] -
                                                   system.particles.pos[j])
            rsq = np.dot(delta, delta)
            if rsq < self._rcut2[itype, jtype]:
                r2inv = 1.0 / rsq
                r6inv = r2inv**3
                v_pot += (r6inv * (self._lj[3][itype, jtype] * r6inv -
                                   self._lj[4][itype, jtype]) -
                          self._offset[itype, jtype])
                forceij = (delta * r2inv * r6inv *
                           (self._lj[1][itype, jtype] * r6inv -
                            self._lj[2][itype, jtype]))
                forces[i] += forceij
                forces[j] -= forceij
                virial += np.outer(forceij, delta)
        return v_pot, forces, virial


class PairLennardJonesCutnp(PairLennardJonesCut):
    """Lennard-Jones 6-12 potential with numpy.

    A Lennard-Jones 6-12 potential with a simple cut-off which can be
    shifted. `PairLennardJonesCutnp` uses numpy for calculations, i.e.
    most operations are recast as numpy.array operations. Otherwise, it
    is similar to `PairLennardJonesCut`.
    """

    def __init__(self, dim=3, shift=True, mixing='geometric',
                 desc='Lennard-Jones pair potential (numpy)'):
        """Initialise the Lennard-Jones potential.

        Parameters
        ----------
        dim : int, optional
            The dimensionality to use.
        shift : boolean, optional
            Determines if the potential should be shifted or not.
        mixing : string, optional
            Describes the mixing rules for the parameters.
        desc : string, optional
            Description of the potential.

        """
        super().__init__(dim=dim, desc=desc,
                         shift=shift, mixing=mixing)

    def potential(self, system):
        """Calculate the potential energy for the Lennard-Jones interaction.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system for which we calculate the potential.

        Returns
        -------
        out : float
            The potential energy as a float.

        """
        particles = system.particles
        box = system.box
        pot = 0.0
        # the particle list may implement a list which we can
        # loop over. This could be some kind of fancy neighbour list
        # here, we ignore this and loop over all pairs using numpy.
        for i, particle_i in enumerate(particles.pos[:-1]):
            itype = particles.ptype[i]
            delta = particle_i - particles.pos[i+1:]
            delta = box.pbc_dist_matrix(delta)
            rsq = np.einsum('ij, ij->i', delta, delta)
            k = np.where(_check_cutoff(self._rcut2, rsq,
                                       particles.ptype[i+1:],
                                       itype))[0]
            if len(k) > 0:  # pylint: disable=len-as-condition
                r6inv = 1.0 / rsq[k]**3
                pot += np.sum(_pot_term(self._lj, self._offset,
                                        r6inv, particles.ptype[k+i+1], itype))
        return pot

    def force(self, system):
        """Calculate the force for the Lennard-Jones interaction.

        We also calculate the virial here, since the force
        is evaluated.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system for which we calculate the force.

        Note
        ----
        The way the "dim" is used may be reconsidered. There is
        already a self.dim parameter for the potential class.

        Returns
        -------
        out[0] : numpy.array
            The force as a numpy.array of the same shape as the
            positions in particles.pos.
        out[1] : numpy.array
            The virial, as a symmetric matrix with dimensions (dim, dim)
            where dim is given by the box.

        """
        particles = system.particles
        forces = np.zeros(particles.pos.shape)
        virial = np.zeros((system.box.dim, system.box.dim))
        for i, particle_i in enumerate(particles.pos[:-1]):
            itype = particles.ptype[i]
            delta = particle_i - particles.pos[i+1:]
            delta = system.box.pbc_dist_matrix(delta)
            rsq = np.einsum('ij, ij->i', delta, delta)
            k = np.where(_check_cutoff(self._rcut2, rsq,
                                       particles.ptype[i+1:],
                                       itype))[0]
            if len(k) > 0:  # pylint: disable=len-as-condition
                r2inv = 1.0 / rsq[k]
                r6inv = r2inv**3
                forcelj = _force_term(self._lj, r2inv, r6inv,
                                      particles.ptype[k+i+1], itype)
                forceij = np.einsum('i,ij->ij', forcelj, delta[k])
                forces[i] += np.sum(forceij, axis=0)
                forces[k+i+1] -= forceij
                virial += np.einsum('ij,ik->jk', forceij, delta[k])
        return forces, virial

    def potential_and_force(self, system):
        """Calculate the potential & force for the Lennard-Jones interaction.

        We also calculate the virial here, since the force is evaluated.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The system for which we calculate the potential and force.

        Note
        ----
        Currently, the virial is only calculated for all the particles.
        It is not calculated as a per atom virial. The virial per
        atom might be useful to obtain a local pressure or stress,
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
            The virial, as a symmetric matrix with dimensions (dim, dim)
            where dim is given by the box.

        """
        particles = system.particles
        pot = 0.0
        forces = np.zeros(particles.pos.shape)
        virial = np.zeros((system.box.dim, system.box.dim))
        for i, particle_i in enumerate(particles.pos[:-1]):
            delta = particle_i - particles.pos[i+1:]
            delta = system.box.pbc_dist_matrix(delta)
            rsq = np.einsum('ij, ij->i', delta, delta)
            k = np.where(_check_cutoff(self._rcut2, rsq,
                                       particles.ptype[i+1:],
                                       particles.ptype[i]))[0]
            if len(k) > 0:  # pylint: disable=len-as-condition
                r2inv = 1.0 / rsq[k]
                r6inv = r2inv**3
                pot += np.sum(_pot_term(self._lj, self._offset,
                                        r6inv, particles.ptype[k+i+1],
                                        particles.ptype[i]))
                forcelj = _force_term(self._lj, r2inv, r6inv,
                                      particles.ptype[k+i+1],
                                      particles.ptype[i])
                forceij = np.einsum('i,ij->ij', forcelj, delta[k])
                forces[i] += np.sum(forceij, axis=0)
                forces[k+i+1] -= forceij
                virial += np.einsum('ij,ik->jk', forceij, delta[k])
        return pot, forces, virial


@np.vectorize
def _pot_term(llj, offset, r6inv, itype, jtype):
    """Lennard Jones potential term."""
    return (r6inv * (llj[3][itype, jtype] * r6inv - llj[4][itype, jtype])
            - offset[itype, jtype])


@np.vectorize
def _force_term(llj, r2inv, r6inv, jtype, itype):
    """Lennard Jones force term."""
    return r2inv * r6inv * (llj[1][itype, jtype] * r6inv -
                            llj[2][itype, jtype])


@np.vectorize
def _check_cutoff(rcut2, rsq, jtype, itype):
    """Check if we are closer than the cut-off."""
    return rsq < rcut2[itype, jtype]
