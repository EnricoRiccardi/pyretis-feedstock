# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This file contains functions that act on (a selection of) particles.

These functions are intended for calculating particle properties as the
kinetic temperature, pressure etc.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

atomic_kinetic_energy_tensor (:py:func:`.atomic_kinetic_energy_tensor`)
    Return the kinetic energy tensor for each atom in a selection
    of particles.

calculate_kinetic_energy (:py:func:`.calculate_kinetic_energy`)
    Return the kinetic energy of a collection of particles.

calculate_kinetic_energy_tensor (:py:func:`.calculate_kinetic_energy_tensor`)
    Return the kinetic energy tensor for a selection of particles.

calculate_kinetic_temperature (:py:func:`calculate_kinetic_temperature`)
    Return the kinetic temperature of a collection of particles.

calculate_linear_momentum (:py:func:`calculate_linear_momentum`)
    Calculates the linear momentum of a collection of particles.

calculate_pressure_from_temp (:py:func:`calculate_pressure_from_temp`)
    Return the scalar pressure using the temperature and the virial.

calculate_pressure_tensor (:py:func:`calculate_pressure_tensor`)
    Return the pressure tensor, obtained from the virial and the kinetic
    energy tensor.

calculate_scalar_pressure (:py:func:`calculate_scalar_pressure`)
    Return the scalar pressure (from the trace of the pressure tensor).

calculate_thermo (:py:func:`calculate_thermo`)
    Calculate and return several "thermodynamic" properties as the
    potential, kinetic and total energies per particle, the temperature,
    the pressure and the momentum.

calculate_thermo_path (:py:func:`calculate_thermo_path`)
    Calculate and return some thermodynamic properties. This method
    is similar to the `calculate_thermo`, however, it is simpler and
    calculates fewer quantities.

kinetic_energy (:py:func:`kinetic_energy`)
    Return the kinetic energy for velocities and masses given as
    numpy arrays.

kinetic_temperature (:py:func:`kinetic_temperature`)
    Return the temperature for velocities and masses given as
    numpy arrays.

reset_momentum (:py:func:`reset_momentum`)
    Set linear momentum (for a selection of particles) to zero.
"""
import numpy as np

np.set_printoptions(legacy='1.21')

__all__ = [
    'atomic_kinetic_energy_tensor',
    'calculate_kinetic_energy',
    'calculate_kinetic_energy_tensor',
    'calculate_kinetic_temperature',
    'calculate_linear_momentum',
    'calculate_pressure_from_temp',
    'calculate_pressure_tensor',
    'calculate_scalar_pressure',
    'calculate_thermo',
    'calculate_thermo_path',
    'kinetic_energy',
    'kinetic_temperature',
    'reset_momentum'
]


def _get_vel_mass(particles, selection=None):
    """Return velocity and mass for a selection.

    This is just for convenience since we are using this
    selection a lot.

    Parameters
    ----------
    particles : object like :py:class:`.Particles`
        This object represents the particles.
    selection : list of integers, optional
        A list with indices of particles to use in the calculation.

    Returns
    -------
    out[0] : numpy.array
        The velocities corresponding to the selection.
    out[1] : numpy.array, optional
        The masses corresponding to the selection.

    """
    if selection is None:
        vel = particles.vel
        mass = particles.mass
    else:
        vel = particles.vel[selection]
        mass = particles.mass[selection]
    return vel, mass


def atomic_kinetic_energy_tensor(particles, selection=None):
    """Return kinetic energy tensors for a particle selection.

    Parameters
    ----------
    particles : object like :py:class:`.Particles`
        This object represents the particles.
    selection : list of integers, optional
        A list with indices of particles to use in the calculation.

    Returns
    -------
    kin : numpy.array
        A numpy array with dimensionality equal to
        (`len(selection)`, `dim`, `dim`) where `dim` is the number of
        dimensions used in the velocities. `kin[i]` contains the kinetic
        energy tensor formed by the outer product of `mol[selection][i]`
        and `vel[selection][i]`. The sum of the tensor should equal the
        output from `calculate_kinetic_energy_tensor`.

    """
    vel, mass = _get_vel_mass(particles, selection=selection)
    mom = vel * mass
    if len(mass) == 1:  # in general: selection != particles.npart
        kin = 0.5 * np.outer(mom, vel)
    else:
        kin = 0.5 * np.einsum('ij,ik->ijk', mom, vel)
    return kin


def calculate_kinetic_energy(particles, selection=None, kin_tensor=None):
    """Return the kinetic energy of a collection of particles.

    Parameters
    ----------
    particles : object like :py:class:`.Particles`
        This object represents the particles.
    selection : list of integers, optional
        A list with indices of particles to use in the calculation.
    kin_tensor : numpy.array, optional
        If kinetic_tensor is not given, the kinetic energy tensor will
        be calculated.

    Returns
    -------
    out[0] : float
        The scalar kinetic energy.
    out[1] : numpy.array
        The kinetic energy tensor.

    """
    if kin_tensor is None:
        kin_tensor = calculate_kinetic_energy_tensor(particles,
                                                     selection=selection)
    return kin_tensor.trace(), kin_tensor


def calculate_kinetic_energy_tensor(particles, selection=None):
    """Return the kinetic energy tensor for a selection of particles.

    The tensor is formed as the outer product of the velocities.

    Parameters
    ----------
    particles : object like :py:class:`.Particles`
        This object represents the particles.
    selection : list of integers, optional
        A list with indices of particles to use in the calculation.

    Returns
    -------
    out : numpy.array
        The kinetic energy tensor. Dimensionality is equal to (dim, dim)
        where dim is the number of dimensions used in the velocities.
        The trace gives the kinetic energy.

    """
    vel, mass = _get_vel_mass(particles, selection=selection)
    _, kin = kinetic_energy(vel, mass)
    return kin


def kinetic_energy(vel, mass):
    """Obtain the kinetic energy for given velocities and masses.

    Parameters
    ----------
    vel : numpy.array
        The velocities
    mass : numpy.array
        The masses. This is assumed to be a column vector.

    Returns
    -------
    out[0] : float
        The kinetic energy
    out[1] : numpy.array
        The kinetic energy tensor.

    """
    mom = vel * mass
    if len(mass) == 1:
        kin = 0.5 * np.outer(mom, vel)
    else:
        kin = 0.5 * np.einsum('ij,ik->jk', mom, vel)
    return kin.trace(), kin


def calculate_kinetic_temperature(particles, boltzmann, dof=None,
                                  selection=None,
                                  kin_tensor=None):
    """Return the kinetic temperature of a collection of particles.

    Parameters
    ----------
    particles : object like :py:class:`.Particles`
        This object represents the particles.
    boltzmann : float
        This is the Boltzmann factor/constant in correct units.
    dof : list of floats, optional
        The degrees of freedom to subtract. Its shape should
        be equal to the number of dimensions.
    selection : list of integers, optional
        A list with indices of particles to use in the calculation.
    kin_tensor : numpy.array optional
        The kinetic energy tensor. If the kinetic energy tensor is not
        given, it will be recalculated here.

    Returns
    -------
    out[0] : numpy.array
        Array with the same size as the kinetic energy. It
        contains the temperature in each spatial dimension.
    out[1] : float
        The temperature averaged over all dimensions.
    out[2] : numpy.array
        The kinetic energy tensor.

    """
    vel, mass = _get_vel_mass(particles, selection=selection)
    npart = len(mass)  # using mass, since selection may be != particles.npart
    ndof = npart * np.ones(vel[0].shape)

    if kin_tensor is None:
        kin_tensor = calculate_kinetic_energy_tensor(particles,
                                                     selection=selection)
    if dof is not None:
        ndof = ndof - dof
    temperature = (2.0 * kin_tensor.diagonal() / ndof) / boltzmann
    return temperature, np.average(temperature), kin_tensor


def kinetic_temperature(vel, mass, boltzmann, dof=None, kin_tensor=None):
    """Return the kinetic temperature given velocities and masses.

    This method does not work on a particle object, but rather with
    numpy arrays. That is, it is intended for use when we can't rely
    on the particle object.

    Parameters
    ----------
    vel : numpy.array
        The velocities
    mass : numpy.array
        The masses. This is assumed to be a column vector.
    boltzmann : float
        This is the Boltzmann factor/constant in correct units.
    dof : list of floats, optional
        The degrees of freedom to subtract. Its shape should
        be equal to the number of dimensions.
    kin_tensor : numpy.array optional
        The kinetic energy tensor. If the kinetic energy tensor is not
        given, it will be recalculated here.

    Returns
    -------
    out[0] : numpy.array
        It contains the temperature in each spatial dimension.
        Array with the same size as the kinetic energy.
    out[1] : float
        The temperature averaged over all dimensions.
    out[2] : numpy.array
        The kinetic energy tensor.

    """
    npart = len(mass)  # using mass, since selection may be != particles.npart
    ndof = npart * np.ones(vel[0].shape)

    if kin_tensor is None:
        _, kin_tensor = kinetic_energy(vel, mass)

    if dof is not None:
        ndof = ndof - dof
    temperature = (2.0 * kin_tensor.diagonal() / ndof) / boltzmann
    return temperature, np.average(temperature), kin_tensor


def calculate_linear_momentum(particles, selection=None):
    """Calculate the linear momentum for a collection of particles.

    Parameters
    ----------
    particles : object like :py:class:`.Particles`
        This object represents the particles.
    selection : list of integers, optional
        A list with indices of particles to use in the calculation.

    Returns
    -------
    out : numpy.array
        The array contains the linear momentum for each dimension.

    """
    vel, mass = _get_vel_mass(particles, selection=selection)
    return np.sum(vel * mass, axis=0)


def calculate_pressure_from_temp(particles, dim, boltzmann, volume,
                                 dof=None):
    """Evaluate the scalar pressure.

    The scalar pressure is here calculated  using the temperature
    and the degrees of freedom.

    Parameters
    ----------
    particles : object like :py:class:`.Particles`
        This object represents the particles.
    dim : int
        This is the dimensionality of the system.
        Typically provided by `system.get_dim()`.
    boltzmann : float
        This is the Boltzmann factor/constant in correct units.
        Typically it can be supplied by `system.get_boltzmann()`.
    volume : float
        This is the volume 'occupied' by the particles. It can typically
        be obtained by a `box.calculate_volume()`
    dof : list of floats, optional
        The degrees of freedom to subtract. Its shape should
        be equal to the number of dimensions.

    Returns
    -------
    out[0] : float
        Pressure times volume.
    out[1] : float
        The pressure.

    Note
    ----
    This function may possibly be removed - it does not appear to be
    very useful right now.

    """
    if dof is None:
        ndof = particles.npart
    else:
        ndof = (particles.npart * dim - np.sum(dof)) / float(dim)
    _, temperature, _ = calculate_kinetic_temperature(particles, boltzmann,
                                                      dof=dof)
    pressvolume = ndof * temperature * boltzmann
    pressvolume += particles.virial.trace() / float(dim)
    press = pressvolume / volume
    return pressvolume, press


def calculate_pressure_tensor(particles, volume, kin_tensor=None):
    """Calculate the pressure tensor.

    The pressure tensor is obtained from the virial the kinetic
    energy tensor.

    Parameters
    ----------
    particles : object like :py:class:`.Particles`
        This object represents the particles.
    volume : float
        This is the volume 'occupied' by the particles. It can typically
        be obtained by a `box.calculate_volume()`.
    kin_tensor : numpy.array, optional
        The kinetic energy tensor. If `kin_tensor` is not given, it will
        be calculated here.

    Returns
    -------
    out : numpy.array
        The symmetric pressure tensor, dimensions (`dim`, `dim`), where
        `dim` = the number of dimensions considered in the simulation.

    """
    if kin_tensor is None:
        kin_tensor = calculate_kinetic_energy_tensor(particles, selection=None)
    pressure = (particles.virial + 2. * kin_tensor) / volume
    return pressure


def calculate_scalar_pressure(particles, volume, dim, press_tensor=None,
                              kin_tensor=None):
    """Evaluate the scalar pressure using the pressure tensor.

    Parameters
    ----------
    particles : object like :py:class:`.Particles`
        This object represents the particles.
    volume : float
        This is the volume 'occupied' by the particles. It can typically
        be obtained by a `box.calculate_volume()`.
    dim : int
        This is the dimensionality of the system. Typically provided by
        `system.get_dim()`
    press_tensor : numpy.array, optional
        If `press_tensor` is not given, the pressure tensor will be
        calculated here.
    kin_tensor : numpy.array, optional
        If `kin_tensor` is not given, the kinetic energy tensor will be
        calculated here.

    Returns
    -------
    out : float
        The scalar pressure averaged over the diagonal components of
        the pressure tensor.

    """
    if press_tensor is None:
        press_tensor = calculate_pressure_tensor(particles, volume,
                                                 kin_tensor=kin_tensor)
    return press_tensor.trace() / float(dim)


def calculate_thermo(system, dof=None, dim=None, volume=None, vpot=None):
    """Calculate and return several thermodynamic properties.

    The calculated properties are the potential, kinetic and total
    energies per particle, the temperature, the pressure and the
    momentum.

    Parameters
    ----------
    system : object like :py:class:`.System`
        This object is used to access the particles and the box.
    dof : list of floats, optional
        The degrees of freedom.
    dim : float, optional
        The dimensionality of, typically provided with a
        `system.get_dim()`.
    volume : float, optional
        This is the volume 'occupied' by the particles. It can typically
        be obtained by a `system.box.calculate_volume()`.
    vpot : float, optional
        The potential energy of the particles.

    Returns
    -------
    out : dict
        This dict contains the float that is calculated in this routine.

    """
    if volume is None:
        volume = system.box.calculate_volume()
    if vpot is None:
        vpot = system.particles.vpot
    if dim is None:
        dim = system.get_dim()
    if dof is None:
        dof = system.temperature['dof']
    kin_tens = calculate_kinetic_energy_tensor(system.particles)
    _, temp, _ = calculate_kinetic_temperature(system.particles,
                                               system.get_boltzmann(),
                                               dof=dof,
                                               kin_tensor=kin_tens)
    press_tens = calculate_pressure_tensor(system.particles, volume,
                                           kin_tensor=kin_tens)
    ekin = kin_tens.trace()
    press = calculate_scalar_pressure(system.particles, volume, dim,
                                      kin_tensor=kin_tens)
    mom = calculate_linear_momentum(system.particles)
    npart = float(system.particles.npart)
    result = {'vpot': vpot / npart, 'ekin': ekin / npart,
              'etot': (ekin + vpot) / npart,
              'temp': temp, 'press': press, 'mom': mom,
              'press-tens': press_tens}
    return result


def calculate_thermo_path(system):
    """Calculate and return several thermodynamic properties.

    The calculated properties are the potential, kinetic and total
    energies for the system and the current temperature. The name here
    ``calculate_thermo_path`` just indicates that this function is
    useful in connection with path sampling simulations, i.e. it
    just calculates a few energy terms.

    Parameters
    ----------
    system : object like :py:class:`.System`
        This object is used to access the particles and the box.

    Returns
    -------
    out : dict
        This dict contains the float that is calculated in this routine.

    """
    particles = system.particles
    kin_tens = calculate_kinetic_energy_tensor(particles)
    _, temp, _ = calculate_kinetic_temperature(particles,
                                               system.get_boltzmann(),
                                               dof=system.temperature['dof'],
                                               kin_tensor=kin_tens)
    ekin = kin_tens.trace()
    vpot = particles.vpot
    return {'vpot': vpot, 'ekin': ekin, 'etot': ekin + vpot, 'temp': temp}


def reset_momentum(particles, selection=None, dim=None):
    """Set the linear momentum of a selection of particles to zero.

    Parameters
    ----------
    particles : object like :py:class:`.Particles`
        This object represents the particles.
    selection : list of integers, optional
        A list with indices of particles to use in the calculation.
    dim : list or None, optional
        If ``dim`` is None, the momentum will be reset for ALL
        dimensions. Otherwise, it will only be applied to the
        dimensions where ``dim`` is True.

    Returns
    -------
    out : None
        Returns `None` and modifies velocities of the selected
        particles.

    """
    vel, mass = _get_vel_mass(particles, selection=selection)
    mom = np.sum(vel * mass, axis=0)
    if dim is not None:
        for i, reset in enumerate(dim):
            if not reset:
                mom[i] = 0
    particles.vel[selection] -= (mom / mass.sum())
