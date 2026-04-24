# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This module handles random number generation.

It derives most of the random number procedures from `RandomState` in
`numpy.random` and defines a class which used `RandomState` to generate
pseudo-random numbers.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RandomGeneratorBase (:py:class:`.RandomGeneratorBase`)
    An abstract base class for random number generators.

RandomGenerator (:py:class:`.RandomGenerator`)
    A class representing a random number generator.

ReservoirSampler (:py:class:`.ReservoirSampler`)
    A class for reservoir sampling.

MockRandomGenerator (:py:class:`.MockRandomGenerator`)
    A non-random number generator which is useful for testing.
"""
from abc import ABCMeta, abstractmethod
import logging
import numpy as np
from numpy.random import RandomState
from pyretis.core.particlefunctions import (calculate_kinetic_temperature,
                                            reset_momentum)
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = [
    'RandomGeneratorBase',
    'RandomGenerator',
    'ReservoirSampler',
    'MockRandomGenerator',
    'create_random_generator'
]


class RandomGeneratorBase(metaclass=ABCMeta):
    """A base class for random number generators.

    This is a base class for random number generators. It does not
    actually implement a generator.

    Attributes
    ----------
    seed : int
        A seed for the generator

    """

    def __init__(self, seed=0):
        """Initialise the random number generator.

        Parameters
        ----------
        seed : int, optional
            An integer used for seeding the generator if needed.

        """
        self.seed = seed

    @abstractmethod
    def rand(self, shape=1):
        """Draw random numbers in [0, 1).

        Parameters
        ----------
        shape : int, optional
            The number of numbers to draw

        Returns
        -------
        out : float
            Pseudo-random number in [0, 1)

        """
        return

    @abstractmethod
    def get_state(self):
        """Return info about the current state."""
        return

    @abstractmethod
    def set_state(self, state):
        """Set state for random generator."""
        return

    @abstractmethod
    def random_integers(self, low, high):
        """Draw random integers in [low, high].

        Parameters
        ----------
        low : int
            This is the lower limit
        high : int
            This is the upper limit

        Returns
        -------
        out : int
            The pseudo-random integers in [low, high].

        """
        return

    @abstractmethod
    def normal(self, loc=0.0, scale=1.0, size=None):
        """Return values from a normal distribution.

        Parameters
        ----------
        loc : float, optional
            The mean of the distribution
        scale : float, optional
            The standard deviation of the distribution
        size : int, tuple of ints, optional
            Output shape, i.e. how many values to generate. Default is
            None which is just a single value.

        Returns
        -------
        out : float, numpy.array of floats
            The random numbers generated.

        """
        return

    @abstractmethod
    def multivariate_normal(self, mean, cov, cho=None, size=1):
        """Draw numbers from a multi-variate distribution.

        Parameters
        ----------
        mean : numpy array (1D, 2)
            Mean of the N-dimensional array
        cov : numpy array (2D, (2, 2))
            Covariance matrix of the distribution.
        cho : numpy.array (2D, (2, 2)), optional
            Cholesky factorization of the covariance. If not given,
            it will be calculated here.
        size : int, optional
            The number of samples to do.

        Returns
        -------
        out : float or numpy.array of floats size
            The random numbers that are drawn here.

        """
        return

    def generate_maxwellian_velocities(self, particles, boltzmann, temperature,
                                       dof, selection=None, momentum=True):
        """Generate velocities from a Maxwell distribution.

        The velocities are drawn to match a given temperature. This
        function can be applied to a subset of the particles.

        The generation is performed in three steps:

        1.  Initial velocities are drawn from a standard normal
            distribution (0, 1).
        2.  Velocities are scaled by ``1.0 / sqrt(mass)`` to account
            for particle masses. The linear momentum is reset if
            requested.
        3.  The velocities are scaled to match the target kinetic
            temperature exactly.

        Parameters
        ----------
        particles : object
            These are the particles to set the velocity of.
        boltzmann : float
            The Boltzmann factor in correct units.
        temperature : float
            The desired temperature. Typically,
            ``system.temperature['set']`` will be used here.
        dof : list of floats or numpy.array
            The degrees of freedom to subtract. Its shape should be
            equal to the number of dimensions.
        selection : list of ints, optional
            A list with indices of the particles to consider. Can be
            used to only apply it to a selection of particles.
        momentum : boolean, optional
            If True, we will reset the linear momentum of the particles
            (or the selected subset) after generating the initial
            velocities.

        Returns
        -------
        None
            The velocities of the selected particles are modified
            in-place.

        See also
        --------
        draw_maxwellian_velocities : For drawing velocities without
            assigning them in-place.

        """
        if selection is None:
            vel, imass = particles.vel, particles.imass
        else:
            vel, imass = particles.vel[selection], particles.imass[selection]
        vel = np.sqrt(imass) * self.normal(loc=0.0, scale=1.0,
                                           size=vel.shape)
        # NOTE: x[None] = x for a numpy.array - this is not valid for a list.
        particles.vel[selection] = vel
        if momentum:
            reset_momentum(particles, selection=selection)

        _, avgtemp, _ = calculate_kinetic_temperature(particles, boltzmann,
                                                      dof=dof,
                                                      selection=selection)
        scale_factor = np.sqrt(temperature/avgtemp)
        particles.vel[selection] *= scale_factor

    def draw_maxwellian_velocities(self, system, sigma_v=None):
        """Draw numbers from a Gaussian distribution.

        Parameters
        ----------
        system : object like :py:class:`.System`
            This is used to determine the temperature parameter(s) and
            the shape (number of particles and dimensionality)
        sigma_v : numpy.array, optional
            The standard deviation in velocity, one for each particle.
            If it's not given it will be estimated.

        """
        if not sigma_v or sigma_v < 0.0:
            kbt = 1.0/system.temperature['beta']
            # sigma_v is (n, 1) matrix
            sigma_v = np.sqrt(kbt*system.particles.imass)
        npart, dim = system.particles.vel.shape
        vel = self.normal(loc=0.0, scale=sigma_v, size=(npart, dim))
        return vel, sigma_v


class RandomGenerator(RandomGeneratorBase):
    """A random number generator from numpy.

    This class that defines a random number generator. It will use
    `numpy.random.RandomState` for the actual generation and we refer
    to the numpy documentation [1]_.

    Attributes
    ----------
    seed : int
        A seed for the pseudo-random generator.
    rgen : object like :py:class:`numpy.random.RandomState`
        This is a container for the Mersenne Twister pseudo-random
        number generator as implemented in numpy [#]_.

    References
    ----------
    .. [#] The NumPy documentation on RandomState,
       http://docs.scipy.org/doc/numpy/reference/generated/numpy.random.RandomState.html

    """

    def __init__(self, seed=0):
        """Initialise the random number generator.

        If a seed is given, the random number generator will be seeded.

        Parameters
        ----------
        seed : int, optional
            An integer used for seeding the generator if needed.

        """
        super().__init__(seed=seed)
        self.rgen = RandomState(seed=seed)

    def choice(self, a, size=None, replace=True, p=None):
        """Chooses random samples.

        Parameters
        ----------
        a : iterable, int
            The group to choose from, an int is converted to range(a).
        size : int, optional
            The number of samples to choose.
        replace : bool, optional
            If to replace samples between draws, default True.
        p : iterable, optional
            The probabilities for every option in a,
            default is an uniform dirstribution.

        Returns
        -------
        choice : array-like
            The picked choices.

        """
        return self.rgen.choice(a, size, replace, p)

    def rand(self, shape=1):
        """Draw random numbers in [0, 1).

        Parameters
        ----------
        shape : int or tuple of ints, optional
            The number of numbers to draw. Default is 1, which
            returns an array of one random number.

        Returns
        -------
        out : float or numpy.array
            Pseudo-random number in [0, 1). If shape is an integer,
            a 1D array is returned. If shape is a tuple, an array
            with the given shape is returned.

        Note
        ----
        This is a convenience wrapper for `self.rgen.random_sample`.

        """
        if shape is None:
            shape = 1
        return self.rgen.rand(shape)

    def get_state(self):
        """Return current state."""
        state = {'seed': self.seed,
                 'state': self.rgen.get_state(),
                 'rgen': 'rgen'}
        return state

    def set_state(self, state):
        """Set state for random generator."""
        return self.rgen.set_state(state['state'])

    def random_integers(self, low, high):
        """Draw random integers in [low, high].

        Parameters
        ----------
        low : int
            This is the lower limit
        high : int
            This is the upper limit

        Returns
        -------
        out : int
            The pseudo-random integers in [low, high].

        Note
        ----
        np.random.randint(low, high) is defined as drawing
        from `low` (inclusive) to `high` (exclusive).

        """
        return self.rgen.randint(low, high + 1)

    def normal(self, loc=0.0, scale=1.0, size=None):
        """Return values from a normal distribution.

        Parameters
        ----------
        loc : float, optional
            The mean of the distribution
        scale : float, optional
            The standard deviation of the distribution
        size : int, tuple of ints, optional
            Output shape, i.e. how many values to generate. Default is
            None which is just a single value.

        Returns
        -------
        out : float, numpy.array of floats
            The random numbers generated.

        """
        return self.rgen.normal(loc=loc, scale=scale, size=size)

    def multivariate_normal(self, mean, cov, cho=None, size=1):
        """Draw numbers from a multi-variate distribution.

        This is an attempt on speeding up the call of
        `RandomState.multivariate_normal` if we need to call it over and
        over again. Such repeated calling will do an SVD repeatedly,
        which is wasteful. In this function, this transform can be
        supplied and it is only estimated if it's not explicitly given.

        Parameters
        ----------
        mean : numpy array (1D, 2)
            Mean of the N-dimensional array
        cov : numpy array (2D, (2, 2))
            Covariance matrix of the distribution.
        cho : numpy.array (2D, (2, 2)), optional
            Cholesky factorization of the covariance. If not given,
            it will be calculated here.
        size : int, optional
            The number of samples to do.

        Returns
        -------
        out : float or numpy.array of floats size
            The random numbers that are drawn here.

        See also
        --------
        numpy.random.multivariate_normal

        """
        if cho is None:
            cho = np.linalg.cholesky(cov)
        norm = self.normal(loc=0.0, scale=1.0, size=2*size)
        norm = norm.reshape(size, 2)
        meanm = np.array([mean, ] * size)
        return meanm + np.dot(norm, cho.T)


class ReservoirSampler:
    """A class representing a reservoir sampler.

    The reservoir sampler will maintain a list of `k` items drawn
    randomly from a set of `N > k` items. The list is created and
    maintained so that we only need to store `k`items This is useful
    when `N` is very large or when storing all `N` items require a lot
    of memory. The algorithm is described by Knuth [#]_ but here we do
    a variation, so that each item may be picked several times.


    Attributes
    ----------
    rgen : object like :py:class:`numpy.random.RandomState`
        This is a container for the Mersenne Twister pseudo-random
        number generator as implemented in numpy, see the documentation
        of `RandomGenerator`.
    items : integer
        The number of items seen so far, i.e. the current `N`.
    reservoir : list
        The items we have stored.
    length : integer
        The maximum number of items to store in the reservoir
        (i.e. `k`).
    ret_idx : integer
        This is the index of the item to return if we are requesting
        items from the reservoir.

    References
    ----------
    .. [#] The Art of Computer Programming.

    """

    def __init__(self, seed=0, length=10, rgen=None):
        """Initialise the reservoir.

        Parameters
        ----------
        seed : int, optional
            An integer used for seeding the generator.
        length : int, optional
            The maximum number of items to store.
        rgen : object like :py:class:`.RandomGenerator`, optional
            In case we want to re-use a random generator object.
            If this is specified, the parameter `seed` is ignored.

        """
        if rgen is not None:
            self.rgen = rgen
        else:
            self.rgen = RandomState(seed=seed)
        self.items = 0
        self.reservoir = []
        self.length = length
        self.ret_idx = 0

    def append(self, new_item):
        """Try to add an item to the reservoir.

        Parameters
        ----------
        new_item : any type
            This is the item we try to add to the reservoir.

        """
        self.items += 1
        if self.items == 1:
            self.reservoir = [new_item for _ in range(self.length)]
        else:
            factor = 1.0/float(self.items)
            for i in range(self.length):
                if self.rgen.rand() < factor:
                    self.reservoir[i] = new_item

    def get_item(self):
        """Return the next item from the reservoir.

        The items are returned in the order they were added to the
        reservoir (with wrapping). If the end of the reservoir is
        reached, a critical message is logged and the index is reset
        to zero.

        Returns
        -------
        out : any
            The next item from the reservoir. The type depends on what
            was stored in the reservoir.

        """
        if self.ret_idx >= self.length:
            self.ret_idx = 0
            msg = ['Out of bounds in the reservoir sampler!',
                   'Please increase the size of the reservoir.']
            msgtxt = '\n'.join(msg)
            logger.critical(msgtxt)
        ret = self.reservoir[self.ret_idx]
        self.ret_idx += 1
        return ret


class MockRandomGenerator(RandomGeneratorBase):
    """A **mock** random generator, useful **only for testing**.

    This class represents a random generator that can be used for
    testing algorithms. It will simply return numbers from a
    small list of pseudo-random numbers. This class is only useful
    for testing algorithms on different systems. It should *NEVER*
    be used for actual production runs!
    """

    def __init__(self, seed=0, norm_shift=False):
        """Initialise the mock random number generator.

        Here, we set up some predefined random number which we will
        use as a pool for the generation.

        Parameters
        ----------
        seed : int, optional
            An integer used for seeding the generator if needed.
        norm_shift: boolean
            If True is will ensure that the fake 'normal distribution'
            is shifted to get the right mean.

        """
        super().__init__(seed=seed)
        self.rgen = [0.78008018, 0.04459916, 0.76596775, 0.97676713,
                     0.53799598, 0.98657116, 0.36343553, 0.55356511,
                     0.03172585, 0.48984682, 0.73416687, 0.98453452,
                     0.55129902, 0.40598753, 0.59448394, 0.26823255,
                     0.31168372, 0.05072849, 0.44876368, 0.94301709]
        self.length = len(self.rgen)
        self.randint = self.seed
        self.norm_shift = norm_shift
        logger.critical('You are using a "mock" random generator!\n')
        if norm_shift:
            logger.critical('Fake-normal is shifted.\n')
            logger.critical('Comparison with TISMOL might fail\n')
        else:
            logger.critical('Fake-normal is not shifted.\n')
            logger.critical('Random numbers not centered around 0.\n')

    def rand(self, shape=1):
        """Draw random numbers in [0, 1).

        Parameters
        ----------
        shape : int, optional
            The number of numbers to draw. Default is 1.

        Returns
        -------
        out : numpy.array
            The pseudo-random numbers in [0, 1).

        Note
        ----
        This method draws numbers from a pre-defined sequence
        (starting with 0.78008018) to ensure deterministic behavior
        during testing. If the end of the sequence is reached, it
        wraps around to the beginning.

        """
        numbers = []
        for _ in range(shape):
            if self.seed >= self.length:
                self.seed = 0
            numbers.append(self.rgen[self.seed])
            self.seed += 1
        return np.array(numbers)

    def get_state(self):
        """Return current state."""
        state = {'seed': self.seed,
                 'state': self.seed,
                 'rgen': 'mock'}
        return state

    def set_state(self, state):
        """Set current state."""
        self.seed = state['state']

    def random_integers(self, low, high):
        """Return random integers in [low, high].

        Parameters
        ----------
        low : int
            The lower limit (inclusive).
        high : int
            The upper limit (inclusive).

        Returns
        -------
        out : int
            A pseudo-random integer in [low, high].

        Note
        ----
        This method uses a number from the pre-defined sequence and
        scales it to the requested range using:
        ``idx = rand() * (high - low + 1)``.

        """
        idx = self.rand()[0]*(high-low+1)
        return int(idx.item()) + low

    def normal(self, loc=0.0, scale=1.0, size=None):
        """Return values from a "normal" distribution.

        Parameters
        ----------
        loc : float, optional
            The mean of the distribution.
        scale : float, optional
            The standard deviation of the distribution.
        size : int or tuple of ints, optional
            Output shape. Default is None, which returns a single
            value.

        Returns
        -------
        out : float or numpy.array
            The generated "random" numbers.

        Note
        ----
        This is part of a mock generator and does NOT provide a true
        normal distribution. It simply returns values from the
        pre-defined sequence, shifted by ``loc - 0.5`` if
        ``norm_shift`` is enabled.

        """
        if self.norm_shift:
            shift = loc-0.5
        else:
            shift = 0.
        if size is None:
            return self.rand(shape=1)+shift
        numbers = np.zeros(size)
        for i in np.nditer(numbers, op_flags=['readwrite']):
            i[...] = self.rand(shape=1)[0] + shift
        return numbers

    def multivariate_normal(self, mean, cov, cho=None, size=1):
        """Draw numbers from a multi-variate distribution.

        This is an attempt on speeding up the call of
        `RandomState.multivariate_normal` if we need to call it over and
        over again. Such repeated calling will do an SVD repeatedly,
        which is wasteful. In this function, this transform can be
        supplied and it is only estimated if it's not explicitly given.

        Parameters
        ----------
        mean : numpy array (1D, 2)
            Mean of the N-dimensional array
        cov : numpy array (2D, (2, 2))
            Covariance matrix of the distribution.
        cho : numpy.array (2D, (2, 2)), optional
            Cholesky factorization of the covariance. If not given,
            it will be calculated here.
        size : int, optional
            The number of samples to do.

        Returns
        -------
        out : float or numpy.array of floats size
            The random numbers that are drawn here.

        See also
        --------
        numpy.random.multivariate_normal

        """
        norm = self.normal(loc=0.0, scale=1.0, size=2*size)
        norm = norm.reshape(size, 2)
        meanm = np.array([mean, ] * size)
        return 0.01*(meanm + norm)

    def choice(self, a, size=None, replace=True, p=None):
        """Choose random samples from a given array.

        Parameters
        ----------
        a : 1-D array-like or int
            If an ndarray, a random sample is generated from its
            elements. If an int, the random sample is generated as if
            a were np.arange(a).
        size : int or tuple of ints, optional
            Output shape. Default is None, in which case a single
            value is returned.
        replace : boolean, optional
            Whether the sample is with or without replacement.
        p : 1-D array-like, optional
            The probabilities associated with each entry in a. If not
            given the sample assumes a uniform distribution over all
            entries in a.

        Returns
        -------
        choice : any or numpy.array
            The generated random samples.

        Note
        ----
        This is part of a mock generator and uses the pre-defined
        sequence to make choices.

        """
        if isinstance(a, int):
            a = list(range(a))
        if p is None:
            p = [1/len(a) for i in a]
        if size is None:
            size = 1
        out = []
        for _ in range(size):
            rand = self.rand()
            prob0 = 0
            # Make sure p sums to 1 even after popping
            p_sum = sum(p)
            p = [i/p_sum for i in p]
            for i, probi in zip(a, p):
                prob0 += probi
                if rand < prob0:
                    out.append(i)
                    if not replace:
                        idx = a.index(i)
                        _ = a.pop(idx)
                        _ = p.pop(idx)
                    break
        return out


class Borg:
    """A class for sharing states of objects."""

    class_state = None
    number_of_borgs = 0

    @classmethod
    def update_state(cls, state):
        """Update the class state and enable sharing of it."""
        if cls.class_state is None:
            logger.debug('Setting state to the shared.')
            cls.class_state = state
        else:
            logger.debug('Reusing shared state.')
        return cls.class_state

    @classmethod
    def reset_state(cls):
        """Remove memory of the shared state."""
        cls.class_state = None

    @classmethod
    def make_new_swarm(cls):
        """Make new swarm (new class)."""
        name = cls.__name__ + str(cls.number_of_borgs)
        cls.number_of_borgs += 1
        new_borg = type(name, cls.__bases__, dict(cls.__dict__))
        new_borg.reset_state()  # pylint: disable=no-member
        return new_borg

    def __init__(self, *args, **kwargs):
        """Enable share of the state.

        This will update the Borg.class_state and
        update self.__dict__ to this object, which enables
        sharing of the state.

        """
        super().__init__(*args, **kwargs)
        self.__dict__ = self.update_state(self.__dict__)


class RandomGeneratorBorg(Borg, RandomGenerator):
    """A class for sharing the state between RandomGenerator objects."""

    def get_state(self):
        """Return current state including the Borg generator type."""
        state = RandomGenerator.get_state(self)
        state['rgen'] = 'rgen-borg'
        return state


class MockRandomGeneratorBorg(Borg, MockRandomGenerator):
    """A class for sharing the state between MockRandomGenerator objects."""

    def get_state(self):
        """Return current state including the Borg generator type."""
        state = MockRandomGenerator.get_state(self)
        state['rgen'] = 'mock-borg'
        return state


def create_random_generator(settings=None):
    """Create a random generator from given settings.

    Parameters
    ----------
    settings : dict, optional
        This is the dict used for creating the random generator.
        Currently, we will actually just look for a seed value.

    Returns
    -------
    out : object like :py:class:`.RandomGenerator`
        The random generator created.

    """
    if settings is None:
        settings = {}

    rgen = settings.get('rgen', 'rgen')
    seed = settings.get('seed', 0)
    logger.debug('Seed for random generator: %s %s', rgen, seed)

    class_map = {
        'rgen': RandomGenerator,
        'rgen-borg': RandomGeneratorBorg,
        'mock': MockRandomGenerator,
        'mock-borg': MockRandomGeneratorBorg
    }
    rgen_class = class_map.get(rgen, RandomGenerator)
    rgen = rgen_class(seed=seed)

    if 'state' in settings:
        rgen.set_state(settings)
        rgen.status = 'restarted'
    else:
        rgen.status = 'new rgen'

    return rgen
