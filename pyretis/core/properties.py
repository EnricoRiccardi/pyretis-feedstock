# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This file contains a class for a generic property."""
import numpy as np

__all__ = ['Property']

np.set_printoptions(legacy='1.21')


class Property:
    """A generic numerical value with standard deviation and average.

    A generic object to store values obtained during a simulation.
    It will maintain the mean and variance as values are added using
    Property.add(val)

    Attributes
    ----------
    desc : string
        Description of the property.
    nval : integer
        Number values stored.
    mean : float
        The current mean.
    delta2 : float
        Helper variable used for calculating the variance.
    variance : float
        The current variance
    val : list
        Store all values.

    Parameters
    ----------
    desc : string, optional
        Used to set the attribute desc.

    Examples
    --------
    >>> from pyretis.core.properties import Property
    >>> ener = Property(desc='Energy of the particle(s)')
    >>> ener.add(42.0)
    >>> ener.add(12.220)
    >>> ener.add(99.22)
    >>> ener.mean

    """

    def __init__(self, desc=''):
        """Initialise the property.

        Parameters
        ----------
        desc : string, optional
            Description of the object.

        """
        self.desc = desc
        self.nval = 0
        self.mean = 0.0
        self.delta2 = 0.0
        self.variance = 0.0
        self.val = []

    def add(self, val):
        """Add a value to the property & update the mean and variance.

        Parameters
        ----------
        val : float or another type (list, numpy.array)
            The value to add.

        Returns
        -------
        out : None
            Returns `None` but updates the mean and variance.

        """
        self.nval += 1
        self.val.append(val)
        self.update_mean_and_variance()

    def update_mean_and_variance(self):
        """Calculate the mean and variance on the fly.

        Source:
        http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance

        Returns
        -------
        out : None
            Returns `None` but updates the mean and variance.

        Note
        ----
        Consider if this should be moved/deleted and just
        replaced with a function from the analysis package.

        """
        val = self.val[-1]  # most recent value
        delta = val - self.mean
        self.mean += delta/float(self.nval)
        self.delta2 += delta * (val - self.mean)
        if self.nval < 2:
            self.variance = float('inf')
        else:
            self.variance = self.delta2/float(self.nval - 1)

    def dump_to_file(self, filename):
        """Dump the contents in `self.val` to a file.

        Parameters
        ----------
        filename : string
            Name/path of file to write.

        Note
        ----
        Consider if this should be moved/deleted and just
        replaced with a function from a more general input-output
        module

        """
        np.savetxt(filename, self.val)

    def __str__(self):
        """Return string description of the property."""
        return self.desc
