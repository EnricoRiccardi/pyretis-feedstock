# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Definition of the base class for the plotter.

This module just defines a base class for plotters. This is just to
ensure that all plotters at least implement some functions and that we
can make use of them.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Plotter (:py:class:`.Plotter`)
    A generic class for creating plots.
"""
from abc import ABCMeta, abstractmethod

__all__ = ['Plotter']


class Plotter(metaclass=ABCMeta):
    """Base class for PyRETIS plotters.

    This class defines a plotter. A plotter is just an object
    that supports certain functions which can be called by analysis
    output functions. It should define certain plots and the
    ``Plotter`` class is an abstract class just to make sure that
    all plotters define the needed plots.

    Attributes
    ----------
    backup : boolean
        Determines if we overwrite old files or try to back them up.
    plotter_type : string
        Defines a name for the plotter, in case we want to identify it.
    out_dir : string, optional
        Defines an output directory for the plotter.

    """

    def __init__(self, backup=True, plotter_type=None, out_dir=None):
        """Initialise the plotting object.

        Parameters
        ----------
        backup : boolean, optional
            Determines if we overwrite old files or not.
        plotter_type : string, optional
            A name for the plotter.
        out_dir : string, optional
            A string which can be used to set an output directory
            for the plotter.

        """
        self.plotter_type = plotter_type
        self.backup = backup in (True, 'yes', 'True')
        self.out_dir = out_dir

    @abstractmethod
    def output_flux(self, results):
        """Plot flux results."""
        return

    @abstractmethod
    def output_energy(self, results, energies):
        """Plot energy results."""
        return

    @abstractmethod
    def output_orderp(self, results, orderdata):
        """Plot order parameter results."""
        return

    @abstractmethod
    def output_path(self, results, path_ensemble):
        """Plot results for each path ensemble."""
        return

    @abstractmethod
    def output_matched_probability(self, path_ensembles, detect, matched):
        """Plot the overall probability for path ensembles."""
        return

    def __str__(self):
        """Print out the basic info."""
        return f'Plotter: {self.plotter_type}'
