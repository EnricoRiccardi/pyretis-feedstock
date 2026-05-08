# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A collection of potential functions.

This package defines some potential functions. These potential functions
can be used to create force fields.


Package structure
-----------------

Modules
~~~~~~~

potentials.py (:py:mod:`pyretis.forcefield.potentials.potentials`)
    This module defines some simple potential functions.

Sub-packages
~~~~~~~~~~~~

pairpotentials (:py:mod:`pyretis.forcefield.potentials.pairpotentials`)
    This package defines different pair interactions, for instance the
    Lennard-Jones 6-12 simple cut potential.

Important classes defined in this package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DoubleWell (:py:class:`.DoubleWell`)
    A double well potential

RectangularWell (:py:class:`.RectangularWell`)
    A rectangular well potential -- useful as a bias potential.

PairLennardJonesCut (:py:class:`.PairLennardJonesCut`)
    The Lennard-Jones potential in pure python.

PairLennardJonesCutnp (:py:class:`.PairLennardJonesCutnp`)
    The Lennard-Jones potential, making use of numpy.

DoubleWellWCA (:py:class:`.DoubleWellWCA`)
    A n-dimensional Double Well potential.
"""
from .potentials import DoubleWell, RectangularWell
from .pairpotentials import (PairLennardJonesCut, PairLennardJonesCutnp,
                             DoubleWellWCA)
