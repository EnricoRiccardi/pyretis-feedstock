# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Package defining classes for pair potentials.

This package defines different pair potentials for use with an internal
calculation in PyRETIS.

Package structure
-----------------

Modules
~~~~~~~

lennardjones.py (:py:mod:`.lennardjones`)
    Potential functions for Lennard-Jones interactions.

pairpotential.py (:py:mod:`.pairpotential`)
    This module defines some helper functions for pair potentials.

wca.py (:py:mod:`.wca`)
    Potential functions for WCA-type interactions.

Important classes defined in this package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PairLennardJonesCut (:py:class:`.PairLennardJonesCut`)
    A class defining a Lennard-Jones potential.

PairLennardJonesCutnp (:py:class:`.PairLennardJonesCutnp`)
    A class defining a Lennard-Jones potential using numpy for the
    evaluation.

DoubleWellWCA (:py:class:`.DoubleWellWCA`)
    This class defines a double well WCA potential.
"""
from .lennardjones import PairLennardJonesCut, PairLennardJonesCutnp
from .wca import DoubleWellWCA
from .pairpotential import generate_pair_interactions
