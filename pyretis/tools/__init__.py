# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Tools which can help with setting up simulations.

This package defines some simple tools which may be useful for
creating simulations.

Package structure
-----------------

Modules
~~~~~~~

lattice.py (:py:mod:`pyretis.tools.lattice`)
    Defines tools for setting up and generating lattice points.

recalculate_order.py (:py:mod:`pyretis.tools.recalculate_order`)
    Methods for recalculating order parameters on external paths.

Important methods defined in this package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

generate_lattice (:py:func:`.generate_lattice`)
    Generates points on a lattice.

recalculate_order (:py:func:`.recalculate_order`)
    Recalculate order parameter(s).
"""
from pyretis.tools.lattice import generate_lattice
from pyretis.tools.recalculate_order import recalculate_order
