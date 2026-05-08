# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""PyRETIS.

This file is part of PyRETIS - a simulation package for rare events.

Copyright (c) 2026, PyRETIS Development Team

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


PyRETIS documentation
---------------------

The documentation for PyRETIS is available either doc-strings provided
with the code and from `the PyRETIS homepage <http://www.pyretis.org>`_.

PyRETIS sub-packages
--------------------

analysis (:py:mod:`pyretis.analysis`)
    Analysis tools for calculating crossing probabilities, rates etc.

core (py:mod:`pyretis.core`)
    Core classes and functions for running the rare event simulations.
    This includes classes defining the system, particles, simulations
    etc.

engines (:py:mod:`pyretis.engines`)
    A package defining engines which can be used to
    evolve the dynamics/solve Newton's equations of motion in time.

forcefield (:py:mod:`pyretis.forcefield`)
    This package defines force fields and potentials functions.

inout (:py:mod:`pyretis.inout`)
    This package defines the input-output operations for PyRETIS.
    This includes generating output from the analysis and reading
    input-files etc.

orderparameter (:py:mod:`pyretis.orderparameter`)
    Definition of classes for order parameters. Defines the base class
    for order parameters.

testing (:py:mod:`pyretis.testing`)
    This package defines common methods which are used in testing.

tools (:py:mod:`pyretis.tools`)
    This package defines some functions which can be useful for
    setting up simple systems, for example, functions for generating
    lattices.

pyvisa (:py:mod:`pyretis.pyvisa`)
    This package contains the compressor, the gui and analysis of
    the package PyVisA. This is component nearly independent from
    PyRETIS and has different prerequisites.
    Note: if PyQt5 is not installed, all PyVisA functions that
    are connected to the GUI are not usable. We thus only check here
    if PyQt5 is present.

"""
# PyRETIS imports:
from .version import VERSION as __version__
from .version import MIN_PYTHON_VERSION
from . import analysis
from . import core
from . import engines
from . import forcefield
from . import info
from . import initiation
from . import inout
from . import orderparameter
from . import setup
from . import testing
from . import tools
# PyVisA imports:
from .pyvisa import HAS_PYVISA_REQ
if HAS_PYVISA_REQ:
    from . import pyvisa  # pragma: no cover
