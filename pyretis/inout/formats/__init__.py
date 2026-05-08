# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This sub-package handles data formats for PyRETIS.

The methods and classes defined in this package will format the data
created by PyRETIS and also enable the reading of such output data.

Further, this package includes some methods for interacting with other
simulation packages.

Package structure
-----------------

Modules
~~~~~~~

formatter.py (:py:mod:`pyretis.inout.formats.formatter`)
    Generic methods for formatting output from PyRETIS. Defines
    the base class for all formatters.

cp2k.py (:py:mod:`pyretis.inout.formats.cp2k`)
    A module defining input/output methods for use with CP2k.

cross.py (:py:mod:`pyretis.inout.formats.cross`)
    Formatting of crossing data from PyRETIS.

energy.py (:py:mod:`pyretis.inout.formats.energy`)
    Formatting of energy data from PyRETIS.

gromacs.py (:py:mod:`pyretis.inout.formats.gromacs`)
    A module defining input/output methods for use with GROMACS.

__init__.py
    This file.

order.py (:py:mod:`pyretis.inout.formats.order`)
    Formatting of order parameter data from PyRETIS.

path.py (:py:mod:`pyretis.inout.formats.path`)
    Formatting of path/trajectory data from PyRETIS.

pathensemble.py (:py:mod:`pyretis.inout.formats.pathensemble`)
    Formatting of path ensemble data from PyRETIS.

snapshot.py (:py:mod:`pyretis.inout.formats.snapshot`)
    Formatting of snapshot data from PyRETIS.

txt_table.py (:py:mod:`pyretis.inout.formats.txt_table`)
    Definition of a table-like format.

xyz.py (:py:mod:`pyretis.inout.formats.xyz`)
    A module defining input/output methods for use with a XYZ format.


Important classes defined in this package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CrossFormatter (:py:class:`.CrossFormatter`)
    A class for formatting crossing data from flux simulations.

EnergyFormatter (:py:class:`.EnergyFormatter`)
    A class for formatting energy data from PyRETIS.

EnergyPathFormatter (:py:class:`.EnergyPathFormatter`)
    A class for formatting energy data for paths.

OrderFormatter (:py:class:`.OrderFormatter`)
    A class for formatting order parameter data.

OrderPathFormatter (:py:class:`.OrderPathFormatter`)
    A class for formatting order parameter data for paths.

PathExtFormatter (:py:class:`.PathExtFormatter`)
    A class for formatting external trajectories.

PathIntFormatter (:py:class:`.PathIntFormatter`)
    A class for formatting internal trajectories.

PathEnsembleFormatter (:py:class:`.PathEnsembleFormatter`)
    A class for formatting path ensemble data.

SnapshotFormatter (:py:class:`.SnapshotFormatter`)
    Generic class for formatting system snapshots (coordinates).

TxtTableFormatter (:py:class:`.TxtTableFormatter`)
    A class for creating a table-like format.

"""
import logging
from pyretis.inout.formats.cross import CrossFormatter, CrossFile
from pyretis.inout.formats.energy import (
    EnergyFormatter,
    EnergyPathFormatter,
    EnergyFile,
    EnergyPathFile,
)
from pyretis.inout.formats.order import (
    OrderFormatter,
    OrderPathFormatter,
    OrderFile,
    OrderPathFile,
)
from pyretis.inout.formats.path import (
    PathExtFormatter,
    PathIntFormatter,
    PathExtFile,
    PathIntFile,
)
from pyretis.inout.formats.pathensemble import (
    PathEnsembleFormatter,
    PathEnsembleFile,
)
from pyretis.inout.formats.snapshot import (
    SnapshotFormatter,
    SnapshotFile,
)
from pyretis.inout.formats.txt_table import (
    TxtTableFormatter,
    PathTableFormatter,
    ThermoTableFormatter,
    RETISResultFormatter,
)
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())
