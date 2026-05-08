# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This package contains methods for generating reports.

The reports will typically summarise the results from different
analysis and present it as a text file, pdf or web-page.

Package structure
-----------------

Modules
~~~~~~~

__init__.py
    This file. Handles imports for PyRETIS and defines a method for
    writing a report to a file.

markup.py (:py:mod:`pyretis.inout.report.markup`)
    This module defines some methods for generating simple tables and
    formatting numbers for rst and latex.

report.py (:py:mod:`pyretis.inout.report.report`)
    General methods for generating reports.

report_md.py (:py:mod:`pyretis.inout.report.report_md`)
    This module defines the molecular dynamics reports. Specifically
    it defines the report that is made based on results from a MD Flux
    simulations.

report_path.py (:py:mod:`pyretis.inout.report.report_path`)
    This module defines the reports for path simulations like TIS and
    RETIS.


Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

generate_report (:py:func:`.generate_report`)
    Method for generating reports.

Folders
~~~~~~~

templates
    A folder containing templates for generating reports.
"""
from .report import generate_report
