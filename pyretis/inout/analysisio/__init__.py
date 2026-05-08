# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This package contains functions for input/output for the analysis.

The functions defined here will typically run the analysis on the
given input and write outputs, typically this will be plots and simple
text files.

Package structure
-----------------

Modules
~~~~~~~

analysisio.py (:py:mod:`.analysisio`)
    Methods that will output results from the analysis functions.
    The methods defined here can also be used to run an analysis on
    output files from PyRETIS.

__init__.py
    This file handles imports for PyRETIS.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

analyse_file (:py:func:`.analyse_file`)
    Method to analyse a file.

run_analysis_files (:py:func:`.run_analysis_files`)
    Method to analyse simulation data in output files. This is a
    post-processing step.
"""
from .analysisio import analyse_file, run_analysis_files
