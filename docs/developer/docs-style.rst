.. _developer-guide-docs-style:

Documentation Terminology
=========================

This page records a few naming and wording conventions for user-facing
documentation. The goal is to keep the documentation consistent across the
user guide, examples, and API landing pages.

Project Names
-------------

Use the configured substitutions in prose:

* Use ``|pyretis|`` for the project name.
* Use ``|pyvisa|`` for the visualization and analysis tool.
* Use ``pyretisrun``, ``pyretisanalyse``, and ``pyvisa`` for command names.

Use the literal package or module name when referring to Python imports, for
example ``pyretis.analysis`` or ``pyretis.orderparameter``.

Scientific and File Terms
-------------------------

Use these forms in prose:

* ``order parameter`` for the concept.
* ``orderparameter`` only for the input section, package, or module name.
* ``Monte Carlo`` for the method.
* ``HDF5`` for the format, and ``.hdf5`` or ``.hdf5.zip`` for extensions.
* ``reStructuredText`` for the markup language, and ``.rst`` for files.
* ``Fortran`` for the language, unless quoting source text that uses another
  spelling.

Style Notes
-----------

Prefer concise, direct instructions. For example, write "Run the analysis with"
instead of "The analysis can be executed using". Keep command examples in
literal blocks and wrap filenames, input keys, and command-line flags in
double backticks.

When a page describes a workflow, start with the result the user is trying to
get, then give the commands or settings needed to get there. If a step depends
on a previous example, link to that example near the start of the page.
