
.. _api-doc:

#####################
The |pyretis| library
#####################

This is the documentation for the |pyretis| library and usage
of the application programming interface (API).
An overview of the main classes and routines are
given in the :ref:`introduction to the library <user-guide-intro-api>` in
the :ref:`user guide <user-guide-index>`.
It might also be a good idea to have a look at the
:ref:`developer guide <developer-guide-index>`
before making changes or reading through the API :ref:`docs <api-doc>`.

The |pyretis| package
=====================

The |pyretis| library contains methods and classes that handle
the different parts of simulations. These are further grouped
into sub-packages:

* :ref:`pyretis.analysis <api-analysis>` for analysing the
  output from simulations.

* :ref:`pyretis.core <api-core>` with the core functions for the
  different methods.

* :ref:`pyretis.engines <api-engines>` for integrating
  the equations of motion and for interfacing
  external MD software.

* :ref:`pyretis.forcefield <api-forcefield>` for defining
  force fields to use in simulations.

* :ref:`pyretis.initiation <api-initiation>` with functions for
  initiation of path sampling simulations.

* :ref:`pyretis.inout <api-inout>` for handling the input
  and output to |pyretis|.

* :ref:`pyretis.orderparameter <api-orderparameter>` for defining
  order parameters.

* :ref:`pyretis.setup <api-setup>` for creating the main objects.

* :ref:`pyretis.simulation <api-simulation>` for running predefined
  simulations.

* :ref:`pyretis.testing <api-testing>` with functions which are used
  in connection with testing.

* :ref:`pyretis.tools <api-tools>` for performing some simple
  tasks useful for setting up simulations.

* :ref:`pyretis.pyvisa <api-pyvisa>` for visualizing and
  analysing the simulation outputs. 


The |pyretis| library can be imported by:

.. code-block:: python

    import pyretis
    print(pyretis.__version__)

Or you can just import specific functions and classes from
the sub-packages:

.. literalinclude:: /_static/api/box_example.py
   :language: python
   :lines: 5-18

The usage of these sub-packages are described in the following:

Contents: The |pyretis| sub-packages
------------------------------------

.. toctree::
    :maxdepth: 2

    pyretis.analysis
    pyretis.core
    pyretis.engines
    pyretis.forcefield
    pyretis.inout
    pyretis.initiation
    pyretis.orderparameter
    pyretis.setup
    pyretis.simulation
    pyretis.testing
    pyretis.tools
    pyretis.pyvisa

In addition |pyretis| has the following submodules:

List of submodules
------------------

There are also some submodules with information about
|pyretis|. These are:

* :ref:`pyretis.info <api-info>`
* :ref:`pyretis.version <api-version>`

.. _api-info:

pyretis.info module
...................

.. automodule:: pyretis.info
    :members:
    :undoc-members:
    :show-inheritance:

.. _api-version:

pyretis.version module
......................

.. automodule:: pyretis.version
    :members:
    :undoc-members:
    :show-inheritance:

|pyretis| executables
---------------------

The main |pyretis| executables can be found in:

.. toctree::
    :maxdepth: 1

    pyretis.bin
