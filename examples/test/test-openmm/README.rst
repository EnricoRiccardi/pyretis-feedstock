Testing the OpenMM engine in PyRETIS
=====================================

This directory contains tests for the PyRETIS OpenMM engine.
The system used throughout is two Argon atoms interacting via a
Lennard-Jones potential, propagated with a Langevin integrator.

.. note::
   These tests require OpenMM to be installed.
   Install with: ``conda install -c conda-forge openmm``

Tests
-----

* ``test-retis``: Runs a RETIS simulation with the OpenMM engine and
  compares the results to previously stored reference output.

Running all tests
-----------------

.. code-block:: bash

   bash run-all-openmm.sh
