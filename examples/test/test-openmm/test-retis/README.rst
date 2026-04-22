Testing the OpenMM engine by running a RETIS simulation
=========================================================

This test runs a RETIS simulation with the OpenMM engine and compares
the results to previously stored reference output.

The folder ``openmm`` contains the simulation settings (``retis.rst``)
and the OpenMM system definition (``openmm_sim.py``), along with a
``results/`` sub-directory holding the reference output files.

System
------

Two Argon atoms interacting via a Lennard-Jones potential, propagated
by a Langevin integrator at 300 K (2 fs time-step, 10 sub-cycles per
PyRETIS step).

Instructions
------------

Run the test with:

.. code-block:: bash

   bash run.sh

The script will

1. Run ``pyretisrun`` inside the ``openmm/`` directory.
2. Compare the output against the reference files in ``openmm/results/``
   using ``compare.py``.
3. Clean up generated files.
