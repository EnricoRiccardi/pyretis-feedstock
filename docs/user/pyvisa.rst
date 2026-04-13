.. _user-guide-pyvisa:

The |pyvisa| application
========================

|pyvisa| is the analysis and visualization tool for |pyretis| simulations.
It consists of two components:

* A **compressor** that reads raw simulation output, checks consistency,
  and stores the data in a compressed ``.hdf5`` binary file.
* A **visualizer** (GUI) that loads compressed or raw data and produces
  interactive plots, cluster analyses, and PCA results.
  The GUI requires PyQt5, which must be installed separately
  (see :ref:`installation <user-guide-install>`).

The general syntax is:

.. code-block:: pyretis

   pyvisa [-h] [-i INPUT] [-V] [-cmp] [-data DATA] [-recalculate] [-oo]

where the arguments are described in :numref:`tableappargument_pyvisa`.

Example use
-----------

* Compress raw simulation output into a ``.hdf5`` file:

  .. code-block:: pyretis

     pyvisa -i out.rst -cmp

* Compress using only order parameter files (faster, smaller output):

  .. code-block:: pyretis

     pyvisa -i out.rst -cmp -oo

* Open the visualization GUI from a compressed file:

  .. code-block:: pyretis

     pyvisa -i simulation.hdf5

* Open the GUI from the original input file, loading a specific ensemble:

  .. code-block:: pyretis

     pyvisa -i out.rst -data 000

* Recalculate order parameters and collective variables:

  .. code-block:: pyretis

     pyvisa -i retis.rst -recalculate

For a guided walkthrough of these features, see
:ref:`examples-pyvisa` and :ref:`examples-pyvisa-analysis`.

Input arguments
---------------

.. _tableappargument_pyvisa:

.. table:: Description of input arguments for |pyvisa|.
   :class: table-striped table-hover

   +-----------------------------------------------+--------------------------------------------------+
   | Argument                                      | Description                                      |
   +===============================================+==================================================+
   | -h, --help                                    | Show the help message and exit.                  |
   +-----------------------------------------------+--------------------------------------------------+
   | -i INPUT, --input INPUT                       | Location of the input file (``.rst`` or          |
   |                                               | ``.hdf5``).                                      |
   +-----------------------------------------------+--------------------------------------------------+
   | -V, --version                                 | Show the version number and exit.                |
   +-----------------------------------------------+--------------------------------------------------+
   | -cmp, --pyvisa_compressor                     | Run the compressor tool to generate a            |
   |                                               | ``.hdf5`` file from raw simulation output.       |
   +-----------------------------------------------+--------------------------------------------------+
   | -data DATA, --pyvisa-data DATA                | Select a subset of data to load. Use ``all``     |
   |                                               | for all data, or a folder name (e.g. ``000``)    |
   |                                               | for a single ensemble.                           |
   +-----------------------------------------------+--------------------------------------------------+
   | -recalculate, --pyvisa-recalculate            | Recalculate order parameter and collective       |
   |                                               | variable data using the current ``orderp.py``.   |
   +-----------------------------------------------+--------------------------------------------------+
   | -oo, --only_order                             | Use only data from ``order.txt`` files when      |
   |                                               | compressing (faster, skips energy data).         |
   +-----------------------------------------------+--------------------------------------------------+
