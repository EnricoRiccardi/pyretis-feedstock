.. _user-guide-analyse:

The |pyretis| analysis application
==================================

The |pyretis| analysis application, ``pyretisanalyse``, is used to
analyse the results from |pyretis| simulations.
The general syntax for executing is:

.. code-block:: pyretis

   pyretisanalyse [-h] [-i INPUT] [-V]


where ``INPUT`` is the input file for the analysis. It can be of two types:
An hdf5 file, that can be used with the PyVisA tool, or the 
standard input file used by |pyretis| 
:ref:`application <user-guide-application>` with the 
settings for the analysis, defined via the 
:ref:`analysis section <user-section-analysis>`.


.. container:: panel panel-info

   .. container:: panel-heading
   
      **Note**

   .. container:: panel-body
   
      The ``INPUT`` file for the analysis is created automatically
      by running the |pyretis| application. This file is named ``out.rst``
      and can be directly used as input for the analysis program.


where ``INPUT`` is the input file for the analysis. This file is
similar to the input file to the |pyretis|
:ref:`application <user-guide-application>` with some differences:

1. It contains the actual number of cycles/steps completed,
   specified by the
   :ref:`endcycle <user-section-simulation-endcycle>`
   keyword.

2. It contains the directory from where the simulation was
   executed, specified by the
   :ref:`exe-path <user-section-simulation-exe-path>` keyword.

3. It contains the number of particles in the simulation,
   specified by the 
   :ref:`npart <user-section-particles-keyword-npart>` keyword.

4. It contains settings for the analysis, defined via the 
   :ref:`analysis section <user-section-analysis>`.


.. container:: panel panel-info

   .. container:: panel-heading
   
      **Note**

   .. container:: panel-body
   
      The ``INPUT`` file for the analysis is created automatically
      by running the |pyretis| application. This file is named ``out.rst``
      and can be directly used as input for the analysis program.


Input arguments
---------------

.. _tableappargument_analyse:

.. table:: Description of input arguments for pyretisanalyse
   :class: table-striped table-hover

   +-------------------------------------+--------------------------------------------------+
   | Argument                            | Description                                      |
   +=====================================+==================================================+
   | -h, --help                          | Show the help message and exit.                  |
   +-------------------------------------+--------------------------------------------------+
   | -i INPUT, --input INPUT             | Location of the input file containing analysis   |
   |                                     | settings.                                        |
   +-------------------------------------+--------------------------------------------------+
   | -V, --version                       | Show the version number and exit.                |
   +-------------------------------------+--------------------------------------------------+
