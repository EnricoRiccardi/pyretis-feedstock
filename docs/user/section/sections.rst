.. _user-section:

Input file sections
===================

The |pyretis| input file described in detail
in :ref:`the user guide <user-guide-input-structure>`. The short version is:

1. The input file is organised into ``sections`` where ``keywords`` are given values:

   .. code-block:: rst

       Section Title
       -------------
       keyword = value

2. Comments are marked with a ``#``.

3. Input is in general **not** case-sensitive **unless** you are
   referring to **files** and **Python classes**.

.. toctree::
    :maxdepth: 1
    :hidden:

    simulation.rst
    system.rst
    box.rst
    particles.rst
    forcefield.rst
    potential.rst
    engine.rst
    orderparameter.rst
    retis.rst
    tis.rst
    initial.rst
    output.rst
    unitsystem.rst
    analysis.rst

Below, we list the different sections that you can make use of
in order to define your simulation:

.. |sec_simulation| replace:: :ref:`simulation <user-section-simulation>`

.. |sec_system| replace:: :ref:`system <user-section-system>`

.. |sec_box| replace:: :ref:`box <user-section-box>`

.. |sec_particles| replace:: :ref:`particles <user-section-particles>`

.. |sec_forcefield| replace:: :ref:`forcefield <user-section-forcefield>`

.. |sec_potential| replace:: :ref:`potential <user-section-potential>`

.. |sec_engine| replace:: :ref:`engine <user-section-engine>`

.. |sec_orderparameter| replace:: :ref:`orderparameter <user-section-orderparameter>`

.. |sec_retis| replace:: :ref:`retis <user-section-retis>`

.. |sec_tis| replace:: :ref:`tis <user-section-tis>`

.. |sec_path| replace:: :ref:`initial-path <user-section-initial-path>`

.. |sec_output| replace:: :ref:`output <user-section-output>`

.. |sec_usystem| replace:: :ref:`unit-system <user-section-unit-system>`

.. _table-sections:

.. table:: Input sections for defining simulations.
   :class: table-hover table-striped

   +----------------------+---------------------------------------------------+
   | Section              | Usage                                             |
   +======================+===================================================+
   | |sec_simulation|     | For defining the simulation we are going to run.  |
   +----------------------+---------------------------------------------------+
   | |sec_system|         | For defining system properties.                   |
   +----------------------+---------------------------------------------------+
   | |sec_box|            | For defining a simulation box.                    |
   +----------------------+---------------------------------------------------+
   | |sec_particles|      | For defining the initial state of particles.      |
   +----------------------+---------------------------------------------------+
   | |sec_forcefield|     | For defining a forcefield.                        |
   +----------------------+---------------------------------------------------+
   | |sec_potential|      | For defining potential functions to use in the    |
   |                      | force field.                                      |
   +----------------------+---------------------------------------------------+
   | |sec_engine|         | For defining the simulation engine.               |
   +----------------------+---------------------------------------------------+
   | |sec_orderparameter| | For defining the order parameter.                 |
   +----------------------+---------------------------------------------------+
   | |sec_retis|          | For defining settings for a RETIS simulation.     |
   +----------------------+---------------------------------------------------+
   | |sec_tis|            | For defining settings for a TIS simulation.       |
   +----------------------+---------------------------------------------------+
   | |sec_path|           | For defining how the initial path is generated.   |
   +----------------------+---------------------------------------------------+
   | |sec_output|         | For defining output settings.                     |
   +----------------------+---------------------------------------------------+
   | |sec_usystem|        | For defining custom unit systems.                 |
   +----------------------+---------------------------------------------------+

In addition, an analysis can be defined using:

.. |sec_analysis| replace:: :ref:`analysis <user-section-analysis>`

.. _table-sections-analysis:

.. table:: Input sections for defining an analysis.
   :class: table-hover table-striped

   +----------------------+---------------------------------------------------+
   | Section              | Usage                                             |
   +======================+===================================================+
   | |sec_analysis|       | For defining an analysis.                         |
   +----------------------+---------------------------------------------------+

Notation for describing keywords
--------------------------------

In each of these sections, the keywords are described using the following
notation:

.. pyretis-keyword:: KEYWORD DATA-TYPE

   Description of the keyword.

   Default:
     Description of default settings.

Here, KEYWORD is the actual keyword that is set, and *DATA-TYPE* is the allowed
parameter type for the particular keyword. The types you may encounter are
described in the :ref:`table below <table_data_types>`.


.. _table_data_types:

.. table:: The different data types encountered in |pyretis|.
   :class: table-hover table-striped

   +-------------------+--------------------------+-----------------------------------+
   |  *DATA-TYPE*      | Description              | Example                           |
   +===================+==========================+===================================+
   | *string*          | A string of characters,  |  ``task = retis``                 |
   |                   | i.e. text.               |                                   |
   +-------------------+--------------------------+-----------------------------------+
   | *integer*         | An integer.              | ``steps = 100``                   |
   +-------------------+--------------------------+-----------------------------------+
   | *float*           | A floating-point number  | ``timestep = 0.002``              |
   +-------------------+--------------------------+-----------------------------------+
   | *boolean*         | A boolean value          | ``shift = True``                  |
   |                   | (``True`` or ``False``). |                                   |
   +-------------------+--------------------------+-----------------------------------+
   |  *dictionary*     | A Python dictionary.     | ``mass = {'Ar': 1.0}``            |
   +-------------------+--------------------------+-----------------------------------+
   | *list*            | A Python list.           | ``interfaces = [0.1, 0.2, 0.3]``  |
   +-------------------+--------------------------+-----------------------------------+
   | *tuple*           | A Python tuple.          | ``index = (7,8)``                 |
   +-------------------+--------------------------+-----------------------------------+
   | *None*            | This represents an       |                                   |
   |                   | optional value           |                                   |
   +-------------------+--------------------------+-----------------------------------+

