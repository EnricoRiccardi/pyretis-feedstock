.. _user-guide-input:

Running |pyretis|
=================

|pyretis| is executed using the :ref:`pyretisrun <user-guide-application>` application and
a |pyretis| input file:

.. code-block:: pyretis

   pyretisrun -i INPUT

where ``INPUT`` is the input file. This will produce output files
which can be analysed using the :ref:`pyretisanalyse <user-guide-analyse>`
application.

In the following, we describe the syntax for the |pyretis| input file.

The |pyretis| input file
------------------------

|pyretis| simulations can be set up and run with a simple input file.
The input file defines a simulation by setting values for
keywords which are organised into sections.
Here we will discuss the following:

.. contents::
   :local:

.. _user-guide-input-structure:

Structure of the input file
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The input file is organised into **sections** and for each
section, **keywords** are used to define your settings.
The syntax for setting keywords is ``keyword = setting``,
and sections are marked by the section title followed by an
underline of dashes ``--------``, for example:

.. code-block:: rst

    Simulation settings
    -------------------
    task = md-nve
    steps = 100


which sets the two keywords ``task`` and ``steps`` for the
:ref:`Simulation <user-section-simulation>` section.
When specifying a section, only the
first word in the section title is used to
identify the title internally in |pyretis|. This means that you can
just as well write the following:

.. code-block:: rst

    Simulation
    ----------
    task = md-nve
    steps = 100

or add any text you like, e.g.:

.. code-block:: rst

    Simulation master plan
    ----------------------
    task = md-nve
    steps = 100

Both of these two examples define keywords for
the :ref:`Simulation <user-section-simulation>` section.

Formatting of keywords
^^^^^^^^^^^^^^^^^^^^^^
The format of the input file is relatively free,
you can for instance order things within sections as
you prefer and the input is in general **case-insensitive**:

.. code-block:: rst

    Simulation settings
    -------------------
    task = md-nve
    units = lj

which is identical to:

.. code-block:: rst

    Simulation settings
    -------------------
    UNITS = lj
    tAsK = md-nve

There are some important exceptions where the keyword settings
are in fact **case-sensitive**:

- When using external Python modules and classes, for instance:

  .. code-block:: rst

      Engine settings
      ---------------
      class = MyExternalClass
      module = filename.py

  Here, the values given for the ``class`` and the ``module`` keywords are
  **case-sensitive**.

- When referring to external files, for instance:

  .. code-block:: rst

      Particles settings
      ------------------
      position = {'input_file': 'myfile.gro'}

  Here, we are referring to a file named ``myfile.gro``,
  and |pyretis| will expect this file to be present with exactly
  this file name.

- When defining your own system of units:

  .. code-block:: rst

      Unit-system
      -----------
      length = (1.0, 'm')

  Here, we are using a unit of 1 meter which is identified with
  a ``'m'`` and not a ``'M'``.

Comments
^^^^^^^^

You can also add comments to structure the input file:

.. code-block:: rst

    Simulation settings
    -------------------
    task = md-nve
    units = lj

    # More settings:

    System settings
    ---------------
    temperature = 1.0

Comments are marked as starting with a ``'#'`` and all following text
will be ignored, i.e.

.. code-block:: rst

    task = md-nve  # set up and run a md-nve simulation not TIS this time.

is effectively the same as writing,

.. code-block:: rst

    task = md-nve

Summary
^^^^^^^

* The input file is organised into ``sections`` where ``keywords`` are set:

  .. code-block:: rst

      SectionTitle
      ------------
      keyword = value

* Comments are marked with a ``#``.

* Input is in general not case-sensitive unless you are referring to
  files and Python classes.


The sections in the input file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following sections are recognised by |pyretis|:

* :ref:`simulation <user-section-simulation>`:
  For defining the simulation we are going to run.

* :ref:`system <user-section-system>`:
  For defining system properties.

* :ref:`box <user-section-box>`:
  For defining a simulation box.

* :ref:`particles <user-section-particles>`:
  For defining the initial state of particles.

* :ref:`forcefield <user-section-forcefield>`:
  For defining a forcefield.

* :ref:`potential <user-section-potential>`:
  For defining potential functions to use in the force field.

* :ref:`engine <user-section-engine>`:
  For defining the simulation engine.

* :ref:`orderparameter <user-section-orderparameter>`:
  For defining the order parameter.

* :ref:`retis <user-section-retis>`:
  For defining settings for a RETIS simulation.

* :ref:`initial-path <user-section-initial-path>`:
  For defining settings to initialize a RETIS simulation.

* :ref:`tis <user-section-tis>`:
  For defining settings for a TIS simulation.

* :ref:`output <user-section-output>`:
  For defining output settings.

* :ref:`unit-system <user-section-unit-system>`:
  For defining custom unit systems.

In addition, there are analysis specific settings which can be set
by making use of the following section(s):

* :ref:`analysis <user-section-analysis>`:
  For defining an analysis.
