.. _user-guide-engine:

Creating custom engines
=======================

The Engines are used to propagate the equations of motion or
alter the state of the system in some other way. In this
section, we describe how additional engines can be added to |pyretis|.
Here, we make a distinction between two types of engines:

* :ref:`Internal engines <user-guide-engine-internal>`
  which directly interact with the underlying numpy.arrays
  representing positions, velocities and forces (see :numref:`fig-engine`-left)

* :ref:`External engines <user-guide-engine-external>` which are used by |pyretis|
  to control the execution of external programs. In this case, the state of the
  particles is represented by a **file containing the positions and velocities**.
  As shown in the illustration below (:numref:`fig-engine`-right) , this can be accessed using
  :py:attr:`.ParticlesExt.config` which is a tuple containing a file name and
  an index. The interpretation is that the current state of the system can be
  found in the given file at the given frame/index.
  Further, the attribute :py:attr:`.ParticlesExt.vel_rev` is used to determine
  if the velocities in the configuration should be reversed or not. In this case
  the engine will typically make use of :py:meth:`.ParticlesExt.set_pos` in
  order to update the state of the system.

.. _fig-engine:

.. figure:: /_static/img/enginescomp.png
    :width: 95%
    :alt: Comparison of internal and external engines.
    :align: center

    Illustration of how **internal engines** (left) and **external engines** (right)
    interact with the system and particles.
    (Left) The internal engine interacts with and alters
    the state of the system, and it can directly access the positions
    and velocities, as numpy.arrays, using :py:attr:`.Particles.pos` and
    :py:attr:`.Particles.vel`.
    (Right) The external engine interacts with and alters
    the state of the system using :py:meth:`.ParticlesExt.set_pos` and
    the state of the system is represented by the attribute
    :py:attr:`.ParticlesExt.config`.

.. contents:: Table of contents
   :local:

.. _user-guide-engine-internal:

Creating a new internal engine
------------------------------

In order to create a new internal engine, we create a new class
for our engine which is sub-classing the generic :py:class:`.EngineBase`.
Often, we wish to create a new integrator for molecular dynamics
and then we can simplify our work by sub-classing the more specific
:py:class:`.MDEngine`.  The steps you need to complete in order to
create a new MDengine is as follows:

1. Create a new Python class which is a sub-class of :py:class:`.MDEngine`

2. Write a method to initialise this class, that is a ``__init__`` method.
   This method will be called when |pyretis| is setting up a new simulation
   and it will be fed variables from the |pyretis| input file.

3. Write a method to perform a single integration step.

4. Modify the input script to use the new integrator.


Example: Creating a new internal LeapFrog integrator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To be specific, we will in the following show
the steps needed to create a new MD engine using the
Leapfrog scheme [1]_:

.. _user-engine-leapfrog-equations:

.. math::

   \mathbf{r} (t + \Delta t) &= \mathbf{r} (t) + \mathbf{v} (t) \Delta t + \tfrac{1}{2}
   \mathbf{a} (t) \Delta t^2 \\
   \mathbf{v} (t + \Delta t) &= \mathbf{v} (t) + \tfrac{1}{2} \left[\mathbf{a} (t) +
   \mathbf{a} (t+\Delta t) \right] \Delta t

where :math:`\mathbf{r}` are the positions, :math:`\mathbf{v}` the velocities,
:math:`\mathbf{a}` the accelerations, :math:`\Delta t` the time step and :math:`t` is the
current time.

Step 1 and 2: Sub-classing MDEngine
...................................

In order to define a new class to use with |pyretis| we
sub-class the :py:class:`.MDEngine` and define the
``__init__`` method:

.. literalinclude:: /_static/engine-examples/leapfrog.py
   :language: python
   :lines: 5-21

Step 3: Adding the integration step method
..........................................

Now, we add the actual integration routine:

.. literalinclude:: /_static/engine-examples/leapfrog.py
   :language: python
   :pyobject: LeapFrog.integration_step

Step 4: Making use of the new integrator
........................................

The new integrator can be used by specifying the following
:ref:`engine section <user-section-engine>`:

.. pyretis-input-example:: Engine
   :class-name: LeapFrog

   .. code-block:: rst

      Engine
      ------
      class = LeapFrog
      module = leapfrog.py
      timestep = 0.002

where the keyword ``module`` specifies the file where you have
stored the new ``LeapFrog`` class.

.. _user-guide-engine-external:

Creating a new external engine
------------------------------


The external engines are used to interface external programs.
In |pyretis| there is a generic :py:class:`.ExternalMDEngine` which
you can make use of in order to set up new custom external engines.
Before we give an overview of the methods you will have to implement
in order to create a new engine, we list some useful methods which
are already present in the generic :py:class:`.ExternalMDEngine`:

* :py:meth:`.ExternalMDEngine.execute_command` which can be used to execute
  system commands.

* Methods to interact with files:

  - :py:meth:`.ExternalMDEngine._movefile`

  - :py:meth:`.ExternalMDEngine._copyfile`

  - :py:meth:`.ExternalMDEngine._removefile`

  - :py:meth:`.ExternalMDEngine._remove_files`

* :py:meth:`.ExternalMDEngine._modify_input` which can be used to modify
  simple input files for the external software.

In order to create a new external engine, there are several methods that you need
to implement. These methods are described in the :ref:`next section <user-engine-external-methods>`.
Before reading this description, please have a look at the figure below which illustrates the
interaction between |pyretis| and the external engine:

.. _user-engine-external-methods:

Implementing methods for a new external engine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following methods must be implemented in the new engine:

* :ref:`__init__ <user-engine-external-init>` which does some set-up for the
  external engine.

* :ref:`step <user-engine-external-step>` which is a method for performing
  a single MD step with the external engine.

* :ref:`_read_configuration <user-engine-external-read-configuration>`
  For reading output (configurations) from the external engine.
  This is used for calculating the order parameter(s).

* :ref:`_reverse_velocities <user-engine-external-reverse-velocities>`
  For reversing velocities in a snapshot. This method
  will typically make use of the
  :ref:`_read_configuration <user-engine-external-read-configuration>`.

* :ref:`_extract_frame <user-engine-external-extract-frame>`
  For extracting a single frame from a trajectory.

* :ref:`_propagate_from <user-engine-external-progagate-from>`
  The method for propagating the equations of motion using
  the external engine.

* :ref:`modify_velocities <user-engine-external-modify-velocities>`
  The method used for generating random velocities for
  shooting points.

Some examples/templates are shown in the following sections.


.. _user-engine-external-init:

The __init__ method
...................

The ``__init__`` method sets up the engine and makes it ready for use.
This can, for instance, involve checking that the required input files are
in place and checking the consistency of the input files.


.. literalinclude:: /_static/engine-examples/newexternal.py
   :language: python
   :lines: 5-63

In addition, this is a good point to add a method to execute
the external program. Here we will make use of
:py:meth:`.ExternalMDEngine.execute_command` and the
:py:attr:`.ExternalMDEngine.exe_dir` property. This property
is set internally by |pyretis| and points to a directory which
we will use when executing the external engine. As a concrete
example, if we are generating a path for ensemble :math:`[1^+]`,
this will point to the sub-directory ``generate`` of the
ensemble-directory ``002``.
Here is a short template for executing external commands:

.. literalinclude:: /_static/engine-examples/newexternal.py
   :language: python
   :pyobject: NewEngine._run_program


.. _user-engine-external-step:

The step method
...............

The step method is used to execute a single step with the external software.
This method is used by the |pyretis| engine in the following way:

.. code-block:: python

   conf = self.step(system, 'output-file')

That is, the method takes in two parameters: a system object and a *string*
and returns a single parameter. This single parameter is just a string
which contains the output configuration created after executing the
external program. Here is a short template for the ``step`` method:

.. literalinclude:: /_static/engine-examples/newexternal.py
   :language: python
   :pyobject: NewEngine.step


.. _user-engine-external-read-configuration:

The _read_configuration method
..............................

The ``_read_configuration`` method is used to get the actual
positions and velocities from a configuration file. As
a short template:

.. literalinclude:: /_static/engine-examples/newexternal.py
   :language: python
   :pyobject: NewEngine._read_configuration


.. _user-engine-external-reverse-velocities:

The _reverse_velocities method
..............................

The ``_reverse_velocities`` method is used to reverse velocities for
a single configuration.

.. literalinclude:: /_static/engine-examples/newexternal.py
   :language: python
   :pyobject: NewEngine._reverse_velocities


.. _user-engine-external-extract-frame:

The _extract_frame method
.........................

The ``_extract_frame`` method is used to extract a single snapshot from
a trajectory file.

.. literalinclude:: /_static/engine-examples/newexternal.py
   :language: python
   :pyobject: NewEngine._extract_frame

.. _user-engine-external-progagate-from:

The _propagate_from method
..........................

The ``_propagate_from`` method is used to propagate the equations of motion
using the external engine. In some cases, this can make use of the
:ref:`step <user-engine-external-step>` method, in other cases this is
really distinct (e.g. in GROMACS this would correspond to extending a
simulation).

.. literalinclude:: /_static/engine-examples/newexternal.py
   :language: python
   :pyobject: NewEngine._propagate_from


.. _user-engine-external-modify-velocities:

The modify_velocities method
............................

This method is used to draw random velocities, typically when
performing the shooting move. The generic ``modify_velocities`` method
is quite involved, but for most cases is sufficient to only implement
aimless shooting. Here is a template:

.. literalinclude:: /_static/engine-examples/newexternal.py
   :language: python
   :pyobject: NewEngine.modify_velocities


References
----------

.. [1]  https://en.wikipedia.org/wiki/Leapfrog_integration
