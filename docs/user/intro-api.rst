.. _user-guide-intro-api:

Introduction to the |pyretis| library
=====================================

In this introduction to the |pyretis| library
the main classes and functions
from the |pyretis| library
will be discussed.
For the full description, we refer to the :ref:`API documentation <api-doc>`.

The |pyretis| library contains methods and classes that handle
the different aspects of a simulation. These are grouped
into sub-packages and the different
sub-packages and the classes and methods defined within them will interact.
As an example, assume that we are performing a RETIS simulation.
In terms of objects from the |pyretis| library,
this simulation can be described as follows:

1. We first define the :py:class:`.System` we are studying. This contains
   information about the :py:class:`.Particles`, :py:class:`.Box` and
   the :py:class:`.ForceField`.
2. Next, an :py:class:`.OrderParameter` is defined
   and the order parameter it is representing can be calculated using a
   :py:class:`.System` object as the argument.
3. The RETIS simulation is handled by the :py:class:`.SimulationRETIS` class
   which will use a specific :py:class:`.EngineBase` like object in order to
   generate several :py:class:`.Path` objects for a collection of
   :py:class:`.PathEnsemble` objects.
   This generation is done by a set of methods defined in the
   modules :py:mod:`pyretis.core.tis` and :py:mod:`pyretis.core.retis`.
4. An analysis is carried out by making use of methods from the
   :py:mod:`analysis <pyretis.analysis>` sub-package.

.. contents:: Table of Contents
   :local:

.. _user-guide-intro-api-core:

The core and forcefield sub-packages
------------------------------------

The :py:mod:`pyretis.core` library defines the core method and
classes and here we will introduce the

* :py:class:`.System` class which defines the system we are investigating.
  This class is actually composed of several other objects, which
  we will also discuss here:

  * :py:class:`.Particles`: A class which represents the particles.

  * :py:class:`.Box`: A class which represents the simulation box.

  * :py:class:`.ForceField`: A class representing the force field.
    (Note: This is defined within :py:mod:`pyretis.forcefield`).

* :py:class:`.Path` class which defines paths/trajectories.

* :py:class:`.PathEnsemble` class which defines path ensembles.

Below is an illustration of how some of these classes are interacting.
The classes (shown as boxes in this figure) will be discussed
more in the following.

.. _figure-relation-base-objects:

.. figure:: /_static/img/api-core-and-forcefield.png
    :width: 50%
    :alt: Illustration of core and force field classes
    :align: center

    Illustration of the relations between classes from the
    core and force field sub-packages. The classes as shown as green
    boxes where the class names are shown with white text.
    Examples of class methods are written as "method()" while
    and class attributes are indicated as just *"attribute"*. Some
    of the attributes are in fact references to other classes and these are
    highlighted with grey boxes and the reference is shown by arrows.
    As can be seen in this figure, the :py:class:`.System` is in practice
    composed of several other classes.


.. _user-guide-intro-api-system:

The System class
^^^^^^^^^^^^^^^^

The :py:class:`.System` class
defines the system we are investigating. It will
typically contain particles, a simulation box and a
force field. This class exposes important parts we
can interact with, in particular, the particles.

Example of creation:

.. code-block:: python

    from pyretis.core import System
    new_system = System(temperature=0.8, units='lj')

This will create an empty system with a set temperature equal to ``0.8`` in
``lj`` units (``lj`` refers to Lennard-Jones :ref:`units <user-guide-units>`).
It is also possible
to specify a box here in case that it needed:

.. code-block:: python

    new_system = System(temperature=0.8, units='lj', box=mybox)

where ``mybox`` can be created as
described :ref:`below <user-guide-intro-api-box>`.

Particles can be added by first creating an object as
described :ref:`below <user-guide-intro-api-particles>`. A short example:

.. code-block:: python

    from pyretis.core import System, Particles
    new_system = System(temperature=0.8, units='lj')
    new_system.particles = Particles()
    new_system.add_particle([0.0, 0.0, 0.0], mass=1.0, name='Ar', ptype=0)

Here, we are setting the :py:attr:`.System.particles` attribute and
using the :py:meth:`.System.add_particle` to add a particle with a given
position, mass, name and type.

.. _user-guide-intro-api-particles:

The Particles class
^^^^^^^^^^^^^^^^^^^

The :py:class:`.Particles` class
represents a collection of particles and in many
ways it can be viewed as a particle list.
Internally in |pyretis|, the particle
list is one of the most important classes. The positions, velocities and
forces are accessed through an instance of this class using
the class attributes :py:attr:`.Particles.pos`, :py:attr:`.Particles.vel` and
:py:attr:`.Particles.force`. Actually, there is an additional particle class
within |pyretis| which is called :py:class:`.ParticlesExt`. This class is
used when |pyretis| is using external engines. It is very similar to the
:py:class:`.Particles` class but it has, in addition, a
:py:attr:`.ParticlesExt.config` which can be used to reference files
which hold the current configuration of the particles.

Here are some examples of interacting with
the :py:class:`.Particles` class,
using :py:meth:`.Particles.add_particle` to add
some particles. Actually, the system will make use of this
method if you are calling :py:meth:`.System.add_particle`.

.. code-block:: python

    from pyretis.core import Particles

    part = Particles(dim=3)

    pos = [0.0, 1.0, 0.0]
    vel = [0.0, 0.0, 0.0]
    force = [0.0, 0.0, 0.0]
    part.add_particle(pos, vel, force, mass=1.0, name='Ar', ptype=0)
    print(part.pos)
    pos = [1.0, 0.0, 0.0]
    part.add_particle(pos, vel, force, mass=1.0, name='Ar', ptype=0)
    print(part.pos)

Here, we can add names to particles using the keyword ``name`` and we
can also specify a particle type using ``ptype``.
The ``name`` can be used to identify/tag
specific particles and is used for output purposes.
Internally, the particle type is more important:
The particle type can
be used to specify parameters for pair interactions which is computed
by the force field.

When we initiate an instance of :py:class:`.Particles`, we define the dimensionality
using the ``dim`` keyword parameter.

.. _user-guide-intro-api-box:

The Box class
^^^^^^^^^^^^^

The :py:class:`.Box` class
defines a simulation box. It is useful in
simulations where we wish to have periodic boundaries. Typically,
we do not interact much with the box beyond creating it.
Boxes are created by passing an (optional) ``cell`` argument which is a list
of floats of form ``[lengthx, lengthy, lengthz]``. If more than
three floats are given, we assume that these represent a flattened version
of the box matrix of the form: ``[xx, yy, zz, xy, xz, yx, yz, xz, zy]``.
At the same time, periodicity can
be specified with the keyword ``periodic`` which is a list of boolean values
that determine if a dimension is periodic or not.  The default is periodic
in all directions.

Some examples:

.. literalinclude:: /_static/api/box_example_create.py
   :language: python
   :lines: 5-12

.. _user-guide-intro-api-forcefield:

The ForceField class
^^^^^^^^^^^^^^^^^^^^

The :py:class:`.ForceField` class
is used to define force fields.
A force field is just a list of functions (and parameters)
which can be used to obtain the force and potential energy.
In general, the force field expect that its constituent potential functions
actually supports calling **three** functions which means that the
potential functions must be slightly more complex than just simple
functions — they need to be classes which
subclass the
:py:class:`.PotentialFunction` class.

If we, for the sake of an example,
let an instance of the :py:class:`.ForceField` class have a constituent potential
function named ``func``, then |pyretis| will assume that it can call:

1. ``func.potential(system)`` to obtain the potential energy.

2. ``func.force(system)`` to obtain the forces and the virial.

3. ``func.potential_and_force(system)`` to obtain the potential energy,
   forces and the virial. Typically, this can be done by just calling
   ``func.potential(system)`` and ``func.force(system)``.

Notice that all these functions should only take in a :py:class:`.System` as
the only parameter.

Let's see an example of how we can set-up a potential function (or class)
and add it to a force field.

First, define the potential function using:

.. literalinclude:: /_static/examples/harmonic1d.py
   :language: python
   :lines: 5-38

Using what we have already discussed above about the :ref:`System <user-guide-intro-api-system>`
and :ref:`Particles <user-guide-intro-api-particles>`
we can make a plot of the potential by adding:

.. literalinclude:: /_static/examples/harmonic1d.py
   :language: python
   :lines: 41-60


.. _user-guide-intro-api-pathe:

The Path and PathEnsemble classes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These two classes are representations of paths and path ensembles.
Paths are essentially trajectories and path ensembles are collections of such
paths. The :py:class:`.PathEnsemble` stores some information about
all paths in the attribute :py:attr:`.PathEnsemble.paths`, but
it will not store the full trajectory (positions, velocities, etc.).
It will, however, keep a reference to the last accepted path in
:py:attr:`.PathEnsemble.last_path`. This :py:class:`.Path` object can,
for instance, be inspected by using the :py:attr:`.Path.phasepoints`
attribute which is a list of system snapshots.

.. _figure-relation-pathensemble:

.. figure:: /_static/img/api-pathensemble.png
    :width: 35%
    :alt: Illustration of the relation PathEnsemble->Path
    :align: center

    Illustration of the relation between the
    :py:class:`.PathEnsemble` and the :py:class:`.Path`.

We give here a short example on how you can interact with a path.
A complicating factor here is that we need to create a random
generator to use with the :py:class:`.Path`. This is done by using
the :py:class:`.RandomGenerator` class and the reader is referred to
the API documentation for information about this.

.. literalinclude:: /_static/examples/pathintro.py
   :language: python
   :lines: 5-38

.. _user-guide-intro-api-simeng:

The simulation and engines sub-packages
---------------------------------------

The :py:mod:`pyretis.simulation` sub-package
defines classes which are used to set-up and define different types
of simulations. Typically, such simulations will need to interact
with and change the state of a given :py:class:`.System`. This
interaction is carried out by a particular engine object which behaves
like :py:class:`.EngineBase`
from the :py:mod:`pyretis.engines` sub-package.
The interaction between these classes are illustrated in the
figure below:

.. _figure-relation-simulation-objects:

.. figure:: /_static/img/api-simulation.png
    :width: 35%
    :alt: Illustration the simulation class.
    :align: center

    Illustration of the relations between the :py:class:`.Simulation`,
    :py:class:`.EngineBase` and :py:class:`.System`. A simulation object will
    typically contain references to a system object and to an engine object. The
    simulation can then use the engine in order to interact with the system. For
    instance can the :py:meth:`.Simulation.step` methods use the
    :py:meth:`.EngineBase.integration_step` in order to update positions, velocities
    and forces in the system.

In this section, we will not give a complete example on how to
create a new simulation class or a new engine class. We refer the
reader to the examples, in particular:

* The section of the user guide which describes
  :ref:`how new external engines <user-guide-engine>` can be interfaced.
* :ref:`The particle swarm optimization example <examples-pso>` in
  which a custom potential and engine is created. Further, in this example
  a customized simulation is set up by making use of :py:meth:`.Simulation.add_task`
* The :ref:`example <examples-external-vvE>` showing how C or Fortran can be
  used to create a custom Velocity Verlet engine.


.. _user-guide-intro-api-simulation:

The Simulation class
^^^^^^^^^^^^^^^^^^^^

The :py:class:`.Simulation` class defines a simulation we can run.
A simulation will typically act on a :py:class:`.System`
and alter its state. We will here just describe the generic
base class :py:class:`.Simulation` and we refer the reader to
the extended :ref:`pyretis API documentation <api-doc>` for
information about specific simulation classes, for instance,
the :py:class:`.SimulationRETIS` class. The most commonly
used methods from the :py:class:`.Simulation` are:

* :py:meth:`run() <.Simulation.run>` which will run a simulation and
  for each step it will yield a dictionary with results.
* :py:meth:`step() <.Simulation.step>` which will run a single step from
  the simulation and return a dictionary with results.
* :py:meth:`is_finished() <.Simulation.is_finished>` which will return ``True``
  if the simulation has ended.
* :py:meth:`add_task() <.Simulation.add_task>` which can be used to add simulation
  tasks to a generic simulation.

Example of creation:

.. code-block:: python

    from pyretis.simulation import Simulation
    new_simulation = Simulation(startcycle=0, endcycle=100)

    for step in new_simulation.run():
        print(step['cycle']['stepno'])


The code block above will create a generic simulation object and run it.
This simulation is not doing anything useful, it will only increment the
current simulation step number from the given ``startcycle`` to the
given ``endcycle``. In order to something more productive, we can attach
tasks to the simulation. This can be done as follows:

.. literalinclude:: /_static/examples/simulationintro.py
   :language: python
   :lines: 5-27

As you can see, a new task is added by defining it as a dictionary:

.. literalinclude:: /_static/examples/simulationintro.py
   :language: python
   :lines: 17-20

The following keywords are used:

* ``func`` which is a reference to a function to use.
* ``args`` which is a list of arguments that should be given
  to the function.
* ``first`` whether the task should be executed on
  the first simulation step (i.e. step 0).
* ``result`` a string which labels the result in the dictionary returned by the
  methods :py:meth:`.Simulation.run()` or :py:meth:`.Simulation.step()`.

Typically, when creating a custom simulation, you will rewrite the
methods :py:meth:`.Simulation.run` and :py:meth:`.Simulation.step` to
fit the custom simulation you are going to perform, rather than adding
tasks. However, for interactive work, short examples etc.,
the :py:meth:`.Simulation.add_task` can be useful.

.. _user-guide-intro-api-engine:

The EngineBase class
^^^^^^^^^^^^^^^^^^^^

The :py:class:`.EngineBase` is a base class defining a generic
engine. In |pyretis| engines are either internal or external. External
engines, e.g. :py:class:`.ExternalMDEngine`,
interfaces external programs while
internal engines, e.g. :py:class:`.MDEngine` interact
directly with a :py:class:`.System` object. Creating an external engine
may be somewhat involved depending on the program you wish to interface.
A description on how to create new external engines can be
found in the :ref:`user guide <user-guide-engine>`.


.. _user-guide-intro-api-orderparameter:

The orderparameter sub-package
------------------------------

The :py:mod:`pyretis.orderparameter` package defines order parameters
to use for path sampling simulation. In |pyretis|, such order parameters
are assumed to be functions of a :py:class:`.System` only. Here,
we give some examples of how a generic order parameter can be
used. For information on how custom order parameters can be
created, we refer to the detailed description in
the :ref:`user guide <user-guide-custom-order>`.

.. _user-guide-intro-api-orderparameterc:

The OrderParameter class
^^^^^^^^^^^^^^^^^^^^^^^^

The most important piece of the :py:class:`.OrderParameter` class
is the actual calculation of the order parameter. This should be
defined in a method like :py:meth:`.OrderParameter.calculate`.
Here, it is assumed that the order parameter takes in an object
like :py:class:`.System`. Since this is described
:ref:`elsewhere <user-guide-custom-order>` we will here just describe
the usage of:

* :py:meth:`.OrderParameter.calculate` which is used to
  calculate the order parameters.
* :py:class:`.CompositeOrderParameter`
  which can be used to combine several collective variables (e.g.
  when you are interested in additional order parameters).

First, we create an order parameter. For simplicity we use the pre-defined
:py:class:`.Position` class:

.. literalinclude:: /_static/examples/orderintro.py
   :language: python
   :lines: 5-8

Next, we can calculate the order parameter as follows
(the system we set up here is just for testing):

.. literalinclude:: /_static/examples/orderintro.py
   :language: python
   :lines: 10-14

Several order parameters can be combined by creating
a :py:class:`.CompositeOrderParameter`. Below is an
example of how this can be used:

.. literalinclude:: /_static/examples/orderintro2.py
   :language: python

.. important:: The first order parameter returned from
   :py:meth:`calculate() <.OrderParameter.calculate>` is taken as the
   progress coordinate used in path sampling simulations.


.. _user-guide-intro-api-analysis:

The analysis sub-package
------------------------

The :py:mod:`pyretis.analysis` sub-package defines tools which
are used in the analysis of simulation output. It defines methods
for running predefined analysis tasks, e.g. :py:func:`.retis_flux`, but
also generic analysis methods, e.g. :py:func:`.running_average`.
Here, we refer the reader to the
:ref:`pyretis API documentation <api-analysis>` for more
information about these methods.

.. _user-guide-intro-api-inout:

The inout sub-package
---------------------

The :py:mod:`pyretis.inout` sub-package contains methods which
|pyretis| make use of in order to interact with you.
These are ways to read input from you and to present output to you.
This sub-package is relatively large and it is in fact made up by
the following sub-packages:

* :py:mod:`pyretis.inout.analysisio` which
  handles the input and output needed for analysis.

* :py:mod:`pyretis.inout.formats` for formatting and
  presenting text-based output.

* :py:mod:`pyretis.inout.plotting` which handles plotting of figures.
  It defines simple things like colors etc.
  for plotting. It also defines functions which can be used for
  specific plotting by the analysis and report tools.

* :py:mod:`pyretis.inout.report` which is used to
  generate reports with results from different simulations.

Again, we refer to the
:ref:`pyretis API documentation <api-inout>` for more
information about these sub-packages.

.. _user-guide-intro-api-initiation:

The initiation sub-package
--------------------------

The :py:mod:`pyretis.initiation` sub-package contains methods to
initialize path ensembles. We refer the reader to the
:ref:`pyretis API documentation <api-initiation>` for more
information about this sub-package.

.. _user-guide-intro-api-setup:

The setup sub-package
---------------------

The setup library contains methods to construct the main
objects of PyRETIS from settings and from restart info.

The :py:mod:`pyretis.setup` which handles creation of objects
from simulation settings.
:ref:`pyretis API documentation <api-setup>` for more
information about this sub-package.

.. _user-guide-intro-api-tools:

The tools sub-package
---------------------

The tools library can be used to generate initial structures for a
simulation. In the tools library the function :py:func:`.generate_lattice`
is defined and it supports the creation of the following lattices where
the shorthand keywords (``sc``, ``sq`` etc.) are used to select a
specific lattice:

- ``sc``: A simple cubic lattice.

- ``sq``: Square lattice (2D) with one atom in the unit cell.

- ``sq2``: Square lattice with two atoms in the unit cell.

- ``bcc``: Body-centered cubic lattice.

- ``fcc``: Face-centered cubic lattice.

- ``hcp``: Hexagonal close-packed lattice.

- ``diamond``: Diamond structure.

A lattice is generated in the following way:

.. code-block:: python

    from pyretis.tools import generate_lattice
    xyz, size = generate_lattice('diamond', [1, 1, 1], lcon=3.57)

Where the first parameter selects the lattice type, the second parameter
selects give the number of repetitions in each direction and the optional
parameter ``lcon`` is the lattice constant. The returned values are
``xyz`` -- the coordinates -- and ``size`` which is the bounding box of the
generated structure. This variable can be used to define a simulation box.
It is also possible to specify a number density:

.. code-block:: python

    from pyretis.tools import generate_lattice
    xyz, size = generate_lattice('diamond', [1, 1, 1], density=1)

If the density is specified, the lattice constant ``lcon`` is deduced:

.. math::

    \text{lcon} = \left(\frac{n}{\rho}\right)^{\frac{1}{d}},

where :math:`n` is the number of particles in the unit cell,
:math:`\rho` the specified number density and :math:`d` the dimensionality.
If we wish to save a generated lattice, this can be done as follows

.. literalinclude:: /_static/api/lattice_example.py
   :language: python
   :lines: 5-15

Some examples of using the |pyretis| library
--------------------------------------------

Here, we show some examples of how we can perform some common tasks
using the |pyretis| library.

Reversing a trajectory
^^^^^^^^^^^^^^^^^^^^^^

|pyretis| will not reverse backward trajectories during a simulation if
you are using an external engine. For visualization purposes, it is very
helpful to reverse these trajectories before viewing them. This can be
accomplished with the |pyretis| library as follows:

* For GROMACS .trr trajectories:

  .. code-block:: python

     from pyretis.inout.formats.gromacs import reverse_trr
     reverse_trr('trajB.trr', 'rev-trajB.trr')

  which will read the trajectory ``trajB.trr`` and store it as ``rev-trajB.trr``.

* For .xyz trajectories:

  .. code-block:: python

     from pyretis.inout.formats.xyz import reverse_xyz_file
     reverse_xyz_file('trajB.xyz', 'rev-trajB.xyz')

  which will read the trajectory ``trajB.xyz`` and store it as ``rev-trajB.xyz``.

Recalculating order parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have an existing trajectory and with to calculate order parameters
for each step in the trajectory, this can be accomplished by making use of
the :ref:`pyretis.tools.recalculate_order <api-tools-recalculate_order>`
module. For example, if you have a trajectory named ``traj.trr`` and
an order parameter defined in an input file ``retis.rst`` this can be
done as follows:

.. code-block:: python

   from pyretis.inout.formats.order import OrderFile
   from pyretis.inout.settings import parse_settings_file
   from pyretis.inout.setup import create_orderparameter
   from pyretis.tools.recalculate_order import recalculate_order
 
   settings = parse_settings_file('retis.rst')
   order_parameter = create_orderparameter(settings)
   options = {'reverse': False, 'maxidx': None, 'minidx': None}
   order = recalculate_order(order_parameter, 'traj.trr', options)
 
   with OrderFile('order.txt', 'w') as outfile:
       outfile.write('# Order parameters')
       for step, data in enumerate(order):
           print(step, data)
           outfile.output(step, data)

This will create a new file ``order.txt`` with the re-calculated order parameters.
The keyword ``reverse`` can be used to tell |pyretis| that you are looking at a
time-reversed trajectory. The keyword ``minidx`` gives the frame number from
which we start calculating (``None`` is equal to reading from the first frame)
and the keyword ``maxidx`` gives the last frame number we will read (setting
it to ``None`` will read until the end of the trajectory file).
