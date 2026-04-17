.. _user-section-engine:

The Engine section
==================

The ``engine`` section specifies the engine to use for the dynamics.

.. pyretis-input-example:: engine

   .. code-block:: rst

      Engine
      ------
      class = VelocityVerlet
      timestep = 0.002

In |pyretis|, the engines are used as numerical integrators for
Newton's equations of motion. The different
engines need different settings
and this is specified in detail below for the
different types:


.. |engine_verlet| replace:: :ref:`Verlet <user-section-engine-verlet>`

.. |engine_langevin| replace:: :ref:`Langevin <user-section-engine-langevin>`

.. |engine_vverlet| replace:: :ref:`Velocity Verlet <user-section-engine-velocity-verlet>`

.. |engine_gromacs| replace:: :ref:`GROMACS <user-section-engine-gromacs>`

.. |engine_cp2k| replace:: :ref:`CP2K <user-section-engine-cp2k>`

.. |engine_lammps| replace:: :ref:`LAMMPS <user-section-engine-lammps>`

.. |engine_user| replace:: :ref:`User-defined <user-section-engine-user-defined>`

.. _table-engine-types:

.. table:: Supported engines in |pyretis|
   :class: table-striped table-hover

   +-------------------+------------------------------------------------------+
   | Engine            | Description                                          |
   +===================+======================================================+
   | |engine_verlet|   | Internal engine, integrate using the Verlet scheme.  |
   +-------------------+------------------------------------------------------+
   | |engine_langevin| | Internal engine, stochastic dynamics.                |
   +-------------------+------------------------------------------------------+
   | |engine_vverlet|  | Internal engine, integration using the Velocity      |
   |                   | verlet scheme.                                       |
   +-------------------+------------------------------------------------------+
   | |engine_gromacs|  | External engine, using GROMACS.                      |
   +-------------------+------------------------------------------------------+
   | |engine_cp2k|     | External engine, using CP2K.                         |
   +-------------------+------------------------------------------------------+
   | |engine_lammps|   | External engine, using LAMMPS.                       |
   +-------------------+------------------------------------------------------+
   | |engine_user|     | Internal/External engine, user-defined.              |
   +-------------------+------------------------------------------------------+


.. _user-section-engine-verlet:

The Verlet engine
-----------------

The Verlet engine integrates the equations of motion
according to the Verlet scheme.
NB: Only MD is allowed for this scheme. The engine is selected by
specifying the class ``Verlet`` and the time step:


.. pyretis-input-example:: Engine
   :class-name: Verlet

   .. code-block:: rst

      Engine
      ------
      class = Verlet # select Verlet integrator
      timestep = 0.002  # time step for the integration

Keywords for the Verlet engine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the Verlet engine, the following keywords can be set:

.. |verlet_class| replace:: :ref:`class <user-section-engine-verlet-class>`

.. |verlet_time| replace:: :ref:`timestep <user-section-engine-verlet-time>`

.. _table-verlet-keywords:

.. table:: Keywords for the Verlet engine.
   :class: table-striped table-hover

   +----------------+---------------------------------------------------+
   | Keyword        | Description                                       |
   +================+===================================================+
   | |verlet_class| | Selects the Verlet engine.                        |
   +----------------+---------------------------------------------------+
   | |verlet_time|  | Defines the time step.                            |
   +----------------+---------------------------------------------------+


.. _user-section-engine-verlet-class:

Keyword class
.............

.. pyretis-keyword:: class Verlet
   :specific: yes

   The class selects the Verlet engine and it should be set to ``Verlet``.

.. _user-section-engine-verlet-time:

Keyword timestep
................

.. pyretis-keyword:: timestep float

   The ``timestep`` keyword defines the time step for the
   engine in **internal units**.

   Default:
       Not any. This keyword must be specified.


.. _user-section-engine-velocity-verlet:

The Velocity Verlet engine
--------------------------

The Velocity Verlet engine integrates the equations of motion
according to the Velocity Verlet scheme.
The engine is selected by specifying the
class ``VelocityVerlet`` and the time step:


.. pyretis-input-example:: Engine
   :class-name: VelocityVerlet

   .. code-block:: rst

      Engine
      ------
      class = VelocityVerlet # select Velocity Verlet integrator
      timestep = 0.002  # time step for the integration


Keywords for the Velocity Verlet engine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the Velocity Verlet engine, the following keywords can be set:

.. |vverlet_class| replace:: :ref:`class <user-section-engine-vverlet-class>`

.. |vverlet_time| replace:: :ref:`timestep <user-section-engine-vverlet-time>`

.. _table-vverlet-keywords:

.. table:: Keywords for the Velocity Verlet engine.
   :class: table-striped table-hover

   +-----------------+---------------------------------------------------+
   | Keyword         | Description                                       |
   +=================+===================================================+
   | |vverlet_class| | Selects the Velocity Verlet engine.               |
   +-----------------+---------------------------------------------------+
   | |vverlet_time|  | Defines the time step.                            |
   +-----------------+---------------------------------------------------+


.. _user-section-engine-vverlet-class:

Keyword class
.............

.. pyretis-keyword:: class VelocityVerlet
   :specific: yes

   The class selects the Velocity Verlet engine and it should be set to ``VelocityVerlet``.

.. _user-section-engine-vverlet-time:

Keyword timestep
................

.. pyretis-keyword:: timestep float

   The ``timestep`` keyword defines the time step for the
   engine in **internal units**.

   Default:
       Not any. This keyword must be specified.


.. _user-section-engine-langevin:

The Langevin engine
-------------------

The engine is selected by specifying the
class ``Langevin``.
This is a stochastic (Brownian) integrator and a description of
the implementation can be found in e.g. [1]_. The integrator
is fully specified as follows:

.. pyretis-input-example:: Engine
   :class-name: Langevin

   .. code-block:: rst

       Engine
       ------
       class = Langevin # select Langevin integrator
       timestep = 0.002  # time step for the integration
       gamma = 0.3  # set gamma value,
       seed = 0  # set seed for random value generator used
       high_friction = False  # are we in the high friction limit?

Keywords for the Langevin engine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the Langevin engine, the following keywords can be set:

.. |langevin_class| replace:: :ref:`class <user-section-engine-langevin-class>`

.. |langevin_time| replace:: :ref:`timestep <user-section-engine-langevin-time>`

.. |langevin_gamma| replace:: :ref:`gamma <user-section-engine-langevin-gamma>`

.. |langevin_seed| replace:: :ref:`seed <user-section-engine-langevin-seed>`

.. |langevin_hf| replace:: :ref:`high_friction <user-section-engine-langevin-hf>`


.. _table-langevin-keywords:

.. table:: Keywords for the Langevin engine.
   :class: table-striped table-hover

   +------------------+-------------------------------------------------------+
   | Keyword          | Description                                           |
   +==================+=======================================================+
   | |langevin_class| | Selects the Langevin engine.                          |
   +------------------+-------------------------------------------------------+
   | |langevin_time|  | Which defines the time step.                          |
   +------------------+-------------------------------------------------------+
   | |langevin_gamma| | Which defines the the friction parameter.             |
   +------------------+-------------------------------------------------------+
   | |langevin_seed|  | Which defines a seed to use for the random number     |
   |                  | generator.                                            |
   +------------------+-------------------------------------------------------+
   | |langevin_hf|    | Which determines if we are in the high friction       |
   |                  | limit or not.                                         |
   +------------------+-------------------------------------------------------+


.. _user-section-engine-langevin-class:

Keyword class
.............

.. pyretis-keyword:: class Langevin
   :specific: yes

   The class selects the Langevin engine and it should be set to ``Langevin``.

.. _user-section-engine-langevin-time:

Keyword timestep
................

.. pyretis-keyword:: timestep float

   The ``timestep`` keyword defines the time step for the
   engine in **internal units**.

   Default:
       Not any. This keyword must be specified.

.. _user-section-engine-langevin-gamma:

Keyword gamma
.............

.. pyretis-keyword:: gamma float

   This is the gamma parameter (:math:`\gamma`,
   see the
   :ref:`integration scheme <user-section-engine-langevin-hf>`)
   for the Langevin integrator.

   Default:
     Not any. This keyword must be specified.

.. _user-section-engine-langevin-seed:

Keyword seed
............

.. pyretis-keyword:: seed integer

   The seed value is an integer that is used to seed the random
   number generator used for the stochastic integration.

   Default:
     The default seed is zero: ``seed = 0``

.. _user-section-engine-langevin-hf:

Keyword high_friction
.....................

.. pyretis-keyword:: high_friction boolean

   This keyword determines how the equations of motion
   are integrated. This is a boolean value and thus it can
   be set to ``True`` or ``False``:

   * If ``high_friction`` is ``True``, we are in the high friction
     limit and the equations of
     motion are integrated according to

     .. math::

        \mathbf{r}(t + \Delta t) = \mathbf{r}(t) + \gamma \Delta t \mathbf{f}(t)/m + \delta \mathbf{r},

     where :math:`\mathbf{f}(t)` is the deterministic force and the
     displacements (:math:`\delta \mathbf{r}`) are drawn from a
     normal distribution,


   * If ``high_friction`` is ``False``, we are in the low friction
     limit and the equations of motion are integrated according to

     .. math::

        \mathbf{r}(t + \Delta t) = \mathbf{r}(t) + c_1 \Delta t  \mathbf{v}(t) + c_2 \mathbf{a}(t) \Delta t^2 + \delta \mathbf{r},

     where,

     .. math::

        \mathbf{v}(r + \Delta t) = c_0 \mathbf{v}(t) + (c_1-c_2) \Delta t \mathbf{a}(t) + c_2 \Delta t \mathbf{a}(t+\Delta t) + \delta \mathbf{v},

     and :math:`c_0 = \text{e}^{-\gamma \Delta t}`,
     :math:`c_1 = (1 - c_0) / ( \gamma \Delta t)` and
     :math:`c_2 = (1 - c_1) / ( \gamma \Delta t)`. In this case,
     :math:`\delta \mathbf{r}` and :math:`\delta \mathbf{v}` are
     obtained as stochastic variables.

   Default:
     The default setting is ``high_friction = False``


.. _user-section-engine-gromacs:

The GROMACS engine
------------------

The GROMACS engine will use GROMACS [2]_ in order
to integrate the equations of motion.
The engine is selected by specifying the
class ``gromacs`` with some additional keywords:

.. pyretis-input-example:: Engine
   :class-name: gromacs

   .. code-block:: rst

      Engine settings
      ---------------
      class = gromacs
      gmx = gmx
      mdrun = gmx mdrun
      input_path = gromacs_input
      timestep = 0.002
      subcycles = 3
      maxwarn = 0
      write_vel = True
      write_force = False
      gmx_format = g96

Currently, there is also an experimental ``gromacs2`` engine
implemented in |pyretis|. This engine can read the output from GROMACS
on the fly, but it has not been extensively tested. It will however
run faster than the ``gromacs`` engine.


Keywords for the GROMACS engine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the GROMACS engine, the following keywords can be set:

.. |gmx_class| replace:: :ref:`class <user-section-engine-gromacs-class>`

.. |gmx_gmx| replace:: :ref:`gmx <user-section-engine-gromacs-gmx>`

.. |gmx_mdrun| replace:: :ref:`mdrun <user-section-engine-gromacs-mdrun>`

.. |gmx_input| replace:: :ref:`input_path <user-section-engine-gromacs-input>`

.. |gmx_time| replace:: :ref:`timestep <user-section-engine-gromacs-time>`

.. |gmx_sub| replace:: :ref:`subcycles <user-section-engine-gromacs-sub>`

.. |gmx_maxwarn| replace:: :ref:`maxwarn <user-section-engine-gromacs-maxwarn>`

.. |gmx_write_vel| replace:: :ref:`write_vel <user-section-engine-gromacs-write_vel>`

.. |gmx_write_force| replace:: :ref:`write_force <user-section-engine-gromacs-write_force>`

.. |gmx_gmx_format| replace:: :ref:`gmx_format <user-section-engine-gromacs-gmx_format>`

.. |gmx_domain_dec| replace:: :ref:`domain_decomp <user-section-engine-gromacs-domain_decomp>`

.. _table-gromacs-keywords:

.. table:: Keywords for the GROMACS engine.
   :class: table-striped table-hover

   +-------------------+------------------------------------------------------+
   | Keyword           | Description                                          |
   +===================+======================================================+
   | |gmx_class|       | Selects the GROMACS engine.                          |
   +-------------------+------------------------------------------------------+
   | |gmx_gmx|         | Defines the command used to execute the GROMACS      |
   |                   | ``gmx`` program.                                     |
   +-------------------+------------------------------------------------------+
   | |gmx_mdrun|       | Defines the command used to execute the GROMACS      |
   |                   | ``mdrun`` program.                                   |
   +-------------------+------------------------------------------------------+
   | |gmx_input|       | Defines the folder containing the input to GROMACS.  |
   +-------------------+------------------------------------------------------+
   | |gmx_time|        | Which defines the time step.                         |
   +-------------------+------------------------------------------------------+
   | |gmx_sub|         | Which defines the number of steps GROMACS will       |
   |                   | execute before we calculate order parameters, i.e.   |
   |                   | the frequency of output from GROMACS.                |
   +-------------------+------------------------------------------------------+
   | |gmx_maxwarn|     | Which sets the ``-maxwarn`` option of GROMACS        |
   |                   | grompp.                                              |
   +-------------------+------------------------------------------------------+
   | |gmx_write_vel|   | Determines if velocities should be written by        |
   |                   | GROMACS or not.                                      |
   +-------------------+------------------------------------------------------+
   | |gmx_write_force| | Determines if forces should be written by GROMACS or |
   |                   | not.                                                 |
   +-------------------+------------------------------------------------------+
   | |gmx_gmx_format|  | Can be used to select the format used for GROMACS    |
   |                   | configurations.                                      |
   +-------------------+------------------------------------------------------+
   | |gmx_domain_dec|  | Can be toggled on if domain decomposition is used    |
   |                   | in GROMACS.                                          |
   +-------------------+------------------------------------------------------+


.. _user-section-engine-gromacs-class:

Keyword class
.............

.. pyretis-keyword:: class gromacs
   :specific: yes

   The class selects the GROMACS engine and it should be set to ``gromacs``.

.. _user-section-engine-gromacs-gmx:

Keyword gmx
...........

.. pyretis-keyword:: gmx string

   This keyword sets the command |pyretis| will use for executing the GROMACS
   ``gmx`` program. This can be used if you have different versions
   of GROMACS installed, for instance, a single precision build named ``gmx`` and
   a double precision build named ``gmx_d``.

   Default:
       Not any. This keyword must be specified.

.. _user-section-engine-gromacs-mdrun:

Keyword mdrun
.............

.. pyretis-keyword:: mdrun string

   This keyword sets the command |pyretis| will use for executing the GROMACS
   ``mdrun`` program. This can be used to execute an MPI compiled
   version of ``mdrun``, by, for instance,
   setting: ``mdrun = mpiexec_mpt mdrun_mpi``. Note that this command is
   specific to the system you are using.

   Default:
       Not any. This keyword must be specified.

.. _user-section-engine-gromacs-input:

Keyword input_path
..................

.. pyretis-keyword:: input_path string

   This keyword sets the directory where |pyretis| will look for input files to
   use with GROMACS. If, for instance, ``input_path = gromacs-input`` is set, then
   |pyretis| will look in the folder ``gromacs-input`` relative to the directory
   |pyretis| is executed in.

   Further, in the folder specified by ``input_path`` the following files **must**
   be present:

   - ``conf.g96``: The initial configuration.
   - ``grompp.mdp``: The simulation settings to use. Note that |pyretis| will alter the
     frequency of output and time step to match the setting given in the input to
     |pyretis| using the :ref:`timestep <user-section-engine-gromacs-time>`
     and :ref:`subcycles <user-section-engine-gromacs-sub>` keywords. Effectively,
     |pyretis| will use ``grompp.mdp`` as a template and create a ``pyretis.mdp``, in
     the same folder, which is used for the GROMACS simulations.
   - ``topol.top``: The topology for your system.


   Default:
       Not any. This keyword must be specified.

.. _user-section-engine-gromacs-time:

Keyword timestep
................

.. pyretis-keyword:: timestep float

   The ``timestep`` keyword defines the time step for the
   engine in **internal units**, which in this case will
   be :ref:`gromacs units <table_unit_systems_time>`.

   Default:
       Not any. This keyword must be specified.

.. _user-section-engine-gromacs-sub:

Keyword subcycles
.................

.. pyretis-keyword:: subcycles integer

   The ``subcycles`` keyword defines the frequency of output from GROMACS
   and thus the frequency by which the order parameter(s) are obtained.

   Default:
       Not any. This keyword must be specified.

.. _user-section-engine-gromacs-maxwarn:

Keyword maxwarn
...............

.. pyretis-keyword:: maxwarn integer

   The ``maxwarn`` keyword set the maximum number of warnings
   the GROMACS grompp command will ignore.

   Default:
       The default value is ``maxwarn = 0``.

.. _user-section-engine-gromacs-write_vel:

Keyword write_vel
.................

.. pyretis-keyword:: write_vel boolean

   The ``write_vel`` keyword determines if GROMACS should write velocities when
   integrating the equations of motion. This can safely be set to ``False`` if
   your order parameter does not depend on the velocities.

   Default:
       The default value is ``write_vel = True``.

.. _user-section-engine-gromacs-write_force:

Keyword write_force
...................

.. pyretis-keyword:: write_force  boolean

   The ``write_force`` keyword determines if GROMACS should write forces when
   integrating the equations of motion. Typically, the forces are not needed
   unless your order parameter depends on them.

   Default:
       The default value is ``write_force = False``.

.. _user-section-engine-gromacs-gmx_format:

Keyword gmx_format
..................

.. pyretis-keyword:: gmx_format string

   The ``gmx_format`` keyword specifies the format to use for GROMACS
   configurations. It can be set to ``gro`` or ``g96``.

   Default:
       The default value is ``gmx_format = gro``.

.. _user-section-engine-gromacs-domain_decomp:

Keyword domain_decomp
.....................

.. pyretis-keyword:: domain_decomp boolean

   The ``domain_decomp`` keyword will execute ``gmx trjconv -pbc whole ...`` 
   to a frame from which the system is propagated. This was implemented as 
   some GROMACS versions did not handle bonds over periodic boundaries correctly.

   Default:
       The default value is ``domain_decomp = False``.


.. _user-section-engine-cp2k:

The ``CP2K`` engine
-------------------

The CP2K engine will use CP2K [3]_ in order
to integrate the equations of motion.
The engine is selected by specifying the
class ``cp2k`` with some additional keywords:

.. pyretis-input-example:: Engine
   :class-name: cp2k

   .. code-block:: rst

      Engine
      ------
      class = cp2k
      cp2k = cp2k
      input_path = cp2k-input
      timestep = 0.002
      subcycles = 5


Keywords for the CP2K engine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the CP2K engine, the following keywords can be set:

.. |cp2k_class| replace:: :ref:`class <user-section-engine-cp2k-class>`

.. |cp2k_cp2k| replace:: :ref:`cp2k <user-section-engine-cp2k-cp2k>`

.. |cp2k_input| replace:: :ref:`input_path <user-section-engine-cp2k-input>`

.. |cp2k_time| replace:: :ref:`timestep <user-section-engine-cp2k-time>`

.. |cp2k_sub| replace:: :ref:`subcycles <user-section-engine-cp2k-sub>`


.. _table-cp2k-keywords:

.. table:: Keywords for the CP2K engine.
   :class: table-striped table-hover

   +--------------+-----------------------------------------------------------+
   | Keyword      | Description                                               |
   +==============+===========================================================+
   | |cp2k_class| | Selects the CP2K engine.                                  |
   +--------------+-----------------------------------------------------------+
   | |cp2k_cp2k|  | Defines the command used to execute CP2K.                 |
   +--------------+-----------------------------------------------------------+
   | |cp2k_input| | Defines the location of input files to use for CP2K.      |
   +--------------+-----------------------------------------------------------+
   | |cp2k_time|  | Which defines the time step.                              |
   +--------------+-----------------------------------------------------------+
   | |cp2k_sub|   | Which defines the number of steps CP2K will execute       |
   |              | before we calculate order parameter(s).                   |
   +--------------+-----------------------------------------------------------+


.. _user-section-engine-cp2k-class:

Keyword class
.............

.. pyretis-keyword:: class cp2k
   :specific: yes

   The class selects the CP2K engine and it should be set to ``cp2k``.

.. _user-section-engine-cp2k-cp2k:

Keyword cp2k
............

.. pyretis-keyword:: cp2k string

   The command used for executing CP2K.

   Default:
       Not any. This keyword must be specified.

.. _user-section-engine-cp2k-input:

Keyword input_path
..................

.. pyretis-keyword:: input_path string

   This keyword sets the directory where |pyretis| will look for input files to
   use with CP2K. If for instance ``input_path = cp2k-input`` is set, then
   |pyretis| will look in the folder ``cp2k-input`` relative to the directory
   |pyretis| is executed in.

   Further, in the folder specified by ``input_path`` the following files **must**
   be present:

   * ``initial.xyz``: A file with the initial configuration.

   * ``cp2k.inp``: A CP2K input file, defining a MD simulation.

   Default:
       Not any. This keyword must be specified.

.. _user-section-engine-cp2k-time:

Keyword timestep
................

.. pyretis-keyword:: timestep float

   The ``timestep`` keyword defines the time step for the
   engine in **internal units**.

   Default:
       Not any. This keyword must be specified.

.. _user-section-engine-cp2k-sub:

Keyword subcycles
.................

.. pyretis-keyword:: subcycles integer

   The ``subcycles`` defines the frequency of the output from CP2K
   and thus the frequency by which the order parameter(s) are obtained.

   Default:
       Not any. This keyword must be specified.



.. _user-section-engine-lammps:

The ``LAMMPS`` engine
---------------------

The LAMMPS engine will use LAMMPS [4]_ in order
to integrate the equations of motion.
A description of how to use the LAMMPS engine with |pyretis| can
be found in :ref:`the user guide on LAMMPS <user-guide-lammps>`.

The engine is selected by specifying the
class ``lammps`` with some additional keywords:

.. pyretis-input-example:: Engine
   :class-name: lammps

   .. code-block:: rst

      Engine
      ------
      class = lammps
      lmp = lmp
      input_path = lammps_input
      subcycles = 1
      extra_files = ['potential.in']


Keywords for the LAMMPS engine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the LAMMPS engine, the following keywords can be set:

.. |lammps_class| replace:: :ref:`class <user-section-engine-lammps-class>`

.. |lammps_lmp| replace:: :ref:`lmp <user-section-engine-lammps-lmp>`

.. |lammps_input| replace:: :ref:`input_path <user-section-engine-lammps-input>`

.. |lammps_sub| replace:: :ref:`subcycles <user-section-engine-lammps-sub>`

.. |lammps_extra| replace:: :ref:`extra_files <user-section-engine-lammps-extra>`


.. _table-lammps-keywords:

.. table:: Keywords for the LAMMPS engine.
   :class: table-striped table-hover

   +----------------+-----------------------------------------------------------+
   | Keyword        | Description                                               |
   +================+===========================================================+
   | |lammps_class| | Selects the LAMMPS engine.                                |
   +----------------+-----------------------------------------------------------+
   | |lammps_lmp|   | Defines the command used to execute LAMMPS.               |
   +----------------+-----------------------------------------------------------+
   | |lammps_input| | Defines the location of input files to use for LAMMPS.    |
   +----------------+-----------------------------------------------------------+
   | |lammps_sub|   | Which defines the number of steps LAMMPS will execute     |
   |                | before we calculate order parameter(s).                   |
   +----------------+-----------------------------------------------------------+
   | |lammps_extra| | Defines the location of input files to use for LAMMPS.    |
   +----------------+-----------------------------------------------------------+


.. _user-section-engine-lammps-class:

Keyword class
.............

.. pyretis-keyword:: class lammps
   :specific: yes

   The class selects the LAMMPS engine and it should be set to ``lammps``.

.. _user-section-engine-lammps-lmp:

Keyword lmp
...........

.. pyretis-keyword:: lmp string

   The command used for executing LAMMPS.

   Default:
       Not any. This keyword must be specified.

.. _user-section-engine-lammps-input:

Keyword input_path
..................

.. pyretis-keyword:: input_path string

   This keyword sets the directory where |pyretis| will look for input files to
   use with LAMMPS. If for instance ``input_path = lammps_input`` is set, then
   |pyretis| will look in the folder ``lammps_input`` relative to the directory
   |pyretis| is executed in.

   Further, in the folder specified by ``input_path`` the following files **must**
   be present:

   * ``system.data``: A LAMMPS data file with the initial configuration.

   * ``lammps.in``: A LAMMPS input file, defining a MD simulation.

   * ``order.in``: A LAMMPS input file, defining the order parameter calculation.

   Default:
       Not any. This keyword must be specified.

.. _user-section-engine-lammps-sub:

Keyword subcycles
.................

.. pyretis-keyword:: subcycles integer

   The ``subcycles`` defines the frequency of the output from LAMMPS,
   and thus the frequency by which the order parameter(s) are obtained.

   Default:
       Not any. This keyword must be specified.

.. _user-section-engine-lammps-extra:

Keyword extra
.............

.. pyretis-keyword:: extra_files list of strings

   The ``extra_files`` defines additional files that might be needed
   to execute LAMMPS. For instance, may the user refer to several
   additional LAMMPS script in the given ``lammps.in`` file. This keyword
   is to make |pyretis| aware of these additional files.

   Default:
       None.




.. _user-section-engine-user-defined:

User-defined engines
--------------------

It is also possible to extend |pyretis| with user-defined
engines, written in for instance Python, FORTRAN or C.
User-defined engines are specified in Python modules that
|pyretis| can load and they are typically specified as
follows:

.. _user-section-engine-user-example:

.. pyretis-input-example:: Engine
   :class-name: a user-defined engine

   .. code-block:: rst

       Engine
       ------
       class = VelocityVerletF
       module = vvengineF.py
       timestep = 0.002


Keywords for user-defined engines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Typically, the user-defined engines requires that, at least,
the following keywords are set:

.. |custom_class| replace:: :ref:`class <user-section-engine-user-class>`
.. |custom_module| replace:: :ref:`module <user-section-engine-user-module>`
.. |custom_time| replace:: :ref:`timestep <user-section-engine-user-time>`

.. _table-keyword-user-engine:

.. table:: Keywords for user-defined engines.
   :class: table-striped table-hover

   +-----------------+--------------------------------------------------------+
   | Keyword         | Description                                            |
   +=================+========================================================+
   | |custom_class|  | Selects the engine.                                    |
   +-----------------+--------------------------------------------------------+
   | |custom_module| | Defines the Python file where the engine class is      |
   |                 | defined.                                               |
   +-----------------+--------------------------------------------------------+
   | |custom_time|   | Defines the timestep for the engine.                   |
   +-----------------+--------------------------------------------------------+


.. _user-section-engine-user-class:

Keyword class
.............

.. pyretis-keyword:: class string

   This keyword selects the engine and it should be set to
   the class name as it is defined in the given
   :ref:`module <user-section-engine-user-module>`.

.. _user-section-engine-user-module:

Keyword module
..............

.. pyretis-keyword:: module string

   This keyword specified the location of the file containing the
   user-defined class. This file must be accessible by |pyretis|.

   Default
       Not any. This keyword must be specified.

.. _user-section-engine-user-time:

Keyword timestep
................

.. pyretis-keyword:: timestep float

   The ``timestep`` keyword defines the time step for the
   engine in **internal units**.

   Default:
       Not any. This keyword must be specified.


Additional keywords
...................

In addition, user-defined keywords can be specified, e.g.:

.. pyretis-input-example:: Engine
   :class-name: a user-defined engine

   .. code-block:: rst

       Engine
       ------
       class = VelocityVerletF
       module = vvintegratorf.py
       timestep = 0.002
       keyword = 100
       otherkeyword = 31


References
----------

.. [1] M. P. Allen and D. J. Tildesley,
       Computer Simulation of Liquids,
       1989, Oxford University Press.

.. [2] http://www.gromacs.org

.. [3] https://www.cp2k.org/

.. [4] https://lammps.sandia.gov/
