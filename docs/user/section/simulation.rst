.. _user-section-simulation:

The simulation section
======================

The ``simulation`` section defines and selects the simulation |pyretis| will run.

.. pyretis-input-example:: Simulation

   .. code-block:: rst

       Simulation
       ----------
       task = retis
       steps = 20000
       interfaces = [-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, 1.0]

For this section, the keyword ``task`` specifies
the type of simulation to run and this will
also determine the keyword that can be set in
this section. Further, a specific task might also require
*additional sections* to be defined
as detailed in the table below:

.. table:: Possible values for the task keyword
   :class: table-striped table-hover

   +-------------+-------------------------------+----------------------------+-------------------+
   |   Task      | Description                   | Keywords                   | Required sections |
   +=============+===============================+============================+===================+
   |  ``retis``  | A replica exchange transition | ``task``, ``steps``,       | ``tis``,          |
   |             | interface sampling            | ``interfaces``             | ``retis``         |
   |             | simulation.                   |                            |                   |
   +-------------+-------------------------------+----------------------------+-------------------+
   |   ``tis``   | A transition interface        | ``task``, ``steps``,       |                   |
   |             | sampling simulation.          | ``interfaces``,            | ``tis``           |
   |             |                               | ``ensemble_number``,       |                   |
   |             |                               | ``detect``                 |                   |
   +-------------+-------------------------------+----------------------------+-------------------+
   | ``repptis`` | A replica exchange partial    | ``task``, ``steps``,       | ``tis``           |
   |             | path transition interface     | ``interfaces``             |                   |
   |             | sampling simulation.          |                            |                   |
   +-------------+-------------------------------+----------------------------+-------------------+
   | ``pptis``   | A partial path transition     | ``task``, ``steps``,       | ``tis``           |
   |             | interface sampling simulation.| ``interfaces``             |                   |
   +-------------+-------------------------------+----------------------------+-------------------+
   | ``explore`` | A free energy surface         | ``task``, ``steps``,       |                   |
   |             | exploration simulation.       | ``interfaces``             | ``tis``           |
   +-------------+-------------------------------+----------------------------+-------------------+
   | ``md-flux`` | A MD FLUX simulation.         | ``task``, ``steps``,       |                   |
   |             |                               | ``interfaces``             |                   |
   +-------------+-------------------------------+----------------------------+-------------------+
   | ``md-nve``  | A MD NVE simulation.          | ``task``, ``steps``        |                   |
   +-------------+-------------------------------+----------------------------+-------------------+


Since the different tasks may require different keywords, these are described below
individually. In addition, there are some settings which typically are
used in relation to the analysis, or when extending simulations. These are
described in the section on :ref:`common keywords <user-section-simulation-common>`.

.. contents::
   :local:


.. _user-section-simulation-retis:

The retis task
--------------

The ``retis`` task defines a replica exchange transition
interface sampling simulation.

.. pyretis-input-example:: Simulation
   :class-name: task retis

   .. code-block:: rst

       Simulation
       ----------
       task = retis
       steps = 20000
       interfaces = [-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, 1.0]

Keywords for the retis task
^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the ``retis`` task, the following keywords can be set:

.. |simretis_steps| replace:: :ref:`steps <user-section-task-retis-steps>`
.. |simretis_interfaces| replace:: :ref:`interfaces <user-section-task-retis-interfaces>`
.. |simretis_priority| replace:: :ref:`priority_shooting <user-section-task-retis-priority>`

.. _table-simulation-retis-keywords:

.. table:: Keywords for the retis task.
   :class: table-striped table-hover

   +-----------------------+--------------------------------------------------+
   | Keyword               | Description                                      |
   +=======================+==================================================+
   | |simretis_steps|      | The number of total retis steps (cycles)         |
   +-----------------------+--------------------------------------------------+
   | |simretis_interfaces| | The location of the interfaces to consider.      |
   +-----------------------+--------------------------------------------------+
   | |simretis_priority|   | Prioritize the ensembles with less cycles.       |
   +-----------------------+--------------------------------------------------+


In addition, this task requires that the
sections :ref:`TIS <user-section-tis>` and
:ref:`RETIS <user-section-retis>` are defined.

.. _user-section-task-retis-steps:

Keyword steps
.............

.. pyretis-keyword:: steps integer

   The ``steps`` keyword defines the number of RETIS cycles to perform.
   Note that it indicates the goal/final number of the simulation cycles.

   Default
     Not any, this keyword must be specified.

.. _user-section-task-retis-interfaces:

Keyword interfaces
..................

.. pyretis-keyword:: interfaces list of floats

   The ``interfaces`` keyword specifies the interfaces to use in the
   path simulation.

   Default
     Not any, this keyword must be specified.

   .. note::
      ``interfaces`` is topology-defining: it is locked once the
      simulation is started.  A restart that supplies a different list
      aborts with a ``ValueError`` — use
      :ref:`load <user-section-initial-path-load>` or
      :ref:`flexible_restart <user-section-initial-flexible-restart>`
      to change the interface set.

.. _user-section-task-retis-priority:

Keyword priority_shooting
.........................

.. pyretis-keyword:: priority_shooting boolean

   If ``True``, ensembles with lower cycle numbers will be prioritized
   until their cycle numbers equals the others.

   This setting simplify the use of |pyretis| in cluster environments
   in which walltime is rather short. Each ensemble can be, therefore,
   investigated also if launching several runs.

   Default
       The default value is ``priority_shooting = False``.

.. _user-section-simulation-tis:

The tis task
------------

The ``tis`` task defines a transition interface sampling
simulation.

.. pyretis-input-example:: Simulation
   :class-name: task tis (multiple path ensembles)

   .. code-block:: rst

       Simulation
       ----------
       task = tis
       steps = 20000
       interfaces = [-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, 1.0]

.. pyretis-input-example:: Simulation
   :class-name: task tis (a single path ensemble)

   .. code-block:: rst

       Simulation
       ----------
       task = tis
       steps = 20000
       interfaces = [-0.9, -0.9, 1.0]
       detect = -0.8
       ensemble_number = 1

Keywords for the tis task
^^^^^^^^^^^^^^^^^^^^^^^^^

For the ``tis`` task, the following keywords can be set:

.. |simtis_steps| replace:: :ref:`steps <user-section-task-tis-steps>`
.. |simtis_interfaces| replace:: :ref:`interfaces <user-section-task-tis-interfaces>`
.. |simtis_priority| replace:: :ref:`priority_shooting <user-section-task-tis-priority>`
.. |simtis_ensemble_number| replace:: :ref:`ensemble_number <user-section-task-tis-ensemble_number>`
.. |simtis_maxlength| replace:: :ref:`maxlegth <user-section-tis-maxlength>`
.. |simtis_detect| replace:: :ref:`detect <user-section-task-tis-detect>`
 
.. _table-simulation-tis-keywords:

.. table:: Keywords for the tis task.
   :class: table-striped table-hover

   +--------------------------+-----------------------------------------------+
   | Keyword                  | Description                                   |
   +==========================+===============================================+
   | |simtis_steps|           | The total number of steps to perform.         |
   +--------------------------+-----------------------------------------------+
   | |simtis_interfaces|      | The interfaces defining the TIS ensemble.     |
   +--------------------------+-----------------------------------------------+
   | |simtis_priority|        | Prioritize the ensembles with less cycles.    |
   +--------------------------+-----------------------------------------------+
   | |simtis_ensemble_number| | For defining the ensemble number considered.  |
   +--------------------------+-----------------------------------------------+
   | |simtis_detect|          | The interface used for detecting successful   |
   |                          | paths.                                        |
   +--------------------------+-----------------------------------------------+
   | |simtis_maxlength|       | The maximum number of step allowed for a      |
   |                          | paths.                                        |
   +--------------------------+-----------------------------------------------+

In addition, this task requires that the
section :ref:`TIS <user-section-tis>` is defined.

.. _user-section-task-tis-steps:

Keyword steps
.............

.. pyretis-keyword:: steps integer

   The ``steps`` keyword defines the number of TIS cycles to perform.

   Default
     Not any, this keyword must be specified.

.. _user-section-task-tis-interfaces:

Keyword interfaces
..................

.. pyretis-keyword:: interfaces list of floats

   The ``interfaces`` keyword specifies the interfaces to use in the
   path simulation. If the number of interfaces given is 3 or less,
   a single TIS simulation will be performed, otherwise, input files
   for several single TIS simulations will be written.
   These simulations can then be run manually.

   Default
     Not any, this keyword must be specified.

   .. note::
      ``interfaces`` is topology-defining: it is locked once the
      simulation is started.  A restart that supplies a different list
      aborts with a ``ValueError`` — use
      :ref:`load <user-section-initial-path-load>` or
      :ref:`flexible_restart <user-section-initial-flexible-restart>`
      to change the interface set.

.. _user-section-task-tis-priority:

Keyword priority_shooting
.........................

.. pyretis-keyword:: priority_shooting boolean

   If ``True``, ensembles with lower cycle numbers will be prioritized
   until their cycle numbers equals the others.

   This setting simplify the use of |pyretis| in cluster environments
   in which walltime is rather short. Each ensemble can be, therefore,
   investigated also if launching several runs.

   Default
       The default value is ``priority_shooting = False``.

.. _user-section-task-tis-ensemble_number:

Keyword ensemble_number
.......................

.. pyretis-keyword:: ensemble_number integer

   A number specifying the path ensemble we are simulating,
   ``ensemble_number = 1`` corresponds to ``[0^+]``, 2 to ``[1^+]``
   and so on. This is only needed for running a single TIS simulation.

   Default
     Not any, this keyword must be specified.

.. _user-section-task-tis-detect:

Keyword detect
..............

.. pyretis-keyword:: detect float

   A number specifying the interface used for detecting successful
   paths. This is only needed for a single TIS simulation.

   Default
     Not any, this keyword must be specified.

.. _user-section-simulation-repptis:

The repptis task
----------------

The ``repptis`` task defines a replica exchange partial path
transition interface sampling simulation. Paths between adjacent
``pptis`` ensembles can be swapped using the replica exchange move. 

.. pyretis-input-example:: Simulation
   :class-name: task repptis

   .. code-block:: rst

       Simulation
       ----------
       task = repptis
       steps = 20000
       interfaces = [-0.5, -0.3, 0, 0.3, 0.5]

Keywords for the repptis task
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the ``repptis`` task, the following keywords can be set:

.. |simrepptis_steps| replace:: :ref:`steps <user-section-task-repptis-steps>`
.. |simrepptis_interfaces| replace:: :ref:`interfaces <user-section-task-repptis-interfaces>`
.. |simrepptis_priority| replace:: :ref:`priority_shooting <user-section-task-repptis-priority>`

.. _table-simulation-repptis-keywords:

.. table:: Keywords for the repptis task.
   :class: table-striped table-hover

   +-------------------------+-------------------------------------------------+
   | Keyword                 | Description                                     |
   +=========================+=================================================+
   | |simrepptis_steps|      | The number of total repptis steps (cycles)      |
   +-------------------------+-------------------------------------------------+
   | |simrepptis_interfaces| | The location of the interfaces to consider.     |
   +-------------------------+-------------------------------------------------+
   | |simrepptis_priority|   | Prioritize the ensembles with less cycles.      |
   +-------------------------+-------------------------------------------------+

In addition, this task requires that the
sections :ref:`TIS <user-section-tis>` and 
:ref:`RETIS <user-section-retis>` are defined. 

Note that (RE)PPTIS simulations make use of the TIS and RETIS sections in the
input RST file. This is because almost all of the keywords are similar.
Therefore, REPPTIS or PPTIS sections do not exist in the RST file, but re-use
the TIS and RETIS sections.

.. _user-section-task-repptis-steps:

Keyword steps
.............

.. pyretis-keyword:: steps integer

   The ``steps`` keyword defines the number of REPPTIS cycles to perform.
   Note that it indicates the goal/final number of the simulation cycles.

   Default
     Not any, this keyword must be specified.

.. _user-section-task-repptis-interfaces:

Keyword interfaces
..................

.. pyretis-keyword:: interfaces list of floats

   The ``interfaces`` keyword specifies the interfaces to use in the
   path simulation.

   Default
     Not any, this keyword must be specified.

.. _user-section-task-repptis-priority:

Keyword priority_shooting
.........................

.. pyretis-keyword:: priority_shooting boolean

   If ``True``, ensembles with lower cycle numbers will be prioritized
   until their cycle numbers equals the others.

   This setting simplify the use of |pyretis| in cluster environments
   in which walltime is rather short. Each ensemble can be, therefore,
   investigated also if launching several runs.

   Default
       The default value is ``priority_shooting = False``.


.. _user-section-simulation-pptis:

The pptis task
--------------
The ``pptis`` task defines a partial path transition interface sampling
simulation. PPTIS ensembles are shorter than TIS ensembles, cutting the paths'
memory requirements. While for (RE)TIS the paths of ``[i^+]``
need to start from state ``A``, cross interface ``i``, and either commit to 
state ``B`` or return to state ``A``, PPTIS paths of ``[i^+-]`` only need to 
start and end at interfaces ``i-1`` or ``i+1``, having crossed interface ``i``.

.. pyretis-input-example:: Simulation
   :class-name: task pptis

   .. code-block:: rst

       Simulation
       ----------
       task = pptis
       steps = 20000
       interfaces = [-0.5, -0.3, 0, 0.3, 0.5]

Keywords for the pptis task
^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the ``pptis`` task, the following keywords can be set:

.. |simpptis_steps| replace:: :ref:`steps <user-section-task-pptis-steps>`
.. |simpptis_interfaces| replace:: :ref:`interfaces <user-section-task-pptis-interfaces>`
.. |simpptis_priority| replace:: :ref:`priority_shooting <user-section-task-pptis-priority>`


.. _table-simulation-pptis-keywords:

.. table:: Keywords for the pptis task.
   :class: table-striped table-hover

   +--------------------------+-----------------------------------------------+
   | Keyword                  | Description                                   |
   +==========================+===============================================+
   | |simpptis_steps|         | The total number of steps to perform.         |
   +--------------------------+-----------------------------------------------+
   | |simpptis_interfaces|    | The interfaces defining the PPTIS ensemble.   |
   +--------------------------+-----------------------------------------------+
   | |simpptis_priority|      | Prioritize the ensembles with less cycles.    |
   +--------------------------+-----------------------------------------------+

In addition, this task requires that the
section :ref:`TIS <user-section-tis>` is defined. 

Note that ``PPTIS``
simulations use the TIS section, and not the PPTIS section (which does not exist!).

.. _user-section-task-pptis-steps:

Keyword steps
.............

.. pyretis-keyword:: steps integer

   The ``steps`` keyword defines the number of PPTIS cycles to perform.

   Default
     Not any, this keyword must be specified.

.. _user-section-task-pptis-interfaces:

Keyword interfaces
..................

.. pyretis-keyword:: interfaces list of floats

   The ``interfaces`` keyword specifies the interfaces to use in the
   path simulation. If the number of interfaces given is 3 or less,
   a single PPTIS simulation will be performed, otherwise, input files
   for several single PPTIS simulations will be written.
   These simulations can then be run manually.

   Default
     Not any, this keyword must be specified.

.. _user-section-task-pptis-priority:

Keyword priority_shooting
.........................

.. pyretis-keyword:: priority_shooting boolean

   If ``True``, ensembles with lower cycle numbers will be prioritized
   until their cycle numbers equals the others.

   This setting simplify the use of |pyretis| in cluster environments
   in which walltime is rather short. Each ensemble can be, therefore,
   investigated also if launching several runs.

   Default
       The default value is ``priority_shooting = False``.



.. _user-section-simulation-explore:

The explore task
----------------

The ``explore`` task is a molecular dynamics simulation for
exploring the free energy landscape of a transition. It uses
the TIS path simulation approach to allow the exploration
of different regions with different moves.

The idea is to generate paths in order to explore the region delimited
by the first and the last given interface. The generated trajectories
will always be accepted.
To generate more often new paths, with the keyword |simtis_maxlength|,
present in ``tis`` input section, it is possible to limitate the length
of the generated trajectory.

This is a simple strategy to try to locate local minima when exploring
an unknown free energy landscape.


.. _user-section-simulation-mdflux:

The md-flux task
----------------

The ``md-flux`` task is a molecular dynamics simulation for
determining the initial flux for a TIS path simulation.

.. pyretis-input-example:: Simulation
   :class-name: task md-flux

   .. code-block:: rst

       Simulation
       ----------
       task = md-flux
       steps = 10000000
       interfaces = [-0.9]

Keywords for the md-flux task
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the ``md-flux`` task, the following keywords can be specified:

.. |simmdflux_steps| replace:: :ref:`steps <user-section-task-mdflux-steps>`
.. |simmdflux_interfaces| replace:: :ref:`interfaces <user-section-task-mdflux-interfaces>`

.. _table-simulation-mdflux-keywords:

.. table:: Keywords for the md-flux task
   :class: table-striped table-hover

   +------------------------+-------------------------------------------------+
   | Keyword                | Description                                     |
   +========================+=================================================+
   | |simmdflux_steps|      | Defines the number of steps to carry out.       |
   +------------------------+-------------------------------------------------+
   | |simmdflux_interfaces| | Defines the interfaces to consider for the      |
   |                        | flux simulation.                                |
   +------------------------+-------------------------------------------------+


.. _user-section-task-mdflux-steps:

Keyword steps
.............

.. pyretis-keyword:: steps integer

   The ``steps`` keyword specifies the number of MD steps to perform.

   Default
      Not any, this keyword must be specified.

.. _user-section-task-mdflux-interfaces:

Keyword interfaces
..................

.. pyretis-keyword:: interfaces list of floats

   The ``interfaces`` keyword specifies for which interfaces the
   initial flux should be obtained. This can be given as a list
   of floats.

   Default
      Not any, this keyword must be specified.


The md-nve task
---------------

The ``md-nve`` task is a NVE molecular dynamics simulation.

.. pyretis-input-example:: Simulation
   :class-name: task md-nve

   .. code-block:: rst

       Simulation
       ----------
       task = md-nve
       steps = 10000000

Keywords for the md-nve task
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following keywords can be specified for the ``md-nve`` task:

.. |simmdnve_steps| replace:: :ref:`steps <user-section-task-mdnve-steps>`

.. _table-simulation-mdnve-keywords:

.. table:: Keywords for the md-nve task
   :class: table-striped table-hover

   +------------------+-------------------------------------------------------+
   | Keyword          | Description                                           |
   +==================+=======================================================+
   | |simmdnve_steps| | The number of overall steps to consider               |
   |                  | for the simulation.                                   |
   +------------------+-------------------------------------------------------+


.. _user-section-task-mdnve-steps:

Keyword steps
.............

.. pyretis-keyword:: steps integer

   The ``steps`` keyword specifies the number of MD steps to perform.

   Default
      Not any, this keyword must be specified.


.. _user-section-simulation-common:

Common keywords
---------------

The following keywords are common to all simulation tasks:

.. |endcycle| replace:: :ref:`endcycle <user-section-simulation-endcycle>`
.. |startcycle| replace:: :ref:`startcycle <user-section-simulation-startcycle>`
.. |exepath| replace:: :ref:`exe-path <user-section-simulation-exe-path>`
.. |restart| replace:: :ref:`restart <user-section-simulation-restart>`
.. |remove_gen| replace:: :ref:`remove_generate <user-section-simulation-remove_gen>`

.. _table-simulation-common-keywords:

.. table:: Common keywords
   :class: table-striped table-hover

   +------------------------+-------------------------------------------------+
   | Keyword                | Description                                     |
   +========================+=================================================+
   | |endcycle|             | Specifies the cycle step the simulation ended.  |
   +------------------------+-------------------------------------------------+
   | |exepath|              | Specifies the directory from where the          |
   |                        | simulation was executed.                        |
   +------------------------+-------------------------------------------------+
   | |restart|              | Specifies the restart file to use.              |
   +------------------------+-------------------------------------------------+
   | |startcycle|           | Specifies the cycle step the simulation should  |
   |                        | start at.                                       |
   +------------------------+-------------------------------------------------+
   | |remove_gen|           | Specifies if the generated files should be      |
   |                        | removed after every ensemble TIS move.          |
   +------------------------+-------------------------------------------------+


.. _user-section-simulation-endcycle:

Keyword endcycle
^^^^^^^^^^^^^^^^

.. pyretis-keyword:: endcycle integer

   The ``endcycle`` keyword specifies the cycle number at which the
   simulation ended. If the simulation was stopped before the specified
   number of steps were reached, the ``endcycle`` will be difference
   from the ``steps`` keyword. This keyword will be set and updated
   by |pyretis| and written to the processed input file.

   Default
      Not any, note that this keyword will be updated and modified by
      |pyretis| as part of the output process. This keyword is **only**
      used by the analysis program.


.. _user-section-simulation-exe-path:

Keyword exe-path
^^^^^^^^^^^^^^^^

.. pyretis-keyword:: exe-path string

   The ``exe-path`` keyword specifies the location from where the simulation
   was executed. This will be set and updated by |pyretis|.


   Default
      Not any, note that this keyword will be updated and modified by
      |pyretis| as part of the output process. This keyword **does not**
      need to be set in the input file by the user since |pyretis| will
      update this setting automatically. 


.. _user-section-simulation-restart:

Keyword restart
^^^^^^^^^^^^^^^

.. pyretis-keyword:: restart string

   The ``restart`` keyword specifies the path to the restart file to use
   for restarting/continuing the selected simulation task. To be used,
   it requires that ``method`` keyword in the ``initial-path`` section
   is set to ``restart``.

   Default
      The default value is ``pyretis.restart``.


.. _user-section-simulation-startcycle:

Keyword startcycle
^^^^^^^^^^^^^^^^^^

.. pyretis-keyword:: startcycle integer

   The ``startcycle`` keyword specifies the cycle number the simulation
   starts at. This can be used, for instance, when extending a simulation
   to tell |pyretis| that the simulation should start at a specified
   step number.

   Default
      The default value is ``0``.

.. _user-section-simulation-remove_gen:

Keyword remove_generate
^^^^^^^^^^^^^^^^^^^^^^^

.. pyretis-keyword:: remove_generate boolean

   The ``remove_generate`` keyword specifies whether the generated files should 
   be removed after every ensemble TIS move. This is useful when running TIS
   simulations in a cluster environment, where the generated files can
   take up a lot of space. If ``remove_generate = True``, the generated files
   will be removed after every ensemble TIS move. If ``remove_generate = False``,
   the generated files will only be removed after an entire cycle has finished.

   Default
      The default value is ``True``.
