.. _user-section-output:

The Output section
==================

The ``output`` section defines settings for how |pyretis| should create
output files.

.. pyretis-input-example:: Output

   .. code-block:: rst

      Output settings
      ---------------
      backup = 'overwrite'
      pathensemble-file = 1
      order-file = 1
      energy-file = 1
      screen = 10

The following keywords can be set for the output section:

.. |output_backup| replace:: :ref:`bakcup <user-section-output-backup>`
.. |output_prefix| replace:: :ref:`prefix <user-section-output-prefix>`
.. |output_energy| replace:: :ref:`energy-file <user-section-output-keyword-energy-file>`
.. |output_cross| replace:: :ref:`cross-file <user-section-output-keyword-cross-file>`
.. |output_order| replace:: :ref:`order-file <user-section-output-keyword-order-file>`
.. |output_path| replace:: :ref:`pathensemble-file <user-section-output-keyword-pathensemble-file>`
.. |output_restart| replace:: :ref:`restart-file <user-section-output-keyword-restart-file>`
.. |output_traj| replace:: :ref:`trajectory-file <user-section-output-keyword-trajectory-file>`
.. |output_screen| replace:: :ref:`screen <user-section-output-keyword-screen>`


.. _table-keywords-file-output:

.. table:: Keywords for output to files.
   :class: table-striped table-hover

   +------------------+-------------------------------------------------------+
   | Keyword          | Description                                           |
   +==================+=======================================================+
   | |output_backup|  | Determines how we deal with already existing files.   |
   +------------------+-------------------------------------------------------+
   | |output_prefix|  | Add a prefix to output files.                         |
   +------------------+-------------------------------------------------------+
   | |output_energy|  | Set frequency of output to the energy file.           |
   +------------------+-------------------------------------------------------+
   | |output_cross|   | Set frequency of output of crossing data.             |
   +------------------+-------------------------------------------------------+ 
   | |output_order|   | Set frequency of output of order parameter data.      |
   +------------------+-------------------------------------------------------+
   | |output_path|    | Set frequency of output of path ensemble data.        |
   +------------------+-------------------------------------------------------+
   | |output_restart| | Set frequency of restart files overwriting.           |
   +------------------+-------------------------------------------------------+
   | |output_traj|    | Set frequency of output of trajectory data.           |
   +------------------+-------------------------------------------------------+


.. _table-keywords-screen-output:

.. table:: Keywords for output to the screen.
   :class: table-striped table-hover

   +------------------+-------------------------------------------------------+
   | Keyword          | Description                                           |
   +==================+=======================================================+
   | |output_screen|  | Set frequency of output to the screen.                |
   +------------------+-------------------------------------------------------+


Keywords controlling output to files
------------------------------------

For the keywords controlling the frequency of the output
to files, e.g. :ref:`energy-file <user-section-output-keyword-energy-file>`,
writing of such output can be disabled by setting a frequency of ``0`` or less.

Note that |pyretis| will create a set of files for each path ensemble in a RETIS (or TIS)
simulation.

.. _user-section-output-backup:

Keyword backup
^^^^^^^^^^^^^^

.. pyretis-keyword:: backup string

   Determines how we deal with already existing files.
   This can be set to one of the following:

   - ``backup = append``: Output will be appended to already existing
     files.

   - ``backup = overwrite``: Output will **overwrite** already existing files.

   - ``backup = backup``: A backup of existing files will be created.
     Existing files will be renamed with a suffix ``_XXX`` where
     ``XXX`` is a number.


   Default
        The default is ``backup = overwrite``.

.. _user-section-output-prefix:

Keyword prefix
^^^^^^^^^^^^^^

.. pyretis-keyword:: prefix string

   The prefix can be used to add names to the output files, that
   is ``prefix = TEXT`` will prepend ``TEXT`` to output files
   created by |pyretis|. This can, for instance, be used to organize
   your output.

   Default
        The default is ``None``, that is: no prefix is added by
        default.


.. _user-section-output-keyword-energy-file:

Keyword energy-file
^^^^^^^^^^^^^^^^^^^

.. pyretis-keyword:: energy-file integer

   Determines how often we should write to the energy file.
   ``energy-file = N`` will make |pyretis| write to the energy
   file every Nth step.

   Default
        The default is ``energy-file = 1``.


.. _user-section-output-keyword-cross-file:

Keyword cross-file
^^^^^^^^^^^^^^^^^^

.. pyretis-keyword:: cross-file integer

   Determines how often |pyretis| will try to write to the crossing
   file. The crossing file will only be written to if an actual
   crossing has been detected. Note that this only applies to
   simulations for determining the initial flux.

   ``cross-file = N`` will make |pyretis| write
   to the crossing file every Nth step. Normally, this keyword should
   be set to ``1``. This does not mean that we will write to the
   crossing file every step, but only that we will check if we
   can write new crossing data every step. The file will only be
   written to if a crossing has been made.

   Default
        The default is ``cross-file = 1``.


.. _user-section-output-keyword-order-file:

Keyword order-file
^^^^^^^^^^^^^^^^^^

.. pyretis-keyword:: order-file integer

   Determines how often we will write out order parameters to
   the order file.

   ``order-file = N`` will make |pyretis| write to the order
   file every Nth step.

   Default
        The default is ``order-file = 1``.


.. _user-section-output-keyword-restart-file:

Keyword restart-file
^^^^^^^^^^^^^^^^^^^^

.. pyretis-keyword:: restart-file integer

   Determines how often we will write out the restart file.

   ``restart-file = N`` will make |pyretis| overwrite the restart files
   every Nth step. The restart files are: ``pyretis.restart`` in the main
   folder and ``ensemble.restart`` in the ensemble's pathfolders for 
   path ensemble simulations.  

   Default
        The default is ``restart-file = 1``.


.. _user-section-output-keyword-pathensemble-file:

Keyword pathensemble-file
^^^^^^^^^^^^^^^^^^^^^^^^^

.. pyretis-keyword:: pathensemble-file integer

   Determines how often we will write out path data for path
   ensembles. This file contains the basic information about the paths
   (length, maximum/minimum order parameter, what kind of move
   we perform, was the path accepted etc.) and is the **most important**
   file for the analysis. The move type codes and path status (acceptance
   or rejection reason) codes recorded in this file are explained in
   :ref:`user-guide-path-types-rejections`.
   Normally, the frequency of output should be
   set to ``1``, but setting ``pathensemble-file = N`` will
   make |pyretis| write to the order file every Nth step.

   .. WARNING::
       Setting a value different from ``1`` is **not** accounted
       for in the analysis and is **not** recommended.

   Default
       The default is ``pathensemble-file = 1``.


.. _user-section-output-keyword-trajectory-file:

Keyword trajectory-file
^^^^^^^^^^^^^^^^^^^^^^^

.. pyretis-keyword:: trajectory-file integer

   Determines how often we should write trajectory data to a
   trajectory output file.

   .. WARNING::
       Setting this to a *small* value (say 1) is **not** recommended
       when using external engines as this may quickly use up all the
       space on your hard drive. Setting ``trajectory-file = N`` means
       that |pyretis| will store a trajectory every Nth step for each
       path ensemble in a RETIS/TIS simulation. You are advised to
       investigate typical sizes for the trajectories before setting
       this frequency output.

   Default
      The default is ``trajectory-file = 100``.


Keywords controlling output to the screen
-----------------------------------------

Note that if you run the |pyretis| application and choose to
use the progress bar, that is:

.. code-block:: pyretis

   pyretisrun -i input.rst -p

then writing of results to the screen will be **disabled**.


.. _user-section-output-keyword-screen:

Keyword screen
^^^^^^^^^^^^^^

.. pyretis-keyword:: screen integer

   Determines how often output from a simulation should be written to
   the screen during a simulation.

   Default
        The default is ``screen = 10``.
