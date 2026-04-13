.. _user-guide-application:

The |pyretis| application
=========================

The |pyretis| application, ``pyretisrun``, is used to
execute simulations defined in input files. The general syntax for
executing is:

.. code-block:: pyretis

   pyretisrun [-h] -i INPUT [-V] [-f LOG_FILE] [-l LOG_LEVEL] [-p]

where the different arguments are described in :numref:`tableappargument`
and :numref:`tableloglevel`.

Example use
-----------

Here are some examples of the use of the application:

* To run the simulation defined in a file named ``retis.rst``, make a log file
  for output named ``retis.log`` and display a progress bar:

  .. code-block:: pyretis

     pyretisrun -i retis.rst -f retis.log -p

* To run the simulation defined in a file named ``test.rst`` and display debug messages:

  .. code-block:: pyretis

     pyretisrun -i test.rst -l DEBUG

* To run the simulation defined in a file named ``input.rst``, display only critical
  messages, write to a log file named ``mylog.log`` and display a progress bar:

  .. code-block:: pyretis

     pyretisrun -i input.rst -l CRITICAL -p -f mylog.log

The long version of the arguments may also be used, for instance, the previous example
can also be written more verbose as

.. code-block:: pyretis

   pyretisrun --input input.rst --log_level CRITICAL --progress --log_file mylog.log



Input arguments
---------------

.. _tableappargument:

.. table:: Description of input arguments for |pyretis|.
   :class: table-striped table-hover

   +-------------------------------------+--------------------------------------------------+
   | Argument                            | Description                                      |
   +=====================================+==================================================+
   | -h, --help                          | Show the help message and exit                   |
   +-------------------------------------+--------------------------------------------------+
   | -i INPUT, --input INPUT             | Location of the input file.                      |
   +-------------------------------------+--------------------------------------------------+
   | -V, --version                       | Show the version number and exit.                |
   +-------------------------------------+--------------------------------------------------+
   | -f LOG_FILE, --log_file LOG_FILE    | Specify location of log file to write.           |
   +-------------------------------------+--------------------------------------------------+
   | -l LOG_LEVEL, --log_level LOG_LEVEL | Specify log level for the log file               |
   |                                     | (see :numref:`tableloglevel`).                   |
   +-------------------------------------+--------------------------------------------------+
   | -p, --progress                      | Display a progress meter instead of text output. |
   +-------------------------------------+--------------------------------------------------+

.. _tableloglevel:

.. table:: Possible settings for the ``LOG_LEVEL``.
   :class: table-striped table-hover

   +-----------+---------------------------------------------------------------+
   | LOG_LEVEL | Description                                                   |
   +===========+===============================================================+
   | CRITICAL  | Only display critical messages.                               |
   +-----------+---------------------------------------------------------------+
   | ERROR     | Display critical and error messages.                          |
   +-----------+---------------------------------------------------------------+
   | WARNING   | Display critical and error messages and warnings.             |
   +-----------+---------------------------------------------------------------+
   | INFO      | Display all of the above and some info messages.              |
   +-----------+---------------------------------------------------------------+
   | DEBUG     | Display all of the above and additional debug messages.       |
   +-----------+---------------------------------------------------------------+

PyRETIS also defines three custom levels used internally for console output
(they are not intended as ``-l`` arguments but are listed here for completeness):

.. table:: Custom |pyretis| log levels.
   :class: table-striped table-hover

   +-----------+-------+---------------------------------------------------------------+
   | LOG_LEVEL | Value | Description                                                   |
   +===========+=======+===============================================================+
   | PROGRESS  | 25    | Progress messages shown on the console (default level).       |
   +-----------+-------+---------------------------------------------------------------+
   | BANNER    | 26    | Logo and header messages displayed in cyan at startup.        |
   +-----------+-------+---------------------------------------------------------------+
   | REFERENCE | 27    | Citation and URL messages displayed at startup/shutdown.      |
   +-----------+-------+---------------------------------------------------------------+

