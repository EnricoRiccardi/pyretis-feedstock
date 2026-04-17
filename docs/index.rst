.. pyretis documentation master file, created by
   sphinx-quickstart on Fri Jun 19 11:01:24 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. title:: Home

#####################
Welcome to |pyretis|!
#####################


.. container:: jumbotron

  .. image:: /_static/img/logo3.png
     :width: 75%

  |pyretis| is a Python library for **rare event molecular simulations**
  with emphasis on methods based on
  :ref:`transition interface sampling <user-guide-tis-theory>`
  and :ref:`replica exchange transition interface sampling <user-guide-retis-theory>`.

  |pyretis| contains |pyvisa|, a Python library for visualizing
  and analysing path sampling results with a GUI interface.

  |pyretis| is :ref:`open source <pyretis-license>`, easy to use
  and can be interfaced with other simulation packages such as GROMACS, LAMMPS
  or CP2K.

  You can use the |pyretis| :ref:`library <api-doc>` to set up tailored
  simulations or you can use a Python flavored :ref:`input file <user-guide-input>`
  to run different kinds of path sampling simulations. Please see the
  :ref:`user guide <user-guide-index>` for information about the usage and
  how to :ref:`obtain and install <user-guide-install>` |pyretis|.

.. toctree::
    :maxdepth: 2
    :hidden:

    about/index.rst
    user/index.rst
    api/pyretis.rst
    examples/index.rst
    developer/index.rst

