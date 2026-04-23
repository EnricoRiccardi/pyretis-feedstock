Extending PyRETIS with FORTRAN
==============================

This folder contains an example of extending PyRETIS with a new
force field which is implemented in FORTRAN.

The FORTRAN code must be compiled before it can be executed and this
is done by running ``make``.

Note that you may have to edit the Makefile depending on your version
of Python. In case you have several Python versions installed, you
may have to specify the precise version of ``f2py`` to use as this should
correspond to the Python version you are currently using. This can for instance
be: ``f2py``, ``f2py3`` or ``f2py3.X`` (where ``X`` denotes the minor version number).

