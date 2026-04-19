.. _user-guide-install:

Obtaining and installing |pyretis|
==================================

|pyretis| and its :ref:`requirements <user-guide-install-requirements>` can be
installed :ref:`using pip <user-guide-install-pip>`
or :ref:`using conda <user-guide-install-conda>`
as described below.
We also describe how it can be installed directly from the source code
:ref:`using git <user-guide-install-develop>`.

|

.. _user-guide-prereq:

Prerequisites
-------------

|pyretis| requires **Python 3.11** or higher. Please make sure that you have
an updated version of Python installed on your system. |pyretis| depends on several
other Python packages, however, these should be installed automatically as part
of the installation of |pyretis|. A  list of the dependencies :ref:`can be found below <user-guide-install-requirements>`.

.. _user-guide-install-pip:

|

Installing via ``pip``
----------------------
The latest version of |pyretis| can be installed using `pip <https://pypi.org/project/pyretis/>`_
with the following command:

.. code-block:: pyretis

    python -m pip install git+https://gitlab.com/pyretis/pyretis.git

A direct installation of PyRETIS 3:

.. code-block:: pyretis

    python -m pip install pyretis
 
|pyretis| offers an analysis and visualization tool, named |pyvisa|. Its GUI requires
PyQt5, `mdtraj <http://mdtraj.org>`_, and `graphviz <https://pypi.org/project/graphviz/>`_.
To install them via pip:

.. code-block:: pyretis

    python -m pip install pyqt5 mdtraj graphviz

**Note:** Please make sure that you are using a Python environment
of version 3.11 or newer. On some systems, **Python 2** may still be the default
version and ``pip`` might actually try to install for **Python 2**.
In this case, use ``pip3`` in the command above.

|

.. _user-guide-install-conda:


Installing via ``conda``
------------------------
|pyretis| can be installed using `conda <https://anaconda.org/conda-forge/pyretis>`_ 
with the following commands:

.. code-block:: pyretis

    conda create --name pyretis python=3.12
    conda activate pyretis
    conda install pyretis -c conda-forge

|pyretis| offers an analysis and visualization tool, named |pyvisa|. Its GUI requires
PyQt5, mdtraj, and graphviz to be executed. To install them via conda:

.. code-block:: pyretis

    conda install pyqt mdtraj python-graphviz -c conda-forge


**Note:** Since |pyretis| will **only work with Python 3.11 or newer**,
please make sure that you are using an environment with a recent version
of Python.

|

Optional: Installing from the |pyretis| git repository
------------------------------------------------------

Previous versions and the latest (*possibly unstable*) sources can be
obtained using ``git``:

.. code-block:: pyretis

   git clone git@gitlab.com:pyretis/pyretis.git

or,

.. code-block:: pyretis

   git clone https://gitlab.com/pyretis/pyretis.git

After cloning the repository, |pyretis| can be
installed via ``pip`` (after navigating to the main source directory):

.. code-block:: pyretis

   pip install .

**Note:** If you want to be able to make modifications locally and to apply them, you can 
install the package with the ``-e`` option. |pyretis| will be executed with the files present
in the folder you created.

To be able to run |pyvisa|'s GUI, you need to install ``pyqt5``, ``mdtraj``,
and ``graphviz``:

.. code-block:: pyretis

   python -m pip install pyqt5 mdtraj graphviz

.. _user-guide-install-develop:

|

Optional: Installing a development version from the |pyretis| git repository
----------------------------------------------------------------------------

After cloning the repository as described above, check out the development
branch you are interested in, e.g.:

.. code-block:: pyretis

   git checkout develop 

Then install the development requirements (these are defined in the file
``requirements-dev.txt``, note that they include pyqt5, which is not 
supported in some environments):

.. code-block:: pyretis

    pip install -r requirements-dev.txt

Finally, install |pyretis| using:

.. code-block:: pyretis

    pip install -e .

.. _user-guide-install-test:

|

Testing your installation
-------------------------

After installing from the |pyretis| source, your installation can be
tested by running the tests from the main directory:

.. code-block:: pyretis

   python -m unittest discover -v -s test

.. _user-guide-install-requirements:

|

Requirements for |pyretis|
--------------------------

In order to run |pyretis|, several Python libraries are needed, for instance
`SciPy <http://www.scipy.org/>`_, `NumPy <http://www.numpy.org/>`_, and `matplotlib <http://matplotlib.org/>`_
(see also the information on `installing the SciPy Stack <http://www.scipy.org/install.html>`_).
A list of the requirements can be found in the file
:download:`requirements.txt </_static/files/requirements.txt>` in the source
code directory. These packages can be installed by:

.. code-block:: pyretis

    python -m pip install -r requirements.txt

after downloading the :download:`requirements.txt </_static/files/requirements.txt>` file.
This should be automatically done if you
are installing |pyretis| using ``pip``/``conda``.
Notes: (1) the analysis and visualization package |pyvisa| requires PyQt5, mdtraj,
and graphviz, which must be installed separately as described above.
