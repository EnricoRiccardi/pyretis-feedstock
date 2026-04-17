.. _developer-guide-index:

###############
Developer Guide
###############

|pyretis| is open source and hosted on `gitlab <https://gitlab.com>`_.
We are happy to include more developers and we want active users.
If you wish to contribute to the |pyretis| project, please read through
the code :ref:`guidelines <developer-guide-guidelines>` and
the short description on :ref:`contributing <developer-guide-contributing>`
given below.

.. _developer-guide-guidelines:

Guidelines
==========

The guidelines can be summarised as follows:

- Check that your code follows
  `pep8 <https://www.python.org/dev/peps/pep-0008/>`_.

- Document your code using `NumPy style docstrings
  <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_

- Use a ``logger`` rather than ``print()`` in libraries.


|pyretis| follows the Python
`pep8 <https://www.python.org/dev/peps/pep-0008/>`_ style guide
(see also `pep8.org <http://pep8.org/>`_)
and new code should be checked with the
`pep8 style guide checker <https://pycodestyle.readthedocs.io/en/latest/>`_
and `pylint <http://www.pylint.org/>`_:

.. code-block:: pyretis

    pycodestyle source_file.py
    pylint source_file.py

or other tools like `PyChecker <http://pychecker.sourceforge.net/>`_ or
`pyflakes <https://pypi.python.org/pypi/pyflakes>`_.
`NumPy's <http://www.numpy.org/>`_ imports can be a bit
tricky understand so you can help pylint by doing

.. code-block:: pyretis

    pylint --extension-pkg-whitelist=numpy source_file.py

There is also the

.. code-block:: python

    import this

which can be useful to remember.

The |pyretis| project is documented using the
`NumPy/SciPy documentation standard <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_
and contributors are requested to familiarise themselves with this style and *use it*.
For user-facing documentation, also see the
:ref:`documentation terminology guide <developer-guide-docs-style>`.
Documentation style can be checked with
`pydocstyle <https://github.com/PyCQA/pydocstyle>`_

.. code-block:: pyretis

    pydocstyle source_file.py


We also try to avoid ``print()`` statements in the libraries and reserve such
statements for console output from command line scripts/programs. Output, like
debug information, warnings etc. can be handled by using
a `logger <https://docs.python.org/3/howto/logging.html>`_.
In a library, we typically just import and set up the logger
by doing:

.. code-block:: python

    import logging
    logger = logging.getLogger(__name__)  # pylint: disable=C0103
    logger.addHandler(logging.NullHandler())

We can then use it as follows:

.. code-block:: python

    # Report some information
    logger.info('Up and running')
    logger.debug('Debug info')
    logger.warning('This is a warning')
    logger.critical('Something critical happened!')
    logger.error('An error occurred!')

Please note that we do not add any particular handlers here - it is up to the **user** to
define what should happen when a logging event occurs.


.. _developer-guide-contributing:

Contributing to |pyretis|
=========================

There is a nice guide on github about `contributing to open source projects
<https://guides.github.com/activities/contributing-to-open-source/>`_.
In |pyretis| we largely follow this approach and the issue tracker is used
for reporting bugs and issues and requesting new features. This
section on contributing is based on
the `description for gitlab <https://gitlab.com/gitlab-org/gitlab-ce/blob/master/CONTRIBUTING.md>`_.

Before contributing,
please read the short guidelines given below
for :ref:`reporting bugs and issues <developer-guide-bugs>`
and for :ref:`requesting new features <developer-guide-new-features>`.

.. _developer-guide-bugs:

Bugs and issues
---------------

If you find a bug in |pyretis|, please `create an issue
<https://gitlab.com/pyretis/pyretis/issues>`_ using the following
template:

.. code-block:: restructuredtext

    Summary
    -------

    One sentence summary of the issue: What goes wrong, what did you expect to happen

    Steps to reproduce
    ------------------

    Describe how the issue can be reproduced.

    Expected behavior
    -----------------

    Describe what behavior you expected instead of what actually happens.

    Relevant logs
    -------------

    Add logs from the output.

    Output of checks
    ----------------

    Be sure that all the tests pass before filing the issue.

    Possible fixes
    --------------

    If you can, link to the line of code that might be responsible for the problem.

If you wish to fix the bug yourself, please follow the approach described for
:ref:`merge requests <developer-guide-merge-requests>` below.

.. _developer-guide-new-features:

New features
------------

If you are missing some functionality in |pyretis| you can create
a new issue in the issue tracker and
label it as a `feature proposal <https://gitlab.com/pyretis/pyretis/issues?label_name=feature+proposal>`_.
If you do not have rights to add the ``feature proposal`` label, you can ask one
of the core members of the |pyretis| project to add it for you.

You can also implement the changes you want yourself as
:ref:`described below <developer-guide-merge-requests>`.
We cannot promise that we will automatically
include your work in |pyretis| but we are happy to have active users
and we will consider your contribution.
So, when you are finished with your work please make a merge request.


.. _developer-guide-merge-requests:

Merge requests
--------------

The general approach for making your bug-fix or new feature available
in the |pyretis| project is as follows:

1. Fork [#]_ |pyretis| into your own personal space.

2. Create a new branch and work on your feature.

3. Create tests for your feature. Especially, if you are fixing a bug,
   create a test that shows where the bug is causing |pyretis| to fail and
   how your code is fixing this.

4. Squash all your commits into one and push to your fork.

5. Submit a merge request to the master branch using the template given below.

6. Wait for review and the following discussion :-).

.. code-block:: restructuredtext

    Summary
    -------

    What does this merge request do?

    Justification
    -------------

    Why is this merge request needed?

    Description of the new code
    ---------------------------

    A short description of the new code.
    Are there points in the code the reviewer needs to double check?

    References
    ----------

    Add references to relevant issue numbers.

After submitting the merge request the code will be reviewed [#]_ by (another)
member of the |pyretis| team.


.. [#] Gitlab documentation on forking: http://doc.gitlab.com/ce/workflow/forking_workflow.html
.. [#] https://github.com/thoughtbot/guides/tree/master/code-review


.. toctree::
    :maxdepth: 2

    docs-style.rst
