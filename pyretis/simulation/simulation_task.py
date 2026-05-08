# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Definition of a class for simulation tasks.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SimulationTask (:py:class:`.SimulationTask`)
    A class representing a simulation task.

SimulationTaskList (:py:class:`.SimulationTaskList`)
    A class for representing a list of simulation tasks. This class
    defines functionality for adding tasks from a dictionary description.

"""
import logging
from pyretis.core.common import inspect_function
from pyretis.inout.simulationio import Task
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


def _check_args(function, given_args=None, given_kwargs=None):
    """Check consistency for function and the given (keyword) arguments.

    Here we assume that the arguments are given in a list and that
    the keyword arguments are given as a dictionary. The function
    `inspect.getargspec` is used to check the input function.

    Parameters
    ----------
    function : callable
        The function we will inspect.
    given_args : list, optional
        A list of the arguments to pass to the function. 'self' will not
        be considered here since it passed implicitly.
    given_kwargs : dict, optional
        A dictionary with keyword arguments.

    Returns
    -------
    out : boolean
        False if there is some inconsistencies, i.e. when the calling
        of the given `function` will probably fail. True otherwise.

    """
    arguments = inspect_function(function)
    args = [arg for arg in arguments['args'] if arg != 'self']
    defaults = list(arguments['kwargs'])
    # first test, do we give correct number of required arguments?
    if given_args is not None:
        given = len(given_args)
    else:
        given = 0
    if len(args) != given:
        msgtxt = 'Wrong number of arguments given'
        logger.warning(msgtxt)
        return False
    # Check kwargs but only check in case some kwargs are given here.
    # If they are not given, we assume that the user knows what's happening
    # and that the default kwargs will be used.
    if given_kwargs is not None:
        if defaults:
            extra = [key for key in given_kwargs if key not in defaults]
            if extra:
                msg = [f'Task Keyword arguments: {defaults}']
                msg += [f'Unexpected keyword argument: {extra}']
                msgtxt = '\n'.join(msg)
                logger.warning(msgtxt)
                return False
        else:
            msgtxt = 'Unexpected keyword argument!'
            logger.warning(msgtxt)
            return False
    return True


class SimulationTask(Task):
    """Representation of simulation tasks.

    This class defines a task object. A task is executed at specific
    points, at regular intervals etc. in a simulation. A task will
    typically provide a result, but it does not need to. It can simply
    just alter the state of the passed argument(s).

    Attributes
    ----------
    function : function
        The function to execute.
    when : dict
        Determines when the task should be executed.
    args : list
        List of arguments to the function.
    kwargs : dict
        The keyword arguments to the function.
    first : boolean
        True if this task should be executed before the first
        step of the simulation.
    result : string
        This is a label for the result created by the task.

    """

    def __init__(self, function, args=None, kwargs=None, when=None,
                 result=None, first=False):
        """Initialise the task.

        Parameters
        ----------
        function : callable
            The function to execute.
        args : list, optional
            List of arguments to the function.
        kwargs : dict, optional
            The keyword arguments to the function.
        when : dict, optional
            Determines if the task should be executed.
        result : string, optional
            This is a label for the result created by the task.
        first : boolean, optional
            True if this task should be executed before the first
            step of the simulation.

        """
        if not callable(function):
            msg = 'The given function for the task is not callable!'
            raise AssertionError(msg)
        ok_to_add = _check_args(function, given_args=args, given_kwargs=kwargs)
        if not ok_to_add:
            msg = 'Wrong arguments or keyword arguments!'
            raise AssertionError(msg)
        super().__init__(when)
        self.function = function
        self.args = args
        self.kwargs = kwargs
        self._result = result
        self.first = first

    def execute(self, step):
        """Execute the task.

        Parameters
        ----------
        step : dict of ints
            The keys are:

            * 'step': the current cycle number.
            * 'startcycle': the cycle number at the start.
            * 'stepno': the number of cycles we have performed so far.

        Returns
        -------
        out : unknown type
            The result of running `self.function`.

        """
        args = self.args
        kwargs = self.kwargs
        if self.execute_now(step):
            if args is None:
                if kwargs is None:
                    return self.function()
                return self.function(**kwargs)
            if kwargs is None:
                return self.function(*args)
            return self.function(*args, **kwargs)
        return None

    @property
    def result(self):
        """Return the result label."""
        return self._result

    def run_first(self):
        """Return True if task should be executed before first step."""
        return self.first

    def task_dict(self):
        """Return a dict representing the task."""
        return {'func': self.function, 'args': self.args,
                'kwargs': self.kwargs, 'when': self.when,
                'result': self.result, 'first': self.first,
                'func-name': self.function.__name__}

    def __call__(self, step):
        """Execute the task.

        Parameters
        ----------
        step : dict of ints
            The keys are:

            * 'step': the current cycle number.
            * 'startcycle': the cycle number at the start.
            * 'stepno': the number of cycles we have performed so far.

        Returns
        -------
        out : unknown type
            The result of `self.execute(step)`.

        """
        return self.execute(step)

    def __str__(self):
        """Output info about the task."""
        msg = ['Task:']
        msg += [f' -> Function name: {self.function.__name__}']
        msg += [f' -> Function args: {self.args}']
        msg += [f' -> Function kwargs: {self.kwargs}']
        msg += [f' -> Execute when: {self.when}']
        msg += [f' -> Execute at first: {self.first}']
        msg += [f' -> Result: {self._result}']
        return '\n'.join(msg)
