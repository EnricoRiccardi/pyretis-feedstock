# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the SimulationTask class."""
import logging
import pytest
from pyretis.simulation.simulation_task import SimulationTask
from .help import turn_on_logging
logging.disable(logging.CRITICAL)


class TestSimulationTask:
    """Run the tests for SimulationTask."""

    def test_create_task(self, caplog):
        """Test that we can create simulation tasks."""

        def some_function():  # pylint: disable=missing-docstring
            return 'Hello!'

        SimulationTask(some_function)
        # Test if we give wrong arguments:
        with pytest.raises(AssertionError):
            SimulationTask(some_function, args=['dummy'])
        # Test if we give a non-callable:
        not_a_function = 100
        with pytest.raises(AssertionError):
            SimulationTask(not_a_function, args=['dummy'])

        def some_function2(test='hello'):  # pylint: disable=missing-docstring
            return test

        inkw = {'test': 'yes'}
        task2 = SimulationTask(some_function2, kwargs=inkw)
        for key, val in inkw.items():
            assert key in task2.kwargs
            assert val == task2.kwargs[key]
        # Test giving too many kwargs
        inkw = {'test': 'yes', 'more': 'indeed'}
        with pytest.raises(AssertionError):
            with turn_on_logging():
                with caplog.at_level(logging.WARNING):
                    SimulationTask(some_function2, kwargs=inkw)
        # Test giving kwargs when None are expected
        inkw = {'missing': 'indeed'}
        with pytest.raises(AssertionError):
            with turn_on_logging():
                with caplog.at_level(logging.WARNING):
                    SimulationTask(some_function, kwargs=inkw)

    def test_execute(self, caplog):
        """Test that we can execute the task."""

        def some_function():  # pylint: disable=missing-docstring
            return [10, 9, 8]

        def some_function2(pos):  # pylint: disable=missing-docstring
            return pos * pos

        def some_function3(pos, exp=10):  # pylint: disable=missing-docstring
            return pos**exp

        def some_function4(exp=10):  # pylint: disable=missing-docstring
            return 2**exp

        task1 = SimulationTask(some_function, when={'every': 2},
                               result='the-stuff', first=True)
        step = {'step': 10, 'start': 3, 'stepno': 7}
        result = task1.execute(step)
        assert result is None
        step = {'step': 10, 'start': 3, 'stepno': 8}
        result = task1.execute(step)
        for i, j in zip(result, [10, 9, 8]):
            assert i == j

        var = 10
        task2 = SimulationTask(some_function2, when={'every': 2},
                               args=[var], result='the-stuff', first=True)
        result = task2.execute(step)
        assert result == 100

        var = 5
        task3 = SimulationTask(some_function3, when={'every': 2},
                               args=[var], kwargs={'exp': 2},
                               result='the-stuff', first=False)
        result = task3.execute(step)
        assert not task3.run_first()
        assert result == 25

        task4 = SimulationTask(some_function4,
                               kwargs={'exp': 2},
                               result='the-stuff', first=True)
        result = task4(step)
        assert result == 4
        assert task4.run_first()
        assert task4.result == 'the-stuff'

    def test_when_change(self, caplog):
        """Test that we can change the "when" property."""
        def function(var):  # pylint: disable=missing-docstring
            return var * 10
        task = SimulationTask(function, args=[2], result='times-10')
        step = {'step': 10, 'start': 3, 'stepno': 7}
        result = task(step)
        assert result == 20
        assert task.when is None
        task.when = None
        assert task.when is None
        task.when = {'every': 1}
        assert task.when == {'every': 1}
        task.when = None
        assert task.when is None
        task.when = {'every': 1, 'all-the-time': 100}
        assert task.when == {'every': 1}

    def test_task_dict(self, caplog):
        """Test that the dictionary returned from a task is ok."""
        def function(var):  # pylint: disable=missing-docstring
            return var * 2
        task = SimulationTask(function, args=[21], result='double')
        task.when = {'every': 2}
        task_dict = task.task_dict()
        assert task_dict['kwargs'] is None
        assert task_dict['when'] == {'every': 2}
        assert not task_dict['first']
        assert task_dict['result'] == 'double'
        assert task_dict['args'] == [21]
        assert task_dict['func'] == function
