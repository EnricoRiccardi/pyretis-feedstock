# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Dummy order parameter for running tests."""
import numpy as np
from pyretis.orderparameter import OrderParameter


class FooOrderParameter(OrderParameter):  # pylint: disable=abstract-method
    """FooOrderParameter(OrderParameter) - Dummy order parameter for tests."""

    def __init__(self, name, desc='Dummy order parameter'):
        super().__init__(description=desc)
        self.name = name

    def calculate(self, system):
        """Return a valid single-frame order-parameter list."""
        return [0.0]


class BarOrderParameter:  # pylint: disable=too-few-public-methods
    """BarOrderParameter - Dummy order parameter for tests."""

    def __init__(self, description='Dummy test order parameter'):
        self.description = description


class BazOrderParameter:  # pylint: disable=too-few-public-methods
    """BazOrderParameter - Dummy order parameter for tests."""

    def __init__(self, description='Dummy test order parameter'):
        self.description = description

    def calculate(self, system):
        """Obtain the order parameter."""


class ScalarOrderParameter(OrderParameter):  # pylint: disable=abstract-method
    """External order parameter returning a scalar."""

    def __init__(self):
        super().__init__(description='Scalar order parameter')

    def calculate(self, system):
        """Return a scalar value for one frame."""
        return 1.25


class ArrayOrderParameter(OrderParameter):  # pylint: disable=abstract-method
    """External order parameter returning a NumPy array."""

    def __init__(self):
        super().__init__(description='Array order parameter')

    def calculate(self, system):
        """Return a flat array for one frame."""
        return np.array([1.25, 2.5])


class MultiFrameOrderParameter(  # pylint: disable=abstract-method
    OrderParameter
):
    """External order parameter incorrectly returning multiple frames."""

    def __init__(self):
        super().__init__(description='Bad multi-frame order parameter')

    def calculate(self, system):
        """Return multiple frames instead of one frame."""
        return [[1.25], [2.5]]
