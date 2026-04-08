# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A test for the potential base class."""
import logging
import pytest
from pyretis.forcefield.potential import PotentialFunction
logging.disable(logging.CRITICAL)


class TestPotentialFunction:
    """Test the PotentialFunction class."""

    def test_potential_function(self):
        """Test PotentialFunction class."""
        pot = PotentialFunction()
        params = {'a': 10}
        pot.set_parameters(params)
        assert not pot.params
