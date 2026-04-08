# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Tests for the properties module."""
import os
import pytest
import numpy as np
from pyretis.core.properties import Property


def test_property_initialization():
    """Test initialization of the Property class."""
    prop = Property(desc="Energy of the system")
    assert prop.desc == "Energy of the system"
    assert prop.val == []
    assert prop.mean == 0.0
    assert prop.variance == 0.0
    assert prop.nval == 0


def test_property_add():
    """Test adding values to a property."""
    prop = Property()
    values = [1.0, 2.0, 3.0, 4.0, 5.0]
    for val in values:
        prop.add(val)

    assert prop.val == values
    assert prop.nval == len(values)
    assert prop.mean == np.mean(values)
    # ddof=1 for sample variance in numpy matches n-1 in the code
    assert prop.variance == pytest.approx(np.var(values, ddof=1))


def test_property_variance_single_value():
    """Test variance with a single value (should be inf)."""
    prop = Property()
    prop.add(10.0)
    assert prop.variance == float('inf')


def test_property_str():
    """Test string representation of a property."""
    prop = Property(desc="system energy")
    assert str(prop) == "system energy"


def test_property_dump_to_file(tmp_path):
    """Test dumping property values to a file."""
    prop = Property()
    values = [1.0, 2.5, 3.14]
    for val in values:
        prop.add(val)

    dump_file = tmp_path / "property.txt"
    prop.dump_to_file(str(dump_file))

    read_values = np.loadtxt(str(dump_file))
    assert np.allclose(read_values, values)
