# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test setup functions."""
import os
import logging
import pytest
from pyretis.setup.common import (check_settings,
                                  create_external)


logging.disable(logging.CRITICAL)
LOCAL_DIR = os.path.abspath(os.path.dirname(__file__))


class TestCheckSettings:
    """Test that check_settings work as intented."""

    def test_check_settings(self):
        """Test that we can import."""
        settings = {'one': 1, 'two': 2, 'three': 3}
        required = ['three', 'two']
        result, not_found = check_settings(settings, required)
        assert result
        assert 0 == len(not_found)
        extra = ['three fiddy']
        required += extra
        result, not_found = check_settings(settings, required)
        assert not result
        for i, j in zip(not_found, extra):
            assert i == j


class TestCreateExternal:
    """Test that we can import and create from other modules."""

    def test_create_from_module(self):
        """Test that we can import."""
        module = os.path.join(LOCAL_DIR, 'fooengine.py')

        settings = {}
        obj = create_external(settings, 'foo', None, None, key_settings=None)
        assert obj is None

        settings = {'foo': {'module': module, 'class': 'FooEngine'}}
        with pytest.raises(ValueError):  # missing an argument:
            obj = create_external(settings, 'foo', None, [],
                                  key_settings=None)
        settings['foo']['timestep'] = 1
        create_external(settings, 'foo', None, [], key_settings=None)

        with pytest.raises(ValueError):  # not callable:
            create_external(settings, 'foo', None, ['foo_bar'],
                            key_settings=None)

        settings = {'foo': {'module': 'three fiddy', 'class': 'NoClassHere'},
                    'simulation': {}}
        with pytest.raises(ValueError):
            create_external(settings, 'foo', None, [], key_settings=None)
