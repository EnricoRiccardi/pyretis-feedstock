# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test generic initiation methods."""
import logging
from io import StringIO
import pytest
from unittest.mock import patch
from pyretis.initiation import (
    initiate_kick,
    initiate_load,
    initiate_restart,
    get_initiation_method,
)


logging.disable(logging.CRITICAL)


class TestInitiate:
    """Run the tests for the reading external trajectories."""

    def test_get_initiation_method(self):
        """Test that we can get the method from settings."""
        cases = [
            {
                'settings': {'initial-path': {'method': 'kick'}},
                'method': initiate_kick,
            },
            {
                'settings': {'initial-path': {'method': 'KICK'}},
                'method': initiate_kick,
            },
            {
                'settings': {'initial-path': {'method': 'load'}},
                'method': initiate_load,
            },
            {
                'settings': {'initial-path': {'method': 'loAD'}},
                'method': initiate_load,
            },
            {
                'settings': {'initial-path': {'method': 'restart'}},
                'method': initiate_restart,
            },
            {
                'settings': {'initial-path': {'method': 'reSTART'}},
                'method': initiate_restart,
            },

        ]
        with patch('sys.stdout', new=StringIO()):
            for case in cases:
                method = get_initiation_method(case['settings'])
                assert method == case['method']
        settings = {'initial-path': {'method': 'santas little failer'}}
        with pytest.raises(ValueError):
            get_initiation_method(settings)
