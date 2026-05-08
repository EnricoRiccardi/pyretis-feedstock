# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test pyretis.inout.screen"""
import unittest
import pytest
from unittest.mock import patch
from io import StringIO
from pyretis.inout.screen import ScreenOutput


class TestScreenOutput:
    """Test the screen outputter."""

    def test_screen_output(self):
        """Test that the screen outputter works as intended."""
        screen = ScreenOutput(formatter=None)
        test = 'This is a test'
        with patch('sys.stdout', new=StringIO()) as txt:
            write = screen.write(test)
            assert write
            assert test == txt.getvalue().strip()
        with patch('sys.stdout', new=StringIO()):
            write = screen.write(None)
            assert not write
        with patch('sys.stdout', new=StringIO()) as txt:
            write = screen.write(test, end='\n')
            assert write
            assert test == txt.getvalue().split('\n')[0]
