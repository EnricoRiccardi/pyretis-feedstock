# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Tests for pyretis.pyvisa.info."""
from pyretis.pyvisa.info import PROGRAM_NAME, URL, GIT_URL, CITE, LOGO


class TestInfo:
    """Tests for info constants."""

    def test_program_name(self):
        """Program name should be a non-empty string."""
        assert isinstance(PROGRAM_NAME, str)
        assert len(PROGRAM_NAME) > 0

    def test_urls(self):
        """URL and GIT_URL should be non-empty strings."""
        assert isinstance(URL, str) and len(URL) > 0
        assert isinstance(GIT_URL, str) and len(GIT_URL) > 0

    def test_cite(self):
        """CITE should contain at least one reference."""
        assert isinstance(CITE, str)
        assert 'doi' in CITE.lower() or 'J. Comput. Chem.' in CITE

    def test_logo(self):
        """LOGO should be a non-empty string."""
        assert isinstance(LOGO, str)
        assert len(LOGO) > 0
