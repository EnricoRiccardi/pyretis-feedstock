# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Dummy engines for tests."""
from pyretis.engines import MDEngine


__all__ = []


class FooEngine(MDEngine):
    """FooIntegrator(MDEngine) - Dummy engine for tests."""

    foo_bar = 1.0

    def __init__(self, timestep, extra=0.0,
                 description='Dummy test engine'):
        super().__init__(timestep, description, dynamics=None)
        self.extra = extra

    def integration_step(self, system):
        """Perform the integration step."""
        raise NotImplementedError


class BarEngine:  # pylint: disable=too-few-public-methods
    """BarEngine - Dummy engine for tests."""

    def __init__(self, description='Dummy test engine'):
        self.description = description


class BazEngine:  # pylint: disable=too-few-public-methods
    """BazEngine - Dummy engine for tests."""

    def __init__(self, description='Dummy test enginer'):
        self.description = description
        self.integration_step = 'fake_step'
