"""Dummy potential for testing."""
from pyretis.forcefield.potential import PotentialFunction


__all__ = ['FooPotential']


class FooPotential(PotentialFunction):
    """FooPotential(PotentialFunction) - Dummy potential for tests."""

    def __init__(self, a=0.0, desc='Dummy potential'):
        super().__init__(dim=1, desc=desc)
        self.params = {'a': a}

    def potential(self, pos):
        """Evaluate the potential."""
        raise NotImplementedError

    def force(self, pos):
        """Evaluate the force."""
        raise NotImplementedError

    def potential_and_force(self, pos):
        """Evaluate the potential and force."""
        raise NotImplementedError


class BarPotential(PotentialFunction):
    """BarPotential(PotentialFunction) - Dummy potential for tests."""

    def __init__(self, a=0.0, desc='Dummy potential'):
        super().__init__(dim=1, desc=desc)
        self.params = {'a': a}

    def potential(self, pos):
        """Evaluate the potential."""
        raise NotImplementedError

    def force(self, pos):
        """Evaluate the force."""
        raise NotImplementedError
