# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Tests for the Plotter base class."""
import sys
import types
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest


MODULE_PATH = (
    Path(__file__).resolve().parents[2]
    / 'pyretis'
    / 'inout'
    / 'plotting'
    / 'plotting.py'
)
sys.modules.setdefault('pyretis', types.ModuleType('pyretis'))
sys.modules.setdefault('pyretis.inout', types.ModuleType('pyretis.inout'))
sys.modules.setdefault(
    'pyretis.inout.plotting',
    types.ModuleType('pyretis.inout.plotting')
)
MODULE_SPEC = spec_from_file_location(
    'pyretis.inout.plotting.plotting',
    MODULE_PATH,
)
PLOTTING_MODULE = module_from_spec(MODULE_SPEC)
sys.modules['pyretis.inout.plotting.plotting'] = PLOTTING_MODULE
assert MODULE_SPEC.loader is not None
MODULE_SPEC.loader.exec_module(PLOTTING_MODULE)
Plotter = PLOTTING_MODULE.Plotter


class ConcretePlotter(Plotter):
    """Concrete test helper that exercises base-class method bodies."""

    def output_flux(self, results):
        """Call the abstract base implementation."""
        return Plotter.output_flux(self, results)

    def output_energy(self, results, energies):
        """Call the abstract base implementation."""
        return Plotter.output_energy(self, results, energies)

    def output_orderp(self, results, orderdata):
        """Call the abstract base implementation."""
        return Plotter.output_orderp(self, results, orderdata)

    def output_path(self, results, path_ensemble):
        """Call the abstract base implementation."""
        return Plotter.output_path(self, results, path_ensemble)

    def output_matched_probability(self, path_ensembles, detect, matched):
        """Call the abstract base implementation."""
        return Plotter.output_matched_probability(
            self, path_ensembles, detect, matched
        )


class IncompletePlotter(Plotter):
    """Intentionally incomplete subclass for abstract-class checks."""

    def output_flux(self, results):
        """Implement only one abstract method."""
        return results


class TestPlotter:
    """Tests for Plotter."""

    def test_plotter_is_abstract(self):
        """The base class and incomplete subclasses stay abstract."""
        with pytest.raises(TypeError):
            Plotter()

        with pytest.raises(TypeError):
            IncompletePlotter()

    @pytest.mark.parametrize(
        ('backup_value', 'expected'),
        [
            (True, True),
            ('yes', True),
            ('True', True),
            (False, False),
            ('no', False),
            ('false', False),
            (None, False),
        ],
    )
    def test_constructor_sets_attributes(self, backup_value, expected):
        """The constructor normalises the backup flag and stores metadata."""
        plotter = ConcretePlotter(
            backup=backup_value,
            plotter_type='unit-test',
            out_dir='plots',
        )

        assert plotter.backup is expected
        assert plotter.plotter_type == 'unit-test'
        assert plotter.out_dir == 'plots'

    def test_string_representation(self):
        """The string representation reports the plotter type."""
        plotter = ConcretePlotter(plotter_type='txt')

        assert str(plotter) == 'Plotter: txt'

    def test_base_method_bodies_return_none(self):
        """Concrete subclasses can call the abstract base methods."""
        plotter = ConcretePlotter()

        assert plotter.output_flux({}) is None
        assert plotter.output_energy({}, []) is None
        assert plotter.output_orderp({}, []) is None
        assert plotter.output_path({}, []) is None
        assert plotter.output_matched_probability([], False, {}) is None
