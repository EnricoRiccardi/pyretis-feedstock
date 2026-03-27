# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the lattice generation tools in pyretis.tools"""
import logging
import pytest
import numpy as np
from pyretis.tools.lattice import generate_lattice
logging.disable(logging.CRITICAL)


class TestLattice:
    """Test the lattice generation in pyretis.tools.lattice."""

    def test_lattice_generation(self):
        """Test that we can generate all lattices."""
        xyz = {'fcc': np.array([[0., 0., 0.], [1., 1., 0.],
                                [0., 1., 1.], [1., 0., 1.]]),
               'bcc': np.array([[0., 0., 0.], [0.5, 0.5, 0.5],
                                [1., 0., 0.], [1.5, 0.5, 0.5]]),
               'diamond': np.array([[0., 0., 0.],
                                    [0., 0.705, 0.705],
                                    [0.705, 0., 0.705],
                                    [0.705, 0.705, 0.],
                                    [0.3525, 0.3525, 0.3525],
                                    [0.3525, 1.0575, 1.0575],
                                    [1.0575, 0.3525, 1.0575],
                                    [1.0575, 1.0575, 0.3525]]),
               'hcp': np.array([[0., 0., 0.], [0.45, 0.45, 0.],
                                [0.45, 0.75, 0.45], [0., 0.3, 0.45]]),
               'sc': np.array([[0., 0., 0.],
                               [0., 0., 3.14159],
                               [0., 0., 6.28318],
                               [0., 0., 9.42477]]),
               'sq': np.array([[0., 0.], [0., 3.], [3., 0.], [3., 3.]]),
               'sq2': np.array([[0., 0.], [1.5, 1.5], [0., 3.], [1.5, 4.5],
                                [3., 0.], [4.5, 1.5], [3., 3.], [4.5, 4.5]])}

        size = {'fcc': np.array([[0.0, 2.0], [0.0, 2.0], [0.0, 2.0]]),
                'bcc': np.array([[0.0, 2.0], [0.0, 1.0], [0.0, 1.0]]),
                'diamond': np.array([[0.0, 1.41], [0.0, 1.41], [0.0, 1.41]]),
                'hcp': np.array([[0.0, 0.9], [0.0, 0.9], [0.0, 0.9]]),
                'sc': np.array([[0.0, 3.14159], [0.0, 3.14159],
                                [0.0, 12.56636]]),
                'sq': np.array([[0.0, 6], [0.0, 6]]),
                'sq2': np.array([[0.0, 6], [0.0, 6]])}

        args = {'bcc': [2, 1, 1],
                'fcc': [1, 1, 1],
                'hcp': [1, 1, 1],
                'diamond': [1, 1, 1],
                'sc': [1, 1, 4],
                'sq': [2, 2],
                'sq2': [2, 2]}

        kwargs = {'bcc': dict(lcon=1.0),
                  'fcc': dict(lcon=2.0),
                  'hcp': dict(lcon=0.9),
                  'sc': dict(lcon=3.14159),
                  'sq': dict(lcon=3),
                  'sq2': dict(lcon=3),
                  'diamond': dict(lcon=1.41)}
        # Test bcc:
        for key in args:
            xyz_g, size_g = generate_lattice(key, args[key], **kwargs[key])
            assert np.allclose(xyz_g, xyz[key], atol=1.0e-8)
            assert np.allclose(np.array(size_g), size[key], atol=1.0e-8)
        # Test if we give almost nothing
        pos, size = generate_lattice(lattice='1d')
        assert pos[0] == 0
        assert size == [[0.0, 1]]

    def test_generate_no_lcon_or_dens(self):
        """Test if the generate lattice methods fails when it should."""
        # Test if we feed it some unknown lattice:
        args = ['undefined lattice', None]
        kwargs = dict(lcon=None, density=None)
        with pytest.raises(ValueError):
            generate_lattice(*args, **kwargs)
        # Test if we feed it some strange arguments:
        args = [['undefined', 'lattice'], None]
        with pytest.raises(ValueError):
            generate_lattice(*args, **kwargs)
        # Test if we do not give lattice parameters:
        args = ['sc', [1, 1, 1]]
        with pytest.raises(ValueError):
            generate_lattice(*args)
        # Test if we give to few repeats
        args = ['sc', [1, 1]]
        kwargs = dict(lcon=1.0)
        with pytest.raises(ValueError):
            generate_lattice(*args, **kwargs)
