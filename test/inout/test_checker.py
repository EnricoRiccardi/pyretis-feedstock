# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the methods in pyretis.setup.checker"""
import os
import logging
import pytest
from io import StringIO
from pyretis.engines.internal import VelocityVerlet
from pyretis.core.particles import ParticlesExt
from pyretis.core.system import System
from pyretis.inout.checker import (check_for_bullshitt,
                                   check_engine,
                                   check_ensemble,
                                   check_interfaces)
from pyretis.inout.settings import fill_up_tis_and_retis_settings
from unittest.mock import patch

logging.disable(logging.CRITICAL)
HERE = os.path.abspath(os.path.dirname(__file__))

SETTINGS = {'simulation': {'task': 'retis',
                           'interfaces': [-1.0, -0.5, 0, 0.5, 1],
                           'exe_path': HERE,
                           'zero_ensemble': True,
                           'flux': True},
            'orderparameter': {'class': 'Position',
                               'dim': 'x',
                               'index': 0,
                               'periodic': False},
            'particles': {'type': 'internal'},
            'system': {'units': 'lj'},
            'engine': {'obj': VelocityVerlet(0.002)},
            'tis': {}}


def make_test_system(conf):
    """Just make a test system with particles."""
    system = System()
    system.particles = ParticlesExt(dim=3)
    system.particles.config = conf
    return system


class TestMethods:
    """Test some of the methods from .checker."""

    def test_check_ensemble(self):
        """Test check_ensemble."""
        play_set = SETTINGS.copy()
        with pytest.raises(ValueError) as err:
            assert not check_ensemble(play_set)
        assert 'No ensemble in settings' in str(err.value)

        fill_up_tis_and_retis_settings(play_set)
        play_set_base = play_set.copy()

        assert check_ensemble(play_set)

        play_set['simulation']['interfaces'][2] = 4
        with pytest.raises(ValueError) as err:
            assert not check_ensemble(play_set)
        assert ('The ensemble interface 0 is not present in the simulation'
                ' interface list') in str(err.value)

        # Not sorted ensembles
        play_set = play_set_base.copy()
        play_set['ensemble'][0]['interface'] = 4
        play_set['ensemble'][1]['interface'] = 3
        with pytest.raises(ValueError) as err:
            assert not check_ensemble(play_set)
        assert 'NOT properly sorted' in str(err.value)

        # Not declared ensembles
        play_set = play_set_base.copy()
        del play_set['ensemble'][0]['interface']
        with pytest.raises(ValueError) as err:
            assert not check_ensemble(play_set)
        assert 'without reference interface' in str(err.value)

        play_set = play_set_base.copy()
        del play_set['ensemble']
        with pytest.raises(ValueError) as err:
            assert not check_ensemble(play_set)
        assert 'No ensemble in settings' in str(err.value)

    def test_check_engine(self):
        """Test check_engine."""
        play_set = SETTINGS.copy()
        play_set['particles'] = {'type': 'external'}
        play_set['engine'] = {'type': 'external',
                              'cp2k': 'cp2k',
                              'input_path': 'No thank you',
                              'timestep': 1,
                              'cp2k_format': 'xyz',
                              'subcycles': 1}
        assert check_engine(play_set)

        del play_set['engine']['input_path']
        assert not check_engine(play_set)

        play_set['engine']['input_path'] = 'Hahaha'
        del play_set['engine']['cp2k_format']
        assert not check_engine(play_set)

        play_set['engine']['cp2k_format'] = 'xyz'
        del play_set['engine']['cp2k']
        play_set['engine']['gmx'] = 'Buuuuu'
        assert not check_engine(play_set)

    def test_check_interfaces(self):
        """Test check_interfaces."""
        settings = {'simulation': {'task': 'retis',
                                   'interfaces': [-1.0, -0.5, 0, 0.5, 1],
                                   'zero_ensemble': True,
                                   'flux': True}}

        assert check_interfaces(settings)
        settings['simulation']['zero_ensemble'] = False
        assert not check_interfaces(settings)

        settings['simulation']['zero_ensemble'] = True
        settings['simulation']['interfaces'] = [1, 2]
        assert not check_interfaces(settings)

        settings['simulation']['interfaces'] = [3, 2, 2.5]
        assert not check_interfaces(settings)

        settings['simulation']['interfaces'] = [-1.0, -0.5, 0, 0.5, 1]
        assert check_interfaces(settings)

    def test_check_for_bullshitt(self):
        """Test that _check for bullshitt finds the inconsistent settings."""
        # Insufficient interfaces for TIS/RETIS:
        settings = {'simulation': {'task': 'tis', 'interfaces': [1]}}
        with patch('sys.stdout', new=StringIO()):
            with pytest.raises(ValueError) as err:
                check_for_bullshitt(settings)
        assert 'Insufficient number of interfaces for tis' in str(err.value)
        # No references:
        settings = {'simulation': {'task': 'tis', 'interfaces': [1, 2, 3]},
                    'ensemble': [{'gino': 'strada'}, {'interface': 1},
                                 {'interface': 3}]}
        with patch('sys.stdout', new=StringIO()):
            with pytest.raises(ValueError) as err:
                check_for_bullshitt(settings)
        assert ('An ensemble has been introduced without references'
                ' (interface in ensemble settings)') in str(err.value)

        settings = {'simulation': {'task': 'tis', 'interfaces': [1, 2, 3]},
                    'ensemble': [{'ensemble_number': 0}, {'interface': 2},
                                 {'interface': 3}]}
        check_for_bullshitt(settings)

        # Not Sorted interfaces:
        settings = {'simulation': {'task': 'tis', 'interfaces': [2, 5, 1]}}
        with patch('sys.stdout', new=StringIO()):
            with pytest.raises(ValueError) as err:
                check_for_bullshitt(settings)
        assert ('Interface lambda positions in the simulation entry are NOT'
                ' sorted (small to large)') in str(err.value)

        # Wrong interfaces:
        settings = {'simulation': {'task': 'tis', 'interfaces': [1, 2, 3]},
                    'ensemble': [{'interface': 1}, {'interface': 2},
                                 {'interface': 4}]}
        with patch('sys.stdout', new=StringIO()):
            with pytest.raises(ValueError) as err:
                check_for_bullshitt(settings)
        assert 'in the simulation interface' in str(err.value)
