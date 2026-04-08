# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the initiate kick family of methods."""
import collections
from io import StringIO
import logging
import pytest
from unittest.mock import patch
import os
import numpy as np
from pyretis.initiation.initiate_kick import (
    fix_path_by_tis,
    generate_initial_path_kick,
    _get_help,
    initiate_kick,
    initiate_kicki,
    initiate_kick_max,
)
from .help import (
    turn_on_logging,
    create_test_retis_simulation,
    create_test_tis_simulation,
)
logging.disable(logging.CRITICAL)
HERE = os.path.abspath(os.path.dirname(__file__))


class MockLog:
    """A mock log object for testing."""
    def __init__(self, caplog):
        self.output = caplog.messages


class TestGenericMethods:
    """Run tests for generic methods defined in the module."""

    def test_get_help(self):
        """Testing the get_help method."""
        interfaces = [-1, 0, 1]
        # Test what happens for unknown start conditions:
        with pytest.raises(ValueError):
            _get_help('r', interfaces)
        # Test what happens for incorrect number of interfaces:
        for i in (0, 1, 2, 4):
            with pytest.raises(ValueError):
                _get_help('R', [1]*i)
        # For correct input we expect to be given two methods back.
        imp_r, done_r = _get_help('R', interfaces)
        # Test that these two functions behave as expected:
        fake_order_value = collections.namedtuple(
            'orderparameter', ['ordermax', 'ordermin']
        )
        # Run some checks for improvement:
        current = fake_order_value([0.4, 0], [-1, 1])
        # Cases are:
        # 1. max > currentmax and min < middle.
        # 2. max > currentmax and min > middle.
        # 3. max < currentmax and min < middle.
        # 4. max < currentmax and min > middle.
        cases = [
            [(1.1, 0), (-1, 0)],
            [(1, 0), (0.5, 0)],
            [(-0.5, 0), (-1, 0)],
            [(0.3, 0), (0.1, 0)],
        ]
        results = [True, False, False, False]

        for case, res in zip(cases, results):
            newp = fake_order_value(case[0], case[1])
            assert imp_r(newp, current) == res
            assert done_r(newp) == res

        imp_l, done_l = _get_help('L', interfaces)
        # Cases are:
        # 1. min < currentmin and max > middle.
        # 2. min < currentmin and max < middle.
        # 3. min > currentmin and max < middle.
        # 4. min > currentmin and max > middle.
        current = fake_order_value([0.4, 0], [-1, 1])
        cases = [
            [(0.5, 0), (-2, 1)],
            [(-0.9, 0), (-1, 1)],
            [(-0.9, 0), (0, 1)],
            [(0.2, 0), (0, 1)],
        ]
        results = [True, False, False, False]
        for case, res in zip(cases, results):
            newp = fake_order_value(case[0], case[1])
            assert imp_l(newp, current) == res
            assert done_l(newp) == res


def _get_kick_initial(log):
    """Extract initial points for kicking from log messages."""
    kick = []
    for msg in log.output:
        if msg.find('Kicking from') != -1:
            kick.append(float(msg.split(':')[-1]))
    return kick


def _run_generate_initial(simulation):
    """Run the generate_initial_path_kick method for a given simulation."""
    for ens, sett in zip(simulation.ensembles,
                         simulation.settings['ensemble']):
        for attempt in generate_initial_path_kick(ens, sett['tis']):
            yield attempt


class TestInitiateKick:
    """Run tests for the initiate kick methods."""

    def _compare_generator_functions(self, gen1, gen2):
        """Might be over-the-top, but let's compare generator functions."""
        assert gen1.__name__ == gen2.__name__
        # Sure, let's us just compare the compiled bytecode:
        assert gen1.gi_code.co_code == gen2.gi_code.co_code

    def test_initiate_kick(self):
        """Test the initiate_kick selector method."""
        keywords = ('initial', 'previous', 'iNitiAL')
        methods = (initiate_kicki, initiate_kick_max, initiate_kicki)
        for key, method in zip(keywords, methods):
            settings = {
                'initial-path': {'kick-from': key}
            }
            init = initiate_kick(None, settings, cycle=0)
            gen = method(None, settings, cycle=0)
            self._compare_generator_functions(init, gen)
        # Test for unknown method:
        settings = {
            'initial-path': {'kick-from': 'Ekki-Ekki-Ekki-Ekki-Ptang'}
        }
        with pytest.raises(ValueError):
            initiate_kick(None, settings, cycle=0)

    def _help_with_initiation(self, simulation, method, correct_kick, caplog):
        """Run initiation and return points for kicking."""
        kick = []
        results = []
        with turn_on_logging():
            with patch('sys.stdout', new=StringIO()):
                caplog.clear()
                caplog.set_level(logging.INFO)
                for result in method(simulation,
                                     simulation.settings,
                                     cycle=0):
                    results.append(result)
                kick = _get_kick_initial(MockLog(caplog))
        assert len(kick) == len(correct_kick)
        for i, j in zip(kick, correct_kick):
            assert i == pytest.approx(j)
        return results

    def test_initiate_kick_max(self, caplog):
        """Test the initiate_kick_max method."""
        simulation = create_test_retis_simulation()
        correct_kick = [-1.0, -0.9016, -0.7663]
        results = self._help_with_initiation(simulation, initiate_kick_max,
                                             correct_kick, caplog)
        kick_re = [correct_kick[0]]
        for result in results:
            path = result[1]
            interfaces = result[-1].interfaces
            assert result[0]
            # Re-calculate the kicking point:
            order = np.array([i.order[0] for i in path.phasepoints])
            diff = order - interfaces[2]
            diff_neg = np.where(diff <= 0)[0]
            small_neg = np.argmax(diff[diff_neg])
            kick_re.append(order[diff_neg[small_neg]])
        # The last one is not used:
        kick_re.pop()
        assert len(kick_re) == len(correct_kick)
        for i, j in zip(correct_kick, kick_re):
            assert i == pytest.approx(j)
        # Check that is gets the right interface for ensemble 0
        simulation.ensembles[0]['path_ensemble'].ensemble_number = 0
        correct_kick = [-0.8893, -0.9016, -0.5203]
        self._help_with_initiation(simulation, initiate_kick_max,
                                   correct_kick, caplog)

    def test_initiate_kicki(self, caplog):
        """Test the initiate_kicki method."""
        simulation = create_test_retis_simulation()
        correct_kick = [-1.0, -1.0, -1.0]
        # Calculate path weight when under the WF HA ensemble:
        for i in range(len(simulation.settings['ensemble'])):
            simulation.settings['ensemble'][i]['tis']['high_accept'] = True
            simulation.settings['ensemble'][i]['tis']['shooting_move'] = 'wf'
        results = self._help_with_initiation(simulation, initiate_kicki,
                                             correct_kick, caplog)
        path_weights = [0.0, 5.0, 88.0]
        for i, result in enumerate(results):
            assert result[0]
            assert result[1].weight == path_weights[i]

    def test_no_0L_kicking(self):
        simulation = create_test_tis_simulation()
        ens = simulation.ensembles[0]
        # Initial system is -1 so set interfaces around that
        ens['path_ensemble'].interfaces = (-1.5, -0.5, -0.5)
        # Set start condition to LR
        ens['path_ensemble'].start_condition = ("R", "L")

        # Break kick_across_middle to make sure it is not called
        def f(*args):
            raise RuntimeError("THIS SHOULD NOT BE CALLED")
        ens['engine'].kick_across_middle = f
        path = next(
            generate_initial_path_kick(ens, simulation.settings['tis'])
        )

        assert path[0]

    def test_0L_kicking(self):
        simulation = create_test_tis_simulation()
        ens = simulation.ensembles[0]
        # Initial system is -1 so set all interfaces right of that
        ens['path_ensemble'].interfaces = (-0.5, 0, 0.5)
        # Set start condition to LR
        ens['path_ensemble'].start_condition = ("R", "L")

        # Break kick_across_middle to make sure it is not called
        def f(ensemble, middle, tis_settings):
            raise RuntimeError(middle)
        ens['engine'].kick_across_middle = f
        # Test if we see the right middle being called
        with pytest.raises(RuntimeError, match="-0.5"):
            next(generate_initial_path_kick(ens, simulation.settings['tis']))

    def test_0R_kicking(self):
        simulation = create_test_tis_simulation()
        ens = simulation.ensembles[0]
        # Initial system is -1 so set all interfaces left of that
        ens['path_ensemble'].interfaces = (-2.5, -2, -1.5)
        # Set start condition to LR
        ens['path_ensemble'].start_condition = ("R", "L")

        # Break kick_across_middle to make sure it is not called
        def f(ensemble, middle, tis_settings):
            raise RuntimeError(middle)
        ens['engine'].kick_across_middle = f
        # Test if we see the right middle being called
        with pytest.raises(RuntimeError, match="-1.5"):
            next(generate_initial_path_kick(ens, simulation.settings['tis']))

    def test_fix_path(self, caplog):
        """Test fix_path_by_tis."""
        simulation = create_test_tis_simulation()
        correct_kick = [-1.0]
        results = self._help_with_initiation(simulation, initiate_kicki,
                                             correct_kick, caplog)
        assert results[0][0]
        initial_path = results[0][1]
        # This initial path has already been accepted, but
        # try to improve it:
        new_path = fix_path_by_tis(
            initial_path,
            simulation.ensembles[0],
            simulation.settings['tis'],
        )
        # This should trigger the "did not improve" branch,
        # but it is still accepted. Test that we did not change
        # the path we gave in:
        assert initial_path is new_path
        # Make the path start end at wrong interface:
        interfaces = simulation.ensembles[0]['path_ensemble'].interfaces
        initial_path.phasepoints[0].order = [interfaces[-1] + 0.2, 0.0]
        initial_path.phasepoints[-1].order = [interfaces[-1] + 0.2, 0.0]
        check = initial_path.check_interfaces(interfaces)
        assert check[0] == 'R'
        assert check[1] == 'R'
        assert check[2] == 'M'
        assert check[3] == [False, True, True]
        # Try to improve this path:
        new_path = fix_path_by_tis(
            initial_path,
            simulation.ensembles[0],
            simulation.settings['tis'],
        )
        # The new path starts at L, ends at L and crosses M:
        assert initial_path is not new_path
        check = new_path.check_interfaces(interfaces)
        assert check[0] == 'L'
        assert check[1] == 'L'
        assert check[2] == 'M'
        # Check that we crossed the left interface, the middle one and back.
        assert check[3] == [True, True, False]
        # Check that the path is accepted:
        assert new_path.status == 'ACC'

    def test_generate_initial_fail_forward(self):
        """Test forward failure of the generate_initial_path_kick method.

        Here, we test how the generate_initial_path_kick method behaves
        when the forward generation fails.

        """
        # Make forward fail, by limiting to a very short path:
        simulation = create_test_tis_simulation(maxlength=3)
        state = simulation.rgen.get_state()
        for i, attempt in enumerate(_run_generate_initial(simulation)):
            assert not attempt[0]
            assert attempt[-1] is None
            assert attempt[1].startswith('Forward path failed:')
            simulation.rgen.set_state(state)
            simulation.ensembles[0]['engine'].reset()
            if i >= 1:
                # Just do two iterations so that we catch the continue
                # statement after the yield.
                break

        # Test exp initialization
        simulation.settings['ensemble'][0]['tis']['shooting_move'] = 'exp'
        for i, attempt in enumerate(_run_generate_initial(simulation)):
            assert attempt[0]
            break

    def test_generate_initial_fail_backward(self):
        """Test backward failure of the generate_initial_path_kick method.

        Here, we test how the generate_initial_path_kick method behaves
        when the backward generation fails.

        """
        # Make backward fail, by using a special engiene.
        simulation = create_test_tis_simulation(
            engine_type='MockEngineVelocitySupremacist',
            maxlength=50,
        )
        state = simulation.rgen.get_state()
        for i, attempt in enumerate(_run_generate_initial(simulation)):
            assert not attempt[0]
            assert attempt[-1] is None
            assert attempt[1].startswith('Backward path failed:')
            simulation.rgen.set_state(state)
            if i >= 1:
                # Just do two iterations so that we catch the continue
                # statement after the yield.
                break

    def test_generate_initial_fail_length(self):
        """Test path paste failure of generate_initial_path_kick.

        Here, we test how the generate_initial_path_kick method behaves
        when the merging of the forward and backward paths fail.

        """
        # Make path paste fail, by limiting to a short path:
        simulation = create_test_tis_simulation(maxlength=60)
        state = simulation.rgen.get_state()
        for i, attempt in enumerate(_run_generate_initial(simulation)):
            assert not attempt[0]
            assert attempt[-1] is None
            assert attempt[1] == 'Initial path was too long.'
            simulation.rgen.set_state(state)
            simulation.ensembles[0]['engine'].reset()
            if i >= 1:
                # Just do two iterations so that we catch the continue
                # statement after the yield.
                break

    def test_generate_initial_fail_other(self):
        """Test generate_initial_path_kick for start/end conditions.

        Here, we test how the generate_initial_path_kick method behaves
        when the start/end conditions of the path ensemble are set
        incorrectly.

        """
        # Make forward fail, by limiting to a short path:
        simulation = create_test_tis_simulation(maxlength=100)
        simulation.ensembles[0]['path_ensemble'].start_condition = 'tiktok'
        state = simulation.rgen.get_state()
        for i, attempt in enumerate(_run_generate_initial(simulation)):
            assert not attempt[0]
            assert attempt[-1] is None
            assert attempt[1] == 'Could not generate initial path will retry!'
            # Reset rgen and the engine that we just do the same again:
            simulation.rgen.set_state(state)
            simulation.ensembles[0]['engine'].reset()
            if i >= 1:
                # Just do two iterations so that we catch the continue
                # statement after the yield.
                break

    def test_generate_initial_fix(self):
        """Test generate_initial_path_kick with TIS fixes.

        Here, we test how this method behaves when we have a path
        that starts and ends at the wrong interface.

        """
        # Make a path that start and ends at the wrong interface,
        # by using an engine which is always increasing the order
        # parameter:
        simulation = create_test_tis_simulation(
            engine_type='MockEngineOneWay'
        )
        for i, attempt in enumerate(_run_generate_initial(simulation)):
            if i == 0:
                assert not attempt[0]
                assert attempt[-1] is None
                assert attempt[1] == 'Trying to fix path by TIS moves'
                # For this to actually work, we need to allow the engine
                # to go backwards as well:
                simulation.ensembles[0]['engine'].timestep =\
                    -simulation.ensembles[0]['engine'].timestep
            elif i == 1:
                assert attempt[0]
                assert attempt[-1] is not None
                # We will end here, as the fix_path_by_tis is tested elsewhere.
            else:
                raise ValueError('Wrong number of iterations.')
