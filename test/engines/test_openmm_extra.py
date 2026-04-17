import unittest
from unittest.mock import MagicMock, patch
import numpy as np
import pyretis.engines.openmm as openmm_mod
from pyretis.engines.openmm import OpenMMEngine


class TestOpenMMExtra(unittest.TestCase):
    def test_openmm_not_installed(self):
        with patch('pyretis.engines.openmm.HAS_OPENMM', False):
            with self.assertRaises(RuntimeError) as cm:
                OpenMMEngine('sim')
            self.assertEqual(str(cm.exception), "OpenMM is not installed")

    def test_init_invalid_string(self):
        # OpenMM must be installed for this test to pass the first check
        with patch('pyretis.engines.openmm.HAS_OPENMM', True):
            with self.assertRaises(RuntimeError) as cm:
                OpenMMEngine('invalid_sim', openmm_module=None)
            self.assertIn("could not be loaded from module",
                          str(cm.exception))

    def test_clean_up_and_dump(self):
        with patch('pyretis.engines.openmm.HAS_OPENMM', True):
            engine = OpenMMEngine(MagicMock())
            # These should not raise and do nothing
            self.assertIsNone(engine.clean_up())
            self.assertIsNone(engine.dump_phasepoint(None))

    def test_modify_velocities_rescale_negative(self):
        with patch('pyretis.engines.openmm.HAS_OPENMM', True):
            engine = OpenMMEngine(MagicMock())
            system = MagicMock()
            system.particles.vpot = 10.0
            system.particles.pos = np.zeros((1, 3))
            system.particles.vel = np.zeros((1, 3))
            ensemble = {'system': system, 'rgen': MagicMock()}

            with patch('pyretis.engines.openmm.calculate_kinetic_energy',
                       return_value=(5.0, None)):
                dek, kin = engine.modify_velocities(
                    ensemble, {'rescale': -1.0})
                self.assertEqual(dek, 0.0)
                self.assertEqual(kin, 5.0)

    def test_modify_velocities_no_aimless(self):
        with patch('pyretis.engines.openmm.HAS_OPENMM', True):
            engine = OpenMMEngine(MagicMock())
            system = MagicMock()
            system.particles.pos = np.array([[0, 0, 0]])
            system.particles.vel = np.array([[1, 1, 1]])
            rgen = MagicMock()
            rgen.draw_maxwellian_velocities.return_value = (
                np.array([[0.1, 0.1, 0.1]]), None)
            ensemble = {'system': system, 'rgen': rgen}

            vel_settings = {'aimless': False, 'sigma_v': 0.1, 'momentum': True}

            with patch('pyretis.engines.openmm.calculate_kinetic_energy',
                       return_value=(5.0, None)):
                with patch('pyretis.engines.openmm.reset_momentum') as mock_r:
                    engine.modify_velocities(ensemble, vel_settings)
                    np.testing.assert_allclose(
                        system.particles.vel,
                        np.array([[1.1, 1.1, 1.1]]))
                    mock_r.assert_called_once()

    def test_propagate_reverse(self):
        with patch('pyretis.engines.openmm.HAS_OPENMM', True):
            # Mock simulation and context
            mock_sim = MagicMock()
            mock_context = MagicMock()
            mock_sim.context = mock_context

            engine = OpenMMEngine(mock_sim)
            engine.integration_step = MagicMock()
            engine.calculate_order = MagicMock(return_value=[0.5])

            path = MagicMock()
            path.maxlen = 10
            system = MagicMock()
            system.particles.pos = np.zeros((1, 3))
            system.particles.vel = np.array([[1.0, 1.0, 1.0]])
            system.box = MagicMock()

            ensemble = {'system': system, 'interfaces': (0, 0.5, 1.0)}

            with patch('pyretis.engines.openmm.calculate_kinetic_energy',
                       return_value=(5.0, None)):
                with patch.object(OpenMMEngine, 'add_to_path') as mock_add:
                    mock_add.side_effect = [
                        ('status', False, False, None),
                        ('status', False, True, None)
                    ]
                    OpenMMEngine.propagate(
                        engine, path, ensemble, reverse=True)
                    self.assertEqual(engine.integration_step.call_count, 1)


if __name__ == '__main__':
    unittest.main()
