# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the UmbrellaWindowSimulation class."""
import logging
import pytest
import numpy as np
from pyretis.core.units import create_conversion_factors
from pyretis.core.box import create_box
from pyretis.core.system import System
from pyretis.core.random_gen import RandomGenerator
from pyretis.core.particles import Particles
from pyretis.forcefield import ForceField
from pyretis.forcefield.potentials import DoubleWell, RectangularWell
from pyretis.simulation.mc_simulation import UmbrellaWindowSimulation
logging.disable(logging.CRITICAL)


def create_test_system():
    """Just set up and create a simple test system."""
    create_conversion_factors('reduced')
    box = create_box(periodic=[False])
    system = System(temperature=500, units='reduced', box=box)
    system.particles = Particles(dim=system.get_dim())
    system.add_particle(name='X', pos=np.array([-0.7]))
    potential_dw = DoubleWell(a=1, b=1, c=0.02)
    potential_rw = RectangularWell()
    system.forcefield = ForceField(
        'Double well with rectangular bias',
        potential=[potential_dw, potential_rw],
        params=[{'a': 1.0, 'b': 1.0, 'c': 0.02}, None]
    )
    return system


class TestUmbrellaWindow:
    """Run the tests for UmbrellaWindow Simulation class."""

    def test_umbrellawindow(self):
        """Test that we can create the simulation object."""
        system = create_test_system()
        rgen = RandomGenerator()
        potential_rw = system.forcefield.potential[1]
        potential_rw.set_parameters({'left': -1.0, 'right': -0.4})
        system.potential()
        ensemble = {'rgen': rgen, 'system': system}
        settings = {'simulation': {'rgen': 'rgen',
                                   'umbrella': [-1.0, -0.4],
                                   'overlap': -0.5,
                                   'maxdx': 0.1}}
        controls = {'mincycle': 10}
        simulation = UmbrellaWindowSimulation(ensemble, settings, controls)
        for _ in simulation.run():
            pass
        assert system.particles.pos[0][0] > -0.5
        # Test restart features:
        restart = simulation.restart_info()
        assert restart['type'] == simulation.simulation_type
        settings['restart'] = restart
        settings['restart']['overlap'] = 'NonSense'
        simulation.cycle['step'] = 999
        simulation = UmbrellaWindowSimulation(ensemble, settings, controls)
        assert simulation.overlap == 'NonSense'
        assert simulation.cycle['step'] == 17
        # Test the creation when the random generator is not given:
        ensemble = {'system': system}
        simulation = UmbrellaWindowSimulation(
            ensemble,
            settings={'simulation': {'umbrella': [-1.0, -0.4],
                                     'overlap': -0.5,
                                     'maxdx': 0.1}},
            controls={})
        assert isinstance(simulation.rgen, type(rgen))

        with pytest.raises(TypeError):
            simulation = UmbrellaWindowSimulation(ensemble,
                                                  settings=None,
                                                  controls={})
