# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Example of interaction with a Simulation."""
from pyretis.core import System
from pyretis.simulation import Simulation


def function(sys):
    """Return the set temperature for the system."""
    return sys.temperature


simulation = Simulation(settings={}, controls={'steps': 10})
system = System(temperature=300)

my_task = {'func': function,
           'args': [system],
           'first': True,
           'result': 'temperature'}

simulation.add_task(my_task)

for result in simulation.run():
    step = result['cycle']['step']
    temp = result['temperature']['set']
    print(step, temp)
