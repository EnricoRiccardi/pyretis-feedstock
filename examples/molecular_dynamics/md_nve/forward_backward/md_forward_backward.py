# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""
Example of running a MD NVE simulation.

In this example, we are running backward and forward using
the internal engines/force fields.
"""

# pylint: disable=invalid-name
import sys
from matplotlib import pyplot as plt
from matplotlib import gridspec
from pyretis.core.units import create_conversion_factors
from pyretis.setup import create_simulation
from pyretis.inout.fileio import FileIO
from pyretis.inout.formats import ThermoTableFormatter
# for plotting:
from pyretis.inout.plotting import mpl_set_style
# simulation settings:
settings = {}
settings['simulation'] = {
    'task': 'md-nve',
    'steps': 30,
}
settings['system'] = {
    'units': 'lj',
    'temperature': 2.5,
    'dimensions': 2
}
settings['engine'] = {
    'class': 'velocityverlet',
    'timestep': 0.002
}
settings['output'] = {
    'backup': 'overwrite',
    'energy-file': 1,
    'screen': 10,
    'trajectory-file': 1
}
settings['potential'] = [
    {'class': 'PairLennardJonesCutnp',
     'parameter': {0: {'sigma': 1, 'epsilon': 1, 'factor': 2.5},
                   1: {'sigma': 1, 'epsilon': 1, 'factor': 2.5}},
     'shift': True, 'dim': 2}
]
settings['particles'] = {
    'position': {'input_file': 'initial.gro'},
    'velocity': {'generate': 'maxwell', 'momentum': True, 'seed': 0}
}

create_conversion_factors(settings['system']['units'])
simulation_nve = create_simulation(settings)

# Set up extra output:
thermo_file = FileIO('thermo.txt', 'w', ThermoTableFormatter(), backup=False)
thermo_file.open()
store_results = []
# Also create some other outputs:
simulation_nve.set_up_output(settings, progress=False)
# Run the simulation forward:
for result in simulation_nve.run():
    stepno = result['cycle']['stepno']
    thermo_file.output(stepno, result['thermo'])
    result['thermo']['stepno'] = stepno
    store_results.append(result['thermo'])
# Run backward:
simulation_nve.ensemble['system'].particles.vel *= -1.0
simulation_nve.extend_cycles(settings['simulation']['steps'] - 1)
for result in simulation_nve.run():
    stepno = result['cycle']['stepno']
    thermo_file.output(stepno, result['thermo'])
    result['thermo']['stepno'] = stepno
    store_results.append(result['thermo'])
thermo_file.close()

mpl_set_style()  # Load PyRETIS plotting style
step = [res['stepno'] for res in store_results]
pot_e = [res['vpot'] for res in store_results]
kin_e = [res['ekin'] for res in store_results]
tot_e = [res['etot'] for res in store_results]
pressure = [res['press'] for res in store_results]
temp = [res['temp'] for res in store_results]
# Plot some energies:
fig1 = plt.figure()
gs = gridspec.GridSpec(2, 2)
ax1 = fig1.add_subplot(gs[:, 0])
ax1.plot(step, pot_e, label='Potential')
ax1.plot(step, kin_e, label='Kinetic')
ax1.plot(step, tot_e, label='Total')
ax1.set_xlabel('Step no.')
ax1.set_ylabel('Energy per particle')
ax1.legend(loc='center left', prop={'size': 'small'})
ax2 = fig1.add_subplot(gs[0, 1])
ax2.plot(step, temp)
ax2.set_ylabel('Temperature')
ax3 = fig1.add_subplot(gs[1, 1])
ax3.plot(step, pressure)
ax3.set_xlabel('Step no.')
ax3.set_ylabel('Pressure')
fig1.tight_layout()
if 'noplot' not in sys.argv[1:]:
    plt.show()
