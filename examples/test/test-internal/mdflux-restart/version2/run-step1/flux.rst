MD flux simulation
==================

Simulation settings
-------------------
task = md-flux
steps = 483
interfaces = [-0.9, -0.85, -0.80]

System settings
---------------
units = lj
dimensions = 1
temperature = 0.07

Engine settings
---------------
class = Langevin
timestep = 0.002
gamma = 0.3
high_friction = False
seed = 0

Particles
---------
position = {'input_file': '../initial.xyz'}
velocity = {'generate': 'maxwell',
            'momentum': False,
            'seed': 0}
mass = {'Ar': 1.0}
name = ['Ar']
ptype = [0]

ForceField settings
-------------------
description = 1D double well potential

Potential
---------
class = DoubleWell
a = 1.0
b = 2.0
c = 0.0

Orderparameter
--------------
class = Position
dim = x
index = 0
periodic = False

Output
------
backup = overwrite
energy-file = 1
order-file = 1
cross-file = 1
trajectory-file = -1
