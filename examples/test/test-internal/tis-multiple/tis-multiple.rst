TIS 1D example
==============

Simulation
----------
task = make-tis-files
steps = 50
interfaces = [-0.9, -0.6, -0.3,  1.0]
zero_ensemble = True
seed = 16
rgen = 'rgen-borg'

System
------
units = reduced
dimensions = 1
temperature = 0.7
rgen = 'rgen-borg'

Box
---
periodic = [False]

Engine
------
class = Langevin
timestep = 0.05
gamma = 0.3
high_friction = False
seed = 16
rgen = 'rgen'

TIS
---
freq =  0.5
maxlength = 20000
aimless = True
allowmaxlength = False
zero_momentum = False
rescale_energy = False
sigma_v = -1
seed = 16
rgen = 'rgen-borg'

Initial-path
------------
method = kick

Particles
---------
position = {'input_file': 'initial.xyz'}
velocity = {'generate': 'maxwell',
            'momentum': False,
            'seed': 0, 'rgen': 'rgen'}
mass = {'Ar': 1.0}
name = ['Ar']
ptype = [0]

Forcefield
----------
description = 1D double well

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
energy-file = 100
order-file = 100
restart-file = 10
trajectory-file = 100
restart-file = 10
