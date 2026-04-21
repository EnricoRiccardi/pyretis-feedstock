RETIS 1D explorer
=================

Simulation
----------
task = explore
steps = 16
interfaces = [-1.9, -0.7, -0.5, 0.3]
rgen = 'rgen'
seed = 0

System
------
units = reduced
dimensions = 1
temperature = 0.7
rgen = 'rgen'

Box
---
periodic = [False]

Engine
------
class = Langevin
timestep = 0.02
gamma = 0.3
high_friction = False
rgen = 'rgen'
seed = 0

TIS settings
------------
freq = 0.5
maxlength = 100
aimless = True
allowmaxlength = False
zero_momentum = False
rescale_energy = False
sigma_v = -1
rgen = 'rgen'
seed = 0

Initial-path settings
---------------------
method = load
load_folder = pippo

Particles
---------
position = {'input_file': 'pippo/initial.xyz'}
velocity = {'generate': 'maxwell',
            'momentum': False,
            'seed': 0}
mass = {'Ar': 1.0}
name = ['Ar']
ptype = [0]

Forcefield settings
-------------------
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
trajectory-file = 50
energy-file = 1
order-file = 1
