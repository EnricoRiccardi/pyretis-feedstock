Retis 1D example
================

Simulation
----------
task = retis
steps = 12
interfaces = [-0.9, -0.6, -0.3, 1.0]

System
------
units = reduced
dimensions = 1
temperature = 0.7

Box
---
periodic = [False]

Engine
------
class = Langevin
timestep = 0.02
gamma = 0.3
high_friction = False
seed = 0

TIS settings
------------
freq = 0.5
maxlength = 5000
aimless = True
allowmaxlength = False
zero_momentum = False
rescale_energy = False
sigma_v = -1
seed = 0

RETIS settings
--------------
swapfreq = 0.5
relative_shoots = None
nullmoves = True
swapsimul = True

Initial-path settings
---------------------
method = kick
kick-from = initial

Particles
---------
position = {'input_file': '../initial.xyz'}
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
trajectory-file = 1
energy-file = 1
order-file = 1
