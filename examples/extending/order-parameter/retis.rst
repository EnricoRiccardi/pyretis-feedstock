RETIS 1D example - custom order parameter
=========================================

Simulation
----------
task = retis
steps = 100
interfaces = [-1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, 0.0]

System
------
units = reduced
dimensions = 1
temperature = 0.07

Box
---
periodic = [False]

Engine
------
class = Langevin
timestep = 0.002
gamma = 0.3
high_friction = False
seed = 0

TIS settings
------------
freq =  0.5
maxlength = 20000
aimless = True
allowmaxlength = False
zero_momentum = False
rescale_energy = False
sigma_v =  -1
seed = 0
initial_path = kick

Initial-path settings
---------------------
method = kick
kick-from = initial

RETIS settings
--------------
swapfreq = 0.5
relative_shoots = None
nullmoves = True
swapsimul = True

Particles
---------
position = {'file': 'initial.xyz'}
velocity = {'generate': 'maxwell',
            'momentum': False,
            'seed': 0}
mass = {'Ar': 1.0}
name = ['Ar']
type = [0]

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
class = PlaneDistanceX
module = orderparameter1.py
index = 0
plane_position = 1.0

Output
------
trajectory-file = 100
energy-file = 100
order-file = 100
write_vel = True
