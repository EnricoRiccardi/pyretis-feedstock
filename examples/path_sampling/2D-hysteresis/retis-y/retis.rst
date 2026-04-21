RETIS 2D hysteresis, y direction
================================

Simulation
----------
task = retis
steps = 1000000
interfaces = [-0.32, -0.21, -0.13, -0.08, -0.03, 0.03, 0.10, 0.40]

System
------
units = reduced
dimensions = 2
temperature = 0.4

Box
---
periodic = [False, False]

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

RETIS settings
--------------
swapfreq = 0.5
relative_shoots = None
nullmoves = True
swapsimul = True

Initial-path
------------
method = kick

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
description = 2D hysteresis

Potential
---------
class = Hyst2D
module = ../potential.py
parameter gamma1 = 1
parameter gamma2 = -10
parameter gamma3 = -10
parameter alpha1 = -30
parameter alpha2 = -3
parameter beta1 = -30
parameter beta2 = -3
parameter x0 = 0.2
parameter y0 = 0.4

Orderparameter
--------------
class = OrderX
module = ../order.py
dim = y
index = 0
periodic = False
inter_a = -0.30
inter_b = 0.37
energy_a = -9.0
energy_b = -9.0
