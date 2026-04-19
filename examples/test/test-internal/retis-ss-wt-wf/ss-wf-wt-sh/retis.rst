Retis 1D example
================

Simulation
----------
task = retis
steps = 15
interfaces = [-0.99, -0.7, -0.6, -0.5] 

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
timestep = 0.05
gamma = 0.2
high_friction = False
seed = 0

TIS settings
------------
freq = 0.0
maxlength = 10000
aimless = True
allowmaxlength = False
zero_momentum = False
rescale_energy = False
sigma_v = -1
seed = 0
shooting_moves = ['ss','wf','wt','sh']
n_jumps = 6
interface_sour = -0.80
high_accept = False

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
position = {'input_file': 'initial.xyz'}
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
backup = overwrite
order-file = 1
trajectory-file = -1
energy-file = 1
