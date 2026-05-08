PPTIS 1D example
================

Simulation
----------
task = repptis
steps = 16  #1000
interfaces = [0., 0.1, 0.2, 0.3, 0.4]
zero_left = -0.1
permeability = True

System
------
units = reduced
dimensions = 1
temperature = 1.

Box
---
periodic = [True]
low = [-1.]
high = [10.]

Engine
------
class = Langevin
timestep = 0.0002
gamma = 5.
high_friction = True
seed = 0

TIS
---
freq =  0.
maxlength = 200000
aimless = True
allowmaxlength = False
zero_momentum = False
rescale_energy = False
sigma_v = -1
seed = 0

RETIS settings
--------------
swapfreq = 0.75
relative_shoots = None
#nullmoves = True
nullmoves = False
swapsimul = True

Initial-path
------------
method = kick
kick-from = previous

Particles
---------
position = {'input_file': 'initial.xyz'}
velocity = {'generate': 'maxwell',
            'momentum': False,
            'seed': 0}
mass = {'Ar': 1.0}
name = ['Ar']
type = [0]

Forcefield
----------
description = 1D flat, no force

Orderparameter
--------------
class = Position
dim = x
index = 0
periodic = True

Collective-variable
-------------------
class = Velocity
dim = x
index = 0

Potential
---------
class = FlatPotential1D
module = flat-potential.py

Output
------
screen = 1000
energy-file = 1000
order-file = 1

