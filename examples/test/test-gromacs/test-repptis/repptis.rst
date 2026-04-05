RETIS test for GROMACS
======================

Simulation
----------
task = repptis
steps = 16
interfaces = [-0.26, -0.24, -0.22, -0.20]
zero_left = -0.28
permeability = True

System
------
units = gromacs

Engine
------
class = GromacsEngine2R
module = gromacs.py
gmx = GMXCOMMAND
mdrun = GMXCOMMAND mdrun -ntmpi 1 -ntomp 1
input_path = ../../shared_input/gromacs
timestep = 0.002
subcycles = 10
gmx_format = gro

TIS
---
freq =  0.0
maxlength = 20000
aimless = True
allowmaxlength = False
zero_momentum = False
rescale_energy = False
sigma_v = -1
seed = 0

Initial-path
------------
method = load
load_folder = load
#load_and_kick = True

RETIS settings
--------------
swapfreq = 0.5
relative_shoots = None
nullmoves = True
swapsimul = True

Orderparameter
--------------
class = RingDiffusion
module = orderp.py

Output settings
---------------
pathensemble-file = 1
screen = 10
order-file = 1
energy-file = 1
trajectory-file = 1
