RETIS test for LAMMPS
=====================

Simulation
----------
task = retis
steps = 10
interfaces = [0.5, 1.5, 3.5]

System
------
units = lj

Engine
------
class = lammps
lmp = lmp_serial
input_path = lammps_input
subcycles = 2

TIS
---
freq = 0.5
maxlength = 2000
aimless = True
allowmaxlength = False
zero_momentum = False
rescale_energy = False
sigma_v = -1
seed = 0

Initial-path
------------
method = kick
kick-from = previous

RETIS settings
--------------
swapfreq = 0.5
relative_shoots = None
nullmoves = True
swapsimul = True

Output settings
---------------
pathensemble-file = 1
screen = 10
order-file = 1
energy-file = 1
trajectory-file = 1
