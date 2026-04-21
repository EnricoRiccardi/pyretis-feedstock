Retis 2-Ar LJ example (OpenMM)
==============================

Simulation
----------
task = retis
steps = 5
interfaces = [0.30, 0.35, 0.40]

System
------
units = gromacs

Engine
------
type = openmm
class = openmm
openmm_simulation = simulation
openmm_module = openmm_sim.py
subcycles = 10

TIS settings
------------
freq = 0.5
maxlength = 20000
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

Orderparameter
--------------
class = Distance
index = (0, 1)
periodic = False

Initial-path
------------
method = kick
kick-from = previous

Output
------
backup = overwrite
pathensemble-file = 1
screen = 10
order-file = 1
energy-file = 1
trajectory-file = 1
