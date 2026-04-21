RETIS CP2K load test
====================

Simulation
----------
task = retis
steps = 5
interfaces = [0.8, 0.85, 0.9, 0.95]
rgen = 'rgen'
seed = 0

System
------
units = cp2k
temperature = 500

Engine settings
---------------
class = cp2k
cp2k = cp2k
input_path = cp2k_input
timestep = 0.5
subcycles = 5
seed = 0

TIS settings
------------
freq = 0.5
maxlength = 20000
aimless = True
allowmaxlength = False
zero_momentum = False
rescale_energy = False
sigma_v = -1
rgen = 'rgen'
seed = 0

Initial-path
------------
method = load
load_folder = pippo
load_and_kick = True

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

Collective-variable
-------------------
class = Distance
index = (0, 1)

Output settings
---------------
backup = overwrite
pathensemble-file = 1
screen = 10
order-file = 1
energy-file = 1
trajectory-file = 1
