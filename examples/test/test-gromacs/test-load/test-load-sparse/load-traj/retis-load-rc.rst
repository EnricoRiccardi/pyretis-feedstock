RETIS GROMACS sparse load (trajectory)
======================================

Simulation
----------
task = retis
steps = 5
interfaces = [2.19, 2.22, 2.63]
zero_left = 2.16
rgen = 'rgen-borg'
seed = 0

System
------
units = gromacs

Engine settings
---------------
class = gromacs2
gmx = GMXCOMMAND
mdrun = GMXCOMMAND mdrun
input_path = gromacs_input
timestep = 0.002
subcycles = 5
gmx_format = gro

TIS settings
------------
freq = 0.5
maxlength = 20000
aimless = True
allowmaxlength = False
zero_momentum = False
rescale_energy = False
sigma_v = -1
rgen = 'rgen-borg'
seed = 0

Initial-path
------------
method = load
load_folder = pippo

RETIS settings
--------------
swapfreq = 0.5
relative_shoots = None
nullmoves = True
swapsimul = True

Orderparameter
--------------
class = Distance
module = orderp.py
index = [0, 3]

Output settings
---------------
pathensemble-file = 1
screen = 10
order-file = 1
energy-file = 1
trajectory-file = 1
