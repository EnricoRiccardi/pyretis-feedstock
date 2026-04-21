RETIS GROMACS sparse load (frames)
==================================

Simulation
----------
task = retis
steps = 5
interfaces = [-0.305, -0.285, -0.280, -0.270]
zero_left = -0.32
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
input_path = ../../../../shared_input/gromacs
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
load_and_kick = True

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
