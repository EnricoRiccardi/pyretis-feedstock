RETIS GROMACS restart test
==========================

Simulation
----------
task = retis
steps = 5
interfaces = [-0.26, -0.24, -0.22, -0.18]
rgen = 'rgen-borg'
seed = 0

System
------
units = gromacs

Engine settings
---------------
class = GromacsEngine2R
module = gromacs.py
gmx = GMXCOMMAND
mdrun = GMXCOMMAND mdrun
input_path = ../../../../shared_input/gromacs
timestep = 0.002
subcycles = 5
gmx_format = g96

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
method = kick
kick_from = previous

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
backup = 'append'
pathensemble-file = 1
screen = 1
order-file = 1
energy-file = 1
trajectory-file = 1
restart-file = 1
