Molecular dynamics example settings
===================================

Simulation
----------
task = md-flux
steps = 5
interfaces = [2.0]

Engine
------
class = cp2k
cp2k = cp2k
input_path = ../../shared_input/cp2k
timestep = 0.5
subcycles = 3
conf = h2.xyz
cp2k_format = xyz

System settings
---------------
units = cp2k
temperature = 500

Orderparameter
--------------
class = PositionVelocity
dim = x
index = 0
periodic = False

Output
------
energy-file = 1
order-file = 1
trajectory-file = -1
screen = 1
