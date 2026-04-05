Molecular dynamics example settings
===================================

Simulation settings
-------------------
task = md
steps = 10

Engine settings
---------------
class = gromacs2
gmx = GMXCOMMAND
mdrun = GMXCOMMAND mdrun
input_path = ../../../shared_input/gromacs
timestep = 0.002
subcycles = 5
gmx_format = g96

System settings
---------------
units = gromacs

Output
------
energy-file = 1
order-file = 1
trajectory-file = -1
screen = 1
