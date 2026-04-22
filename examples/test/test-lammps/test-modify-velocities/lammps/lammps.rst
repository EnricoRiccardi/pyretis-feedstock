LAMMPS velocity modification settings
======================================

Simulation settings
-------------------
task = md
steps = 1

Engine settings
---------------
class = lammps
lmp = lmp_serial
input_path = ../lammps_input
subcycles = 1

System settings
---------------
units = lj

Output
------
energy-file = 1
order-file = 1
trajectory-file = -1
screen = 1
