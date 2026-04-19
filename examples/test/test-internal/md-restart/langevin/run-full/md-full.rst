Molecular dynamics example settings
===================================

Simulation settings
-------------------
task = md-nve
steps = 100

Engine settings
---------------
class = Langevin
timestep = 0.02
gamma = 0.3 
seed = 0
high_friction = False

System settings
---------------
temperature = 2.0
units = lj

Particles
---------
position = {'generate': 'fcc',
            'repeat': [2, 2, 2],
            'density': 0.9}

velocity = {'generate': 'maxwell',
            'temperature': 2.0,
            'momentum': True,
            'seed': 0}

mass = {'Ar': 1.0}
name = ['Ar']
ptype = [0]

Forcefield settings
-------------------
description = Lennard Jones test

Potential
---------
class = PairLennardJonesCutnp
shift = True
parameter 0 = {'sigma': 1.0, 'epsilon': 1.0, 'rcut': 2.5}

Output
------
prefix = md-full-
backup = overwrite
energy-file = 1
order-file = 1
cross-file = 1
trajectory-file = 1
restart-file = 1
