Retis 2-Ar LJ example (OpenMM)
==============================

Simulation settings
-------------------
task = 'retis'
steps = 5
interfaces = [0.3, 0.35, 0.4]
exe_path = '/home/enrico/Research/pyretis/examples/test/test-pyvisa/test-openmm-pyvisa'
rgen = 'rgen'
priority_shooting = False
flux = True
zero_ensemble = True
endcycle = 0

System settings
---------------
units = 'gromacs'
dimensions = 3
temperature = 1.0

Engine settings
---------------
type = 'openmm'
class = 'openmm'
openmm_simulation = 'simulation'
openmm_module = 'openmm_sim.py'
subcycles = 10
exe_path = '/home/enrico/Research/pyretis/examples/test/test-pyvisa/test-openmm-pyvisa'
rgen = 'rgen'

Particles settings
------------------
type = 'internal'

Orderparameter settings
-----------------------
class = 'Distance'
index = (0, 1)
periodic = False
name = 'Order Parameter'

Output settings
---------------
backup = 'overwrite'
pathensemble-file = 1
screen = 10
order-file = 1
energy-file = 1
trajectory-file = 1
cross-file = 1
restart-file = 1

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
high_accept = False
rgen = 'rgen'
shooting_move = 'sh'
shooting_moves = []
mirror_freq = 0
target_freq = 0
target_indices = []

Initial-path settings
---------------------
method = 'kick'
kick-from = 'previous'

RETIS settings
--------------
swapfreq = 0.5
relative_shoots = None
nullmoves = True
swapsimul = True

Ensemble
--------
heading = {'text': 'Retis 2-Ar LJ example (OpenMM)\n=============================='}
simulation task = retis
simulation steps = 5
simulation interfaces = [-inf, 0.3, 0.3]
simulation exe_path = /home/enrico/Research/pyretis/examples/test/test-pyvisa/test-openmm-pyvisa
simulation rgen = rgen
simulation priority_shooting = False
simulation flux = True
simulation zero_ensemble = True
system units = gromacs
system dimensions = 3
system temperature = 1.0
engine type = openmm
engine class = openmm
engine openmm_simulation = simulation
engine openmm_module = openmm_sim.py
engine subcycles = 10
engine exe_path = /home/enrico/Research/pyretis/examples/test/test-pyvisa/test-openmm-pyvisa
engine rgen = rgen
tis freq = 0.5
tis maxlength = 20000
tis aimless = True
tis allowmaxlength = False
tis zero_momentum = False
tis rescale_energy = False
tis sigma_v = -1
tis seed = 0
tis high_accept = False
tis rgen = rgen
tis shooting_move = sh
tis shooting_moves = []
tis mirror_freq = 0
tis target_freq = 0
tis target_indices = []
tis ensemble_number = 0
tis detect = 0.35
retis swapfreq = 0.5
retis relative_shoots = None
retis nullmoves = True
retis swapsimul = True
orderparameter class = Distance
orderparameter index = (0, 1)
orderparameter periodic = False
orderparameter name = Order Parameter
initial-path method = kick
initial-path kick-from = previous
output backup = overwrite
output pathensemble-file = 1
output screen = 10
output order-file = 1
output energy-file = 1
output trajectory-file = 1
output cross-file = 1
output restart-file = 1
particles type = internal
analysis blockskip = 1
analysis bins = 100
analysis maxblock = 1000
analysis maxordermsd = -1
analysis ngrid = 1001
analysis plot plotter = mpl
analysis plot output = png
analysis plot style = pyretis
analysis report = ['latex', 'rst', 'html']
analysis skipcross = 1000
analysis txt-output = txt.gz
analysis tau_ref_bin = []
analysis skip = 0
interface = 0.3

Ensemble
--------
heading = {'text': 'Retis 2-Ar LJ example (OpenMM)\n=============================='}
simulation task = retis
simulation steps = 5
simulation interfaces = [0.3, 0.3, 0.4]
simulation exe_path = /home/enrico/Research/pyretis/examples/test/test-pyvisa/test-openmm-pyvisa
simulation rgen = rgen
simulation priority_shooting = False
simulation flux = True
simulation zero_ensemble = True
system units = gromacs
system dimensions = 3
system temperature = 1.0
engine type = openmm
engine class = openmm
engine openmm_simulation = simulation
engine openmm_module = openmm_sim.py
engine subcycles = 10
engine exe_path = /home/enrico/Research/pyretis/examples/test/test-pyvisa/test-openmm-pyvisa
engine rgen = rgen
tis freq = 0.5
tis maxlength = 20000
tis aimless = True
tis allowmaxlength = False
tis zero_momentum = False
tis rescale_energy = False
tis sigma_v = -1
tis seed = 0
tis high_accept = False
tis rgen = rgen
tis shooting_move = sh
tis shooting_moves = []
tis mirror_freq = 0
tis target_freq = 0
tis target_indices = []
tis ensemble_number = 1
tis detect = 0.35
retis swapfreq = 0.5
retis relative_shoots = None
retis nullmoves = True
retis swapsimul = True
orderparameter class = Distance
orderparameter index = (0, 1)
orderparameter periodic = False
orderparameter name = Order Parameter
initial-path method = kick
initial-path kick-from = previous
output backup = overwrite
output pathensemble-file = 1
output screen = 10
output order-file = 1
output energy-file = 1
output trajectory-file = 1
output cross-file = 1
output restart-file = 1
particles type = internal
analysis blockskip = 1
analysis bins = 100
analysis maxblock = 1000
analysis maxordermsd = -1
analysis ngrid = 1001
analysis plot plotter = mpl
analysis plot output = png
analysis plot style = pyretis
analysis report = ['latex', 'rst', 'html']
analysis skipcross = 1000
analysis txt-output = txt.gz
analysis tau_ref_bin = []
analysis skip = 0
interface = 0.3

Ensemble
--------
heading = {'text': 'Retis 2-Ar LJ example (OpenMM)\n=============================='}
simulation task = retis
simulation steps = 5
simulation interfaces = [0.3, 0.35, 0.4]
simulation exe_path = /home/enrico/Research/pyretis/examples/test/test-pyvisa/test-openmm-pyvisa
simulation rgen = rgen
simulation priority_shooting = False
simulation flux = True
simulation zero_ensemble = True
system units = gromacs
system dimensions = 3
system temperature = 1.0
engine type = openmm
engine class = openmm
engine openmm_simulation = simulation
engine openmm_module = openmm_sim.py
engine subcycles = 10
engine exe_path = /home/enrico/Research/pyretis/examples/test/test-pyvisa/test-openmm-pyvisa
engine rgen = rgen
tis freq = 0.5
tis maxlength = 20000
tis aimless = True
tis allowmaxlength = False
tis zero_momentum = False
tis rescale_energy = False
tis sigma_v = -1
tis seed = 0
tis high_accept = False
tis rgen = rgen
tis shooting_move = sh
tis shooting_moves = []
tis mirror_freq = 0
tis target_freq = 0
tis target_indices = []
tis ensemble_number = 2
tis detect = 0.4
retis swapfreq = 0.5
retis relative_shoots = None
retis nullmoves = True
retis swapsimul = True
orderparameter class = Distance
orderparameter index = (0, 1)
orderparameter periodic = False
orderparameter name = Order Parameter
initial-path method = kick
initial-path kick-from = previous
output backup = overwrite
output pathensemble-file = 1
output screen = 10
output order-file = 1
output energy-file = 1
output trajectory-file = 1
output cross-file = 1
output restart-file = 1
particles type = internal
analysis blockskip = 1
analysis bins = 100
analysis maxblock = 1000
analysis maxordermsd = -1
analysis ngrid = 1001
analysis plot plotter = mpl
analysis plot output = png
analysis plot style = pyretis
analysis report = ['latex', 'rst', 'html']
analysis skipcross = 1000
analysis txt-output = txt.gz
analysis tau_ref_bin = []
analysis skip = 0
interface = 0.35

Analysis settings
-----------------
blockskip = 1
bins = 100
maxblock = 1000
maxordermsd = -1
ngrid = 1001
plot = {'output': 'png', 'plotter': 'mpl', 'style': 'pyretis'}
report = ['latex', 'rst', 'html']
skipcross = 1000
txt-output = 'txt.gz'
tau_ref_bin = []
skip = 0