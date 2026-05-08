RETIS GROMACS sparse load (trajectory)
======================================

Simulation settings
-------------------
task = 'retis'
steps = 5
interfaces = [2.19, 2.22, 2.63]
zero_left = 2.16
rgen = 'rgen-borg'
seed = 0
exe_path = '/home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa'
priority_shooting = False
flux = True
zero_ensemble = True
endcycle = 5

System settings
---------------
units = 'gromacs'
dimensions = 3
temperature = 1.0

Engine settings
---------------
class = 'gromacs2'
gmx = 'gmx_d'
mdrun = 'gmx_d mdrun'
input_path = 'gromacs_input'
timestep = 0.002
subcycles = 5
gmx_format = 'gro'
exe_path = '/home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa'
rgen = 'rgen'
type = 'external'
input_files = {'conf': '/home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa/gromacs_input/conf.gro',
               'index': '/home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa/gromacs_input/index.ndx',
               'input_o': '/home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa/gromacs_input/grompp.mdp',
               'topology': '/home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa/gromacs_input/topol.top'}

Particles settings
------------------
type = 'external'

Orderparameter settings
-----------------------
class = 'Distance'
module = 'orderp.py'
index = [0, 3]

Output settings
---------------
pathensemble-file = 1
screen = 10
order-file = 1
energy-file = 1
trajectory-file = 1
backup = 'append'
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
rgen = 'rgen-borg'
seed = 0
high_accept = False
shooting_move = 'sh'
shooting_moves = []
mirror_freq = 0
target_freq = 0
target_indices = []

Initial-path settings
---------------------
method = 'kick'
load_folder = 'pippo'

RETIS settings
--------------
swapfreq = 0.5
relative_shoots = None
nullmoves = True
swapsimul = True

Ensemble
--------
heading = {'text': 'RETIS GROMACS sparse load (trajectory)\n======================================'}
simulation task = retis
simulation steps = 5
simulation interfaces = [2.16, 2.19, 2.19]
simulation zero_left = 2.16
simulation rgen = rgen-borg
simulation seed = 0
simulation exe_path = /home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa
simulation priority_shooting = False
simulation flux = True
simulation zero_ensemble = True
system units = gromacs
system dimensions = 3
system temperature = 1.0
engine class = gromacs2
engine gmx = gmx_d
engine mdrun = gmx_d mdrun
engine input_path = gromacs_input
engine timestep = 0.002
engine subcycles = 5
engine gmx_format = gro
engine exe_path = /home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa
engine rgen = rgen
engine type = external
engine input_files conf = /home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa/gromacs_input/conf.gro
engine input_files input_o = /home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa/gromacs_input/grompp.mdp
engine input_files topology = /home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa/gromacs_input/topol.top
engine input_files index = /home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa/gromacs_input/index.ndx
tis freq = 0.5
tis maxlength = 20000
tis aimless = True
tis allowmaxlength = False
tis zero_momentum = False
tis rescale_energy = False
tis sigma_v = -1
tis rgen = rgen-borg
tis seed = 0
tis high_accept = False
tis shooting_move = sh
tis shooting_moves = []
tis mirror_freq = 0
tis target_freq = 0
tis target_indices = []
tis ensemble_number = 0
tis detect = 2.22
initial-path method = kick
initial-path load_folder = pippo
retis swapfreq = 0.5
retis relative_shoots = None
retis nullmoves = True
retis swapsimul = True
orderparameter class = Distance
orderparameter module = orderp.py
orderparameter index = [0, 3]
output pathensemble-file = 1
output screen = 10
output order-file = 1
output energy-file = 1
output trajectory-file = 1
output backup = append
output cross-file = 1
output restart-file = 1
particles type = external
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
interface = 2.19

Ensemble
--------
heading = {'text': 'RETIS GROMACS sparse load (trajectory)\n======================================'}
simulation task = retis
simulation steps = 5
simulation interfaces = [2.19, 2.19, 2.63]
simulation zero_left = 2.16
simulation rgen = rgen-borg
simulation seed = 0
simulation exe_path = /home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa
simulation priority_shooting = False
simulation flux = True
simulation zero_ensemble = True
system units = gromacs
system dimensions = 3
system temperature = 1.0
engine class = gromacs2
engine gmx = gmx_d
engine mdrun = gmx_d mdrun
engine input_path = gromacs_input
engine timestep = 0.002
engine subcycles = 5
engine gmx_format = gro
engine exe_path = /home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa
engine rgen = rgen
engine type = external
engine input_files conf = /home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa/gromacs_input/conf.gro
engine input_files input_o = /home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa/gromacs_input/grompp.mdp
engine input_files topology = /home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa/gromacs_input/topol.top
engine input_files index = /home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa/gromacs_input/index.ndx
tis freq = 0.5
tis maxlength = 20000
tis aimless = True
tis allowmaxlength = False
tis zero_momentum = False
tis rescale_energy = False
tis sigma_v = -1
tis rgen = rgen-borg
tis seed = 0
tis high_accept = False
tis shooting_move = sh
tis shooting_moves = []
tis mirror_freq = 0
tis target_freq = 0
tis target_indices = []
tis ensemble_number = 1
tis detect = 2.22
initial-path method = kick
initial-path load_folder = pippo
retis swapfreq = 0.5
retis relative_shoots = None
retis nullmoves = True
retis swapsimul = True
orderparameter class = Distance
orderparameter module = orderp.py
orderparameter index = [0, 3]
output pathensemble-file = 1
output screen = 10
output order-file = 1
output energy-file = 1
output trajectory-file = 1
output backup = append
output cross-file = 1
output restart-file = 1
particles type = external
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
interface = 2.19

Ensemble
--------
heading = {'text': 'RETIS GROMACS sparse load (trajectory)\n======================================'}
simulation task = retis
simulation steps = 5
simulation interfaces = [2.19, 2.22, 2.63]
simulation zero_left = 2.16
simulation rgen = rgen-borg
simulation seed = 0
simulation exe_path = /home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa
simulation priority_shooting = False
simulation flux = True
simulation zero_ensemble = True
system units = gromacs
system dimensions = 3
system temperature = 1.0
engine class = gromacs2
engine gmx = gmx_d
engine mdrun = gmx_d mdrun
engine input_path = gromacs_input
engine timestep = 0.002
engine subcycles = 5
engine gmx_format = gro
engine exe_path = /home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa
engine rgen = rgen
engine type = external
engine input_files conf = /home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa/gromacs_input/conf.gro
engine input_files input_o = /home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa/gromacs_input/grompp.mdp
engine input_files topology = /home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa/gromacs_input/topol.top
engine input_files index = /home/enrico/Research/master/examples/test/test-pyvisa/test-gromacs-retis-pyvisa/gromacs_input/index.ndx
tis freq = 0.5
tis maxlength = 20000
tis aimless = True
tis allowmaxlength = False
tis zero_momentum = False
tis rescale_energy = False
tis sigma_v = -1
tis rgen = rgen-borg
tis seed = 0
tis high_accept = False
tis shooting_move = sh
tis shooting_moves = []
tis mirror_freq = 0
tis target_freq = 0
tis target_indices = []
tis ensemble_number = 2
tis detect = 2.63
initial-path method = kick
initial-path load_folder = pippo
retis swapfreq = 0.5
retis relative_shoots = None
retis nullmoves = True
retis swapsimul = True
orderparameter class = Distance
orderparameter module = orderp.py
orderparameter index = [0, 3]
output pathensemble-file = 1
output screen = 10
output order-file = 1
output energy-file = 1
output trajectory-file = 1
output backup = append
output cross-file = 1
output restart-file = 1
particles type = external
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
interface = 2.22

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