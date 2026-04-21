RETIS 1D SS-WT-WF test
======================

Simulation settings
-------------------
task = 'retis'
steps = 15
interfaces = [-0.99, -0.7, -0.6, -0.5]
exe_path = '/home/enrico/Research/pyretis/examples/test/test-internal/retis-ss-wt-wf/ss-wf-wt-sh'
rgen = 'rgen'
priority_shooting = False
flux = True
zero_ensemble = True

System settings
---------------
units = 'reduced'
dimensions = 1
temperature = 0.7

Engine settings
---------------
class = 'Langevin'
timestep = 0.05
gamma = 0.2
high_friction = False
seed = 0
exe_path = '/home/enrico/Research/pyretis/examples/test/test-internal/retis-ss-wt-wf/ss-wf-wt-sh'
rgen = 'rgen'
type = 'internal'

Box settings
------------
periodic = [False]

Particles settings
------------------
position = {'input_file': 'initial.xyz'}
mass = {'Ar': 1.0}
name = ['Ar']
ptype = [0]
type = 'internal'

Forcefield settings
-------------------
description = '1D double well'

Potential
---------
class = 'DoubleWell'
a = 1.0
b = 2.0
c = 0.0

Orderparameter settings
-----------------------
class = 'Position'
dim = 'x'
index = 0
periodic = False
name = 'Order Parameter'

Output settings
---------------
backup = 'overwrite'
order-file = 1
trajectory-file = -1
energy-file = 1
cross-file = 1
pathensemble-file = 1
restart-file = 1
screen = 10

TIS settings
------------
freq = 0.0
maxlength = 10000
aimless = True
allowmaxlength = False
zero_momentum = False
rescale_energy = False
sigma_v = -1
seed = 0
shooting_moves = ['ss', 'wf', 'wt', 'sh']
n_jumps = 6
interface_sour = -0.8
high_accept = False
rgen = 'rgen'
shooting_move = 'sh'
mirror_freq = 0
target_freq = 0
target_indices = []

Initial-path settings
---------------------
method = 'kick'
kick-from = 'initial'

RETIS settings
--------------
swapfreq = 0.5
relative_shoots = None
nullmoves = True
swapsimul = True

Ensemble
--------
heading = {'text': 'Retis 1D example\n================'}
simulation task = retis
simulation steps = 15
simulation interfaces = [-inf, -0.99, -0.99]
simulation exe_path = /home/enrico/Research/pyretis/examples/test/test-internal/retis-ss-wt-wf/ss-wf-wt-sh
simulation rgen = rgen
simulation priority_shooting = False
simulation flux = True
simulation zero_ensemble = True
system units = reduced
system dimensions = 1
system temperature = 0.7
box periodic = [False]
engine class = Langevin
engine timestep = 0.05
engine gamma = 0.2
engine high_friction = False
engine seed = 0
engine exe_path = /home/enrico/Research/pyretis/examples/test/test-internal/retis-ss-wt-wf/ss-wf-wt-sh
engine rgen = rgen
engine type = internal
tis freq = 0.0
tis maxlength = 10000
tis aimless = True
tis allowmaxlength = False
tis zero_momentum = False
tis rescale_energy = False
tis sigma_v = -1
tis seed = 0
tis shooting_moves = ['ss', 'wf', 'wt', 'sh']
tis n_jumps = 6
tis interface_sour = -0.8
tis high_accept = False
tis rgen = rgen
tis shooting_move = ss
tis mirror_freq = 0
tis target_freq = 0
tis target_indices = []
tis ensemble_number = 0
tis detect = -0.7
retis swapfreq = 0.5
retis relative_shoots = None
retis nullmoves = True
retis swapsimul = True
initial-path method = kick
initial-path kick-from = initial
particles position input_file = initial.xyz
particles mass Ar = 1.0
particles name = ['Ar']
particles ptype = [0]
particles type = internal
forcefield description = 1D double well
potential0 class = DoubleWell
potential0 a = 1.0
potential0 b = 2.0
potential0 c = 0.0
orderparameter class = Position
orderparameter dim = x
orderparameter index = 0
orderparameter periodic = False
orderparameter name = Order Parameter
output backup = overwrite
output order-file = 1
output trajectory-file = -1
output energy-file = 1
output cross-file = 1
output pathensemble-file = 1
output restart-file = 1
output screen = 10
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
interface = -0.99

Ensemble
--------
heading = {'text': 'Retis 1D example\n================'}
simulation task = retis
simulation steps = 15
simulation interfaces = [-0.99, -0.99, -0.5]
simulation exe_path = /home/enrico/Research/pyretis/examples/test/test-internal/retis-ss-wt-wf/ss-wf-wt-sh
simulation rgen = rgen
simulation priority_shooting = False
simulation flux = True
simulation zero_ensemble = True
system units = reduced
system dimensions = 1
system temperature = 0.7
box periodic = [False]
engine class = Langevin
engine timestep = 0.05
engine gamma = 0.2
engine high_friction = False
engine seed = 0
engine exe_path = /home/enrico/Research/pyretis/examples/test/test-internal/retis-ss-wt-wf/ss-wf-wt-sh
engine rgen = rgen
engine type = internal
tis freq = 0.0
tis maxlength = 10000
tis aimless = True
tis allowmaxlength = False
tis zero_momentum = False
tis rescale_energy = False
tis sigma_v = -1
tis seed = 0
tis shooting_moves = ['ss', 'wf', 'wt', 'sh']
tis n_jumps = 6
tis interface_sour = -0.8
tis high_accept = False
tis rgen = rgen
tis shooting_move = wf
tis mirror_freq = 0
tis target_freq = 0
tis target_indices = []
tis ensemble_number = 1
tis detect = -0.7
retis swapfreq = 0.5
retis relative_shoots = None
retis nullmoves = True
retis swapsimul = True
initial-path method = kick
initial-path kick-from = initial
particles position input_file = initial.xyz
particles mass Ar = 1.0
particles name = ['Ar']
particles ptype = [0]
particles type = internal
forcefield description = 1D double well
potential0 class = DoubleWell
potential0 a = 1.0
potential0 b = 2.0
potential0 c = 0.0
orderparameter class = Position
orderparameter dim = x
orderparameter index = 0
orderparameter periodic = False
orderparameter name = Order Parameter
output backup = overwrite
output order-file = 1
output trajectory-file = -1
output energy-file = 1
output cross-file = 1
output pathensemble-file = 1
output restart-file = 1
output screen = 10
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
interface = -0.99

Ensemble
--------
heading = {'text': 'Retis 1D example\n================'}
simulation task = retis
simulation steps = 15
simulation interfaces = [-0.99, -0.7, -0.5]
simulation exe_path = /home/enrico/Research/pyretis/examples/test/test-internal/retis-ss-wt-wf/ss-wf-wt-sh
simulation rgen = rgen
simulation priority_shooting = False
simulation flux = True
simulation zero_ensemble = True
system units = reduced
system dimensions = 1
system temperature = 0.7
box periodic = [False]
engine class = Langevin
engine timestep = 0.05
engine gamma = 0.2
engine high_friction = False
engine seed = 0
engine exe_path = /home/enrico/Research/pyretis/examples/test/test-internal/retis-ss-wt-wf/ss-wf-wt-sh
engine rgen = rgen
engine type = internal
tis freq = 0.0
tis maxlength = 10000
tis aimless = True
tis allowmaxlength = False
tis zero_momentum = False
tis rescale_energy = False
tis sigma_v = -1
tis seed = 0
tis shooting_moves = ['ss', 'wf', 'wt', 'sh']
tis n_jumps = 6
tis interface_sour = -0.8
tis high_accept = False
tis rgen = rgen
tis shooting_move = wt
tis mirror_freq = 0
tis target_freq = 0
tis target_indices = []
tis ensemble_number = 2
tis detect = -0.6
retis swapfreq = 0.5
retis relative_shoots = None
retis nullmoves = True
retis swapsimul = True
initial-path method = kick
initial-path kick-from = initial
particles position input_file = initial.xyz
particles mass Ar = 1.0
particles name = ['Ar']
particles ptype = [0]
particles type = internal
forcefield description = 1D double well
potential0 class = DoubleWell
potential0 a = 1.0
potential0 b = 2.0
potential0 c = 0.0
orderparameter class = Position
orderparameter dim = x
orderparameter index = 0
orderparameter periodic = False
orderparameter name = Order Parameter
output backup = overwrite
output order-file = 1
output trajectory-file = -1
output energy-file = 1
output cross-file = 1
output pathensemble-file = 1
output restart-file = 1
output screen = 10
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
interface = -0.7

Ensemble
--------
heading = {'text': 'Retis 1D example\n================'}
simulation task = retis
simulation steps = 15
simulation interfaces = [-0.99, -0.6, -0.5]
simulation exe_path = /home/enrico/Research/pyretis/examples/test/test-internal/retis-ss-wt-wf/ss-wf-wt-sh
simulation rgen = rgen
simulation priority_shooting = False
simulation flux = True
simulation zero_ensemble = True
system units = reduced
system dimensions = 1
system temperature = 0.7
box periodic = [False]
engine class = Langevin
engine timestep = 0.05
engine gamma = 0.2
engine high_friction = False
engine seed = 0
engine exe_path = /home/enrico/Research/pyretis/examples/test/test-internal/retis-ss-wt-wf/ss-wf-wt-sh
engine rgen = rgen
engine type = internal
tis freq = 0.0
tis maxlength = 10000
tis aimless = True
tis allowmaxlength = False
tis zero_momentum = False
tis rescale_energy = False
tis sigma_v = -1
tis seed = 0
tis shooting_moves = ['ss', 'wf', 'wt', 'sh']
tis n_jumps = 6
tis interface_sour = -0.8
tis high_accept = False
tis rgen = rgen
tis shooting_move = sh
tis mirror_freq = 0
tis target_freq = 0
tis target_indices = []
tis ensemble_number = 3
tis detect = -0.5
retis swapfreq = 0.5
retis relative_shoots = None
retis nullmoves = True
retis swapsimul = True
initial-path method = kick
initial-path kick-from = initial
particles position input_file = initial.xyz
particles mass Ar = 1.0
particles name = ['Ar']
particles ptype = [0]
particles type = internal
forcefield description = 1D double well
potential0 class = DoubleWell
potential0 a = 1.0
potential0 b = 2.0
potential0 c = 0.0
orderparameter class = Position
orderparameter dim = x
orderparameter index = 0
orderparameter periodic = False
orderparameter name = Order Parameter
output backup = overwrite
output order-file = 1
output trajectory-file = -1
output energy-file = 1
output cross-file = 1
output pathensemble-file = 1
output restart-file = 1
output screen = 10
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
interface = -0.6

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