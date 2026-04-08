#!/usr/bin/env bash

set -e
# Disable HWLOC hardware detection components that may hang on some systems
# when X11 display sockets are in a broken state (e.g. full accept queue).
export HWLOC_COMPONENTS=-gl,x11,opencl,cuda
make clean
gmx=${1:-gmx_d}
echo "Using gmx=$gmx"
replace="s/GMXCOMMAND/$gmx/g"
sed -e "$replace" gromacs.rst > gromacs-run.rst

pyretisrun -i gromacs-run.rst # -p

$gmx energy -f pyretis-gmx.edr -o gmx-energy.xvg <<EOF
Potential
Kinetic-En.
Total-Energy
Temperature
EOF

python ../compare_energies.py energy.txt gmx-energy.xvg plot

cp ../../gmx/make_mdp.py .
python make_mdp.py gromacs-run.rst mdout.mdp gromacs-run.mdp
rm make_mdp.py
$gmx grompp -f gromacs-run.mdp -p ../../../shared_input/gromacs/topol.top -c ../../../shared_input/gromacs/conf.g96 -o gromacs-run.tpr
$gmx mdrun -s gromacs-run.tpr -deffnm gromacs-run

$gmx energy -f gromacs-run.edr -o gmx-gromacs-run.xvg <<EOF
Potential
Kinetic-En.
Total-Energy
Temperature
EOF

python ../compare_energies.py energy.txt gmx-gromacs-run.xvg plot
make clean
