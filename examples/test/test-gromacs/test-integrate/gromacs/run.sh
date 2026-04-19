#!/usr/bin/env bash

set -e
# Disable HWLOC hardware detection components that may hang on some systems
# when X11 display sockets are in a broken state (e.g. full accept queue).
export HWLOC_COMPONENTS=-gl,x11,opencl,cuda
make clean
gmx=${1:-gmx_d}
logfile="gromacs-run.log"
rm -f "$logfile"
echo "Using gmx=$gmx"
echo "Writing GROMACS command output to $logfile"
replace="s/GMXCOMMAND/$gmx/g"
sed -e "$replace" gromacs.rst > gromacs-run.rst

pyretisrun -i gromacs-run.rst # -p

$gmx energy -f pyretis-gmx.edr -o gmx-energy.xvg >>"$logfile" 2>&1 <<EOF
Potential
Kinetic-En.
Total-Energy
Temperature
EOF

python ../compare_energies.py energy.txt gmx-energy.xvg plot \
    >>"$logfile" 2>&1

cp ../../gmx/make_mdp.py .
python make_mdp.py gromacs-run.rst mdout.mdp gromacs-run.mdp
rm make_mdp.py
$gmx grompp -f gromacs-run.mdp -p ../../../shared_input/gromacs/topol.top \
    -c ../../../shared_input/gromacs/conf.g96 -o gromacs-run.tpr \
    >>"$logfile" 2>&1
$gmx mdrun -s gromacs-run.tpr -deffnm gromacs-run >>"$logfile" 2>&1

$gmx energy -f gromacs-run.edr -o gmx-gromacs-run.xvg >>"$logfile" 2>&1 <<EOF
Potential
Kinetic-En.
Total-Energy
Temperature
EOF

python ../compare_energies.py energy.txt gmx-gromacs-run.xvg plot \
    >>"$logfile" 2>&1
make clean
