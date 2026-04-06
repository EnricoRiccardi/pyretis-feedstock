#!/usr/bin/env bash
set -e
# Disable HWLOC hardware detection components that may hang on some systems
# when X11 display sockets are in a broken state (e.g. full accept queue).
export HWLOC_COMPONENTS=-gl,x11,opencl,cuda
make clean
gmx=${1:-gmx_d}
echo "Using gmx=$gmx"
replace="s/GMXCOMMAND/$gmx/g"
sed -e "$replace" engine.rst > engine-run.rst
python test_gromacs.py 1 plot
python test_gromacs.py 2 plot
rm engine-run.rst
make clean
