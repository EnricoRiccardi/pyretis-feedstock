#!/usr/bin/env bash
set -e

# Disable HWLOC hardware detection components that may hang on some systems
# when X11 display sockets are in a broken state (e.g. full accept queue).
export HWLOC_COMPONENTS=-gl,x11,opencl,cuda

make clean
gmx=${1:-gmx_d}
echo "Using gmx=$gmx"
replace="s/GMXCOMMAND/$gmx/g"

sed -e $replace repptis.rst > repptis-run.rst
cp ../gmx/gromacs.py .
cp ../gmx/orderp.py .
pyretisrun -i repptis-run.rst -l DEBUG
python compare.py
rm repptis-run.rst
rm gromacs.py
rm orderp.py

make clean
