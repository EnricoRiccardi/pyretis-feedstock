#!/usr/bin/env bash
set -e
basedir=$(pwd)
# Disable HWLOC hardware detection components that may hang on some systems
# when X11 display sockets are in a broken state (e.g. full accept queue).
export HWLOC_COMPONENTS=-gl,x11,opencl,cuda
make clean
gmx=${1:-gmx_d}
echo "Using gmx=$gmx"
replace="s/GMXCOMMAND/$gmx/g"

gmxversion=$($gmx --version | grep -i "gromacs version")
echo "$gmxversion"

sed -e $replace retis-load-rc.rst > retis-load-rc-run.rst
cp ../../../gmx/orderp.py .
pyretisrun -i retis-load-rc-run.rst -p
rm retis-load-rc-run.rst
rm orderp.py
make clean
