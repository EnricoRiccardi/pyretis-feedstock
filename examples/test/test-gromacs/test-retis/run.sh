#!/usr/bin/env bash
set -e
# Disable HWLOC hardware detection components that may hang on some systems
# when X11 display sockets are in a broken state (e.g. full accept queue).
export HWLOC_COMPONENTS=-gl,x11,opencl,cuda
make clean
gmx=${1:-gmx_d}
echo "Using gmx=$gmx"
replace="s/GMXCOMMAND/$gmx/g"

cd gromacs1
sed -e $replace retis.rst > retis-run.rst
cp ../../gmx/gromacs.py .
cp ../../gmx/orderp.py .
pyretisrun -i retis-run.rst -p -l DEBUG
rm retis-run.rst
rm gromacs.py
rm orderp.py
cd ..

cd gromacs2
sed -e $replace retis.rst > retis-run.rst
cp ../../gmx/gromacs.py .
cp ../../gmx/orderp.py .
pyretisrun -i retis-run.rst -p -l DEBUG
rm retis-run.rst
rm gromacs.py
rm orderp.py
cd ..

gmxversion=$($gmx --version | grep -i "gromacs version")
python compare.py "$gmxversion"
make clean
