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

cd run-5
sed -e $replace retis.rst > retis-run.rst
cp ../../../gmx/gromacs.py .
cp ../../../gmx/orderp.py .
pyretisrun -i retis-run.rst -p
cd ..

cd run-initialise
sed -e $replace retis.rst > retis-run.rst
cp ../../../gmx/gromacs.py .
cp ../../../gmx/orderp.py .
pyretisrun -i retis-run.rst -p
# pyretis_gmx_rnd_state=$(realpath rnd.state)
cd ..

cp ../../gmx/copy_last_path.py .
python copy_last_path.py run-initialise run-load/initial_path
rm copy_last_path.py
cd run-load
sed -e $replace retis.rst > retis-run.rst
cp ../../../gmx/gromacs.py .
cp ../../../gmx/orderp.py .
# cp $pyretis_gmx_rnd_state pyretis_gmx_rnd.state
pyretisrun -i retis-run.rst -p
cd ..

# PyRETIS 3 does not support this part of the test.
# cp ../../gmx/compare.py .
# python compare.py run-5 run-load --path_skip 8
# rm compare.py
rm */retis-run.rst
rm */gromacs.py
rm */orderp.py 
make clean
