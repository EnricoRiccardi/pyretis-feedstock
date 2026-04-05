#!/usr/bin/env bash
shopt -s extglob
set -e
# Disable HWLOC hardware detection components that may hang on some systems
# when X11 display sockets are in a broken state (e.g. full accept queue).
export HWLOC_COMPONENTS=-gl,x11,opencl,cuda
make clean
gmx=${1:-gmx_d}
echo "Using gmx=$gmx"
replace="s/GMXCOMMAND/$gmx/g"

cd run-gromacs-ini
cp ../../gmx/gromacs.py .
cp ../../gmx/orderp.py .
cp ../../../shared_input/gromacs gromacs_input_ini -r
sed -e $replace retis.rst > retis-run.rst
pyretisrun -i retis-run.rst -p
cd ..

# Test source dir change
mv run-gromacs-ini/gromacs_input_ini run-gromacs1/gromacs_input1

cp -r run-gromacs-ini/!(*.rst) run-gromacs1
cp -r run-gromacs-ini/!(*.rst) run-gromacs2


cd run-gromacs1
sed -e $replace retis.rst > retis-run.rst
pyretisrun -i retis-run.rst -p
cd ..
cd run-gromacs2
sed -e $replace retis.rst > retis-run.rst
pyretisrun -i retis-run.rst -p
cd ..

cp ../gmx/compare.py .
python compare.py run-gromacs1 run-gromacs2 --energy_skip 'vpot'
rm compare.py
rm -- */gromacs.py
rm -- */orderp.py
make clean
