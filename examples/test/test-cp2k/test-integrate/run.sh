#!/usr/bin/env bash
set -e
# Disable HWLOC hardware detection components that may hang on some systems
# when X11 display sockets are in a broken state (e.g. full accept queue).
export HWLOC_COMPONENTS=-gl,x11,opencl,cuda
make clean
pyretisrun -i cp2k.rst -p
python compare_energies.py energy.txt pyretis-cp2k-1.ener
mkdir -p cp2k-run
python make_inp.py cp2k.rst ../../shared_input/cp2k/cp2k-run.inp cp2k-run/cp2k.inp
cp ../../shared_input/cp2k/h2.xyz cp2k-run/
cd cp2k-run
cp2k -i cp2k.inp
cd ..
python compare_cp2k_energies.py pyretis-cp2k-1.ener cp2k-run/TEST-1.ener
make clean
