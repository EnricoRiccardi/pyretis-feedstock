#!/usr/bin/env bash
set -e
export MPLBACKEND=Agg
basedir=$(pwd)

cd examples/test/test-internal/
./run-all-internal.sh
cd "$basedir"

cd examples/test/test-gromacs/
./run-all-gromacs.sh
cd "$basedir"

cd examples/test/test-cp2k/
./run-all-cp2k.sh
cd "$basedir"

cd examples/test/test-lammps/
./run-all-lammps.sh
cd "$basedir"

cd examples/test/test-pyvisa/
./run-all-pyvisa.sh
cd "$basedir"
