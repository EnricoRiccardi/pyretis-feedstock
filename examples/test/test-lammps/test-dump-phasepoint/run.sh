#!/usr/bin/env bash
set -e
export MPLBACKEND=Agg
make clean
python test_lammps.py
cd lammps
./run.sh
cd ..
make clean
