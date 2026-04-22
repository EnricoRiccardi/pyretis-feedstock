#!/usr/bin/env bash
set -e
export MPLBACKEND=Agg
make -C .. clean
pyretisrun -i lammps.rst
make -C .. clean
