#!/usr/bin/env bash
set -e
make clean
cd retis-full
pyretisrun -i retis.rst # -p
cd ..
cd retis-6
pyretisrun -i retis.rst # -p
cd ..
cd retis-6-12
cp -r ../retis-6/0* .
pyretisrun -i retis.rst # -p
cd ..
cp ../../test-gromacs/gmx/compare.py .
python compare.py retis-6-12 retis-full --traj_skip
rm compare.py
make clean
