#!/usr/bin/env bash
make clean
gmx=${1:-gmx}
echo "Using gmx=$gmx"
replace="s/GMXCOMMAND/$gmx/g"
sed -e "$replace" engine.rst > engine-run.rst
python test_gromacs.py 2
