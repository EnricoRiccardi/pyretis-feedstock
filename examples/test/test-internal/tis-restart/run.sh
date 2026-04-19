#!/usr/bin/env bash
set -e
make clean
cd run-full
pyretisrun -i tis-001.rst # -p
pyretisanalyse -i tis-001.rst # -p
cd ..
cd run-6
pyretisrun -i tis-001.rst # -p
cd ..
cd run-6-12
cp -r ../run-6/001 .
pyretisrun -i tis-001.rst # -p
cd ..
python compare.py 
make clean
