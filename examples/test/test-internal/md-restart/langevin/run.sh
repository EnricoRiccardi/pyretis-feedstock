#!/usr/bin/env bash
set -e
CLEAN_MAKEFILE=../Makefile
make -f "${CLEAN_MAKEFILE}" clean
cd run-full
pyretisrun -i md-full.rst #-p
cd ..
cd run-10
pyretisrun -i md-10.rst #-p
cd ..
cd run-10-100
pyretisrun -i md-10-100.rst # -p
cd ..
cp ../compare.py .
python compare.py
make -f "${CLEAN_MAKEFILE}" clean
rm compare.py
