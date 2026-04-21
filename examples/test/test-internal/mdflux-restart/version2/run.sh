#!/usr/bin/env bash
set -e
cp ../Makefile .
make clean
cp ../initial.xyz .
cd run-full
# todo -p by using the progress bar, the test fails.
# a smaller number of cycles get printed out.
pyretisrun -i flux.rst # -p
cd ..
cd run-step1
pyretisrun -i flux.rst # -p
cd ..
cd run-step2 
pyretisrun -i flux.rst # -p
cd ..
cp ../compare.py .
python compare.py
rm initial.xyz
rm compare.py
make clean
rm Makefile
