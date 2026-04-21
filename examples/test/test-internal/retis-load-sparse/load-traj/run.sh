#!/usr/bin/env bash
set -e
make clean
pyretisrun -i retis.rst -p
python compare.py
make clean
