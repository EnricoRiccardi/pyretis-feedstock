#!/usr/bin/env bash
set -e
make clean
pyretisrun -i retis.rst -p
pyretisanalyse -i retis.rst
python compare.py
make clean
