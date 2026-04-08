#!/usr/bin/env bash
set -e
make clean
pyretisrun -i repptis.rst -p
pyretisanalyse -i repptis.rst
python compare.py
make clean
