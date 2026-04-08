#!/usr/bin/env bash
set -e
# Disable HWLOC hardware detection components that may hang on some systems
# when X11 display sockets are in a broken state (e.g. full accept queue).
export HWLOC_COMPONENTS=-gl,x11,opencl,cuda
make clean
pyretisrun -i retis.rst -p
python compare.py
make clean
