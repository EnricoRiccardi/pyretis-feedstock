#!/usr/bin/env bash
# Run the OpenMM RETIS test and compare results to reference.
#
# IMPORTANT: This test requires OpenMM (openmm package) to be installed.
# Install with:  conda install -c conda-forge openmm
set -e
export HWLOC_COMPONENTS=-gl,x11,opencl,cuda
export MPLBACKEND=Agg

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

make -C "$SCRIPT_DIR" clean

cd "$SCRIPT_DIR/openmm"
pyretisrun -i retis.rst -p -l DEBUG
cd "$SCRIPT_DIR"

python compare.py

make -C "$SCRIPT_DIR" clean
