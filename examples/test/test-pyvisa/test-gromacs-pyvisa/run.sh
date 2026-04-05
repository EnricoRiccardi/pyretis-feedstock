#!/usr/bin/env bash
set -e
# Disable HWLOC hardware detection components that may hang on some systems
# when X11 display sockets are in a broken state (e.g. full accept queue).
export HWLOC_COMPONENTS=-gl,x11,opencl,cuda
export MPLBACKEND=Agg
cp -nr ../../test-gromacs/test-load/test-load-sparse/load-traj/* . 

pyvisa -i retis-load-rc.rst -recalculate -data ../test-gromacs-pyvisa 
pyvisa -i retis-load-rc.rst -recalculate -data pippo 

find . -mindepth 1 -not -name 'run.sh' -delete
