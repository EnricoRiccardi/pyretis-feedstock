#!/usr/bin/env bash
set -e
# Disable HWLOC hardware detection components that may hang on some systems
# when X11 display sockets are in a broken state (e.g. full accept queue).
export HWLOC_COMPONENTS=-gl,x11,opencl,cuda
export MPLBACKEND=Agg
# Clean any stale files left over from a previous failed run so the source
# rst is guaranteed fresh.
find . -mindepth 1 -not -name 'run.sh' -delete
cp -r ../../test-gromacs/test-load/test-load-sparse/load-traj/* .

# Generate trajectory data first so pyvisa -recalculate has something to
# work on (matches the pattern used in test-gromacs-retis-pyvisa).
gmx=${1:-gmx_d}
sed -e "s/GMXCOMMAND/$gmx/g" retis-load-rc.rst > retis-load-rc-run.rst
pyretisrun -i retis-load-rc-run.rst

# Recalculate collective variables on the freshly produced ensemble data.
pyvisa -i retis-load-rc-run.rst -recalculate

find . -mindepth 1 -not -name 'run.sh' -delete
