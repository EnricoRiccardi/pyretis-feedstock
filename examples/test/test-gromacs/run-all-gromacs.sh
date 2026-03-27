#!/usr/bin/env bash
set -e

basedir=$(pwd)
# Disable HWLOC hardware detection components that may hang on some systems
# when X11 display sockets are in a broken state (e.g. full accept queue).
export HWLOC_COMPONENTS=-gl,x11,opencl,cuda

declare -a tests=(
                  "test-gromacs"
                  "test-gromacs-gromacs2"
                  "test-integrate/gromacs"
                  "test-integrate/gromacs2"
		  "test-load/test-initialise"
                  "test-load/test-load"
                  "test-load/test-load-sparse/load-traj"
                  "test-load/test-load-sparse/load-frames"
                  "test-repptis"
                  "test-restart/test-continue"
                  "test-restart/test-initialise"
                  "test-restart/test-restart"
                  "test-retis"
                  )

for i in "${tests[@]}"
do
    echo "Running in: $i"
    cd "$i"
    ./run.sh
    cd "$basedir"
done
