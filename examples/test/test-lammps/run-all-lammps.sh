#!/usr/bin/env bash
set -e

export MPLBACKEND=Agg
basedir=$(pwd)
# declare -a tofix=(
#                   "test-dump-phasepoint"
#                   "test-integrate"
#                   "test-propagate"
#                   "test-modify-velocities"
#                  )
declare -a tests=(
                  "lammps_testing"
                  )

for i in "${tests[@]}"
do
    echo "Running in: $i"
    cd "$i"
    ./run.sh
    make clean
    cd "$basedir"
done

