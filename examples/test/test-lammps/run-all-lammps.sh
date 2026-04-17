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

echo ""
echo ""
echo ""
echo ""
echo ""
echo "============================================================"
echo "  LAMMPS TEST SUITE"
echo "  Total tests: ${#tests[@]}"
echo "  Started:     $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo ""

for i in "${tests[@]}"
do
    echo "------------------------------------------------------------"
    echo "  [START] $i"
    echo "          $(date '+%Y-%m-%d %H:%M:%S')"
    echo "------------------------------------------------------------"
    start=$(date +%s)
    cd "$i"
    ./run.sh
    make clean
    cd "$basedir"
    end=$(date +%s)
    runtime=$((end-start))
    echo "------------------------------------------------------------"
    echo "  [ END ] $i"
    echo "          $(date '+%Y-%m-%d %H:%M:%S')  (elapsed: ${runtime}s)"
    echo "------------------------------------------------------------"
    echo ""
done

echo "============================================================"
echo "  LAMMPS TEST SUITE -- DONE"
echo "  Finished:    $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo ""

