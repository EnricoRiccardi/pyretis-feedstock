#!/usr/bin/env bash
set -e

basedir=$(pwd)
declare -a tests=(
                  "test-gromacs-pyvisa"
                  "test-gromacs-retis-pyvisa"
                  "test-internal-pyvisa"
                  "test-cp2k-pyvisa"
                  "test-lammps-pyvisa"
                  "test-openmm-pyvisa"
                  )

echo ""
echo ""
echo ""
echo ""
echo ""
echo ""
echo "============================================================"
echo "  PYVISA TEST SUITE"
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
echo "  PYVISA TEST SUITE -- DONE"
echo "  Finished:    $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo ""
