#!/usr/bin/env bash
set -e

basedir=$(pwd)
declare -a tests=(
                  "test-gromacs-pyvisa"
		  "test-gromacs-retis-pyvisa"
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

# Source of the gromacs load-traj fixture that the pyvisa tests reuse.
LOADTRAJ="$basedir/../test-gromacs/test-load/test-load-sparse/load-traj"

for i in "${tests[@]}"
do
    echo "------------------------------------------------------------"
    echo "  [START] $i"
    echo "          $(date '+%Y-%m-%d %H:%M:%S')"
    echo "------------------------------------------------------------"
    start=$(date +%s)
    cd "$i"
    # test-gromacs-pyvisa runs the upstream load-traj run.sh verbatim,
    # so its prerequisites (Makefile, retis-load-rc.rst, gromacs_input/,
    # orderp.py, pippo/, …) need to be staged in the test directory and
    # ../compare.py needs to point at the shared sparse-load wrapper.
    # test-gromacs-retis-pyvisa stages everything itself inside its
    # own run.sh.
    if [ "$i" = "test-gromacs-pyvisa" ]; then
        find . -mindepth 1 -not -name 'run.sh' -delete
        cp -r "$LOADTRAJ"/Makefile "$LOADTRAJ"/retis-load-rc.rst \
            "$LOADTRAJ"/orderp.py "$LOADTRAJ"/gromacs_input \
            "$LOADTRAJ"/pippo .
    fi
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
