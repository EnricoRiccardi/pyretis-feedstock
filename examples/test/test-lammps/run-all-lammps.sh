#!/usr/bin/env bash
set -e

export MPLBACKEND=Agg
basedir=$(pwd)

declare -a full_tests=(
                       "test-dump-phasepoint"
                       "test-integrate"
                       "test-modify-velocities"
                       "test-propagate"
                       "test-retis"
                      )

declare -a smoke_tests=(
                        "test-integrate"
                        "test-propagate"
                        "test-retis"
                       )

case "${PYRETIS_LAMMPS_SCOPE:-full}" in
    full)
        tests=("${full_tests[@]}")
        ;;
    smoke)
        tests=("${smoke_tests[@]}")
        ;;
    *)
        echo "Unknown PYRETIS_LAMMPS_SCOPE=${PYRETIS_LAMMPS_SCOPE}" >&2
        echo "Expected 'full' or 'smoke'." >&2
        exit 2
        ;;
esac

echo ""
echo ""
echo ""
echo ""
echo ""
echo "============================================================"
echo "  LAMMPS TEST SUITE  (scope: ${PYRETIS_LAMMPS_SCOPE:-full})"
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
echo "  LAMMPS TEST SUITE -- DONE"
echo "  Finished:    $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo ""
