#!/usr/bin/env bash
set -e

basedir=$(pwd)
# Disable HWLOC hardware detection components that may hang on some systems
# when X11 display sockets are in a broken state (e.g. full accept queue).
export HWLOC_COMPONENTS=-gl,x11,opencl,cuda
export MPLBACKEND=Agg
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/pyretis-matplotlib}"
mkdir -p "$MPLCONFIGDIR"

declare -a full_tests=(
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

declare -a smoke_tests=(
                        "test-gromacs"
                        "test-integrate/gromacs"
                        "test-integrate/gromacs2"
                        "test-load/test-load-sparse/load-traj"
                        "test-load/test-load-sparse/load-frames"
                        "test-restart/test-continue"
                        )

case "${PYRETIS_GROMACS_SCOPE:-full}" in
    full)
        tests=("${full_tests[@]}")
        ;;
    smoke)
        tests=("${smoke_tests[@]}")
        ;;
    *)
        echo "Unknown PYRETIS_GROMACS_SCOPE=${PYRETIS_GROMACS_SCOPE}" >&2
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
echo "  GROMACS TEST SUITE  (scope: ${PYRETIS_GROMACS_SCOPE:-full})"
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
echo "  GROMACS TEST SUITE -- DONE"
echo "  Finished:    $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo ""
