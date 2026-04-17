#!/usr/bin/env bash
set -e

basedir=$(pwd)

declare -a tests=("test-cp2k"
                  "test-integrate"
                  "test-retis-load")

echo ""
echo ""
echo ""
echo ""
echo ""
echo "============================================================"
echo "  CP2K TEST SUITE"
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
echo "  CP2K TEST SUITE -- DONE"
echo "  Finished:    $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo ""

