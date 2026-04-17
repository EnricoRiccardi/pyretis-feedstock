#!/usr/bin/env bash
set -e

basedir=$(pwd)

declare -a tests=(
                  "compare-internal-with-lammps"
                  "explorer"
                  "partial-path"
                  "mdflux-restart/version1"
                  "mdflux-restart/version2"
                  "md-restart/langevin"
                  "md-restart/velocity-verlet"
                  "retis"
                  "retis-load-sparse/load-traj"
                  "retis-load-sparse/load-frames"
                  "retis-restart"
                  "retis-ss-wt-wf"
                  "tis-multiple"
                  "tis-restart"
              )

echo ""
echo ""
echo ""
echo ""
echo ""
echo "============================================================"
echo "  INTERNAL ENGINE TEST SUITE"
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
    echo "$runtime $i" >> time_spent.txt
done

echo "============================================================"
echo "  INTERNAL ENGINE TEST SUITE -- DONE"
echo "  Finished:    $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo ""
