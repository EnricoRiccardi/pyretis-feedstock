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

for i in "${tests[@]}"
do
    start=$(date +%s)
    echo "Running in: $i"
    cd "$i"
    ./run.sh
    cd "$basedir"
    end=$(date +%s)
    runtime=$((end-start))
    echo "$runtime $i" >> time_spent.txt
done
