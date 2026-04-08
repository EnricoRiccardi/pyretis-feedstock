#!/usr/bin/env bash
set -e

basedir=$(pwd)

declare -a tests=("test-cp2k"
                  "test-integrate"
                  "test-retis-load")

for i in "${tests[@]}"
do
    echo "Running in: $i"
    cd "$i"
    ./run.sh
    make clean
    cd "$basedir"
done

