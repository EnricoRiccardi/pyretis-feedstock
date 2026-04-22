#!/usr/bin/env bash
set -e
export MPLBACKEND=Agg

lmp=${LMP:-lmp_serial}
echo "Using lmp=$lmp"
replace="s/lmp_serial/$lmp/g"

if [ "${1:-}" = "generate" ]; then
    # Generate reference results and store them in results.tgz
    echo "Generating reference results..."
    make clean
    cd lammps1
    sed -e "$replace" retis.rst > retis-run.rst
    pyretisrun -i retis-run.rst -p -l DEBUG
    rm retis-run.rst
    tar czf results.tgz results/
    echo "Reference results written to lammps1/results.tgz"
    echo "Commit this file to version control."
    cd ..
    make clean
else
    # Run simulation and compare against reference results
    make clean
    cd lammps1
    sed -e "$replace" retis.rst > retis-run.rst
    pyretisrun -i retis-run.rst -p -l DEBUG
    rm retis-run.rst
    cd ..
    python compare.py
    make clean
fi
