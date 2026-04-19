#!/usr/bin/env bash
set -e
export MPLBACKEND=Agg
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/pyretis-matplotlib}"
# The complete GROMACS matrix is intentionally available, but dominates
# runtime. Use PYRETIS_GROMACS_SCOPE=full ./test-heavy.sh to run it.
export PYRETIS_GROMACS_SCOPE="${PYRETIS_GROMACS_SCOPE:-smoke}"
mkdir -p "$MPLCONFIGDIR"
basedir=$(pwd)

print_heading() {
    local title="$1"
    printf '\n\n'
    printf '================================================================\n'
    printf '  %s\n' "$title"
    printf '================================================================\n\n'
}

run_suite() {
    local title="$1"
    local suite_dir="$2"
    local runner="$3"

    print_heading "$title"
    (
        cd "${basedir}/${suite_dir}"
        "$runner"
    )
}

run_suite "Internal PyRETIS example tests" \
          "examples/test/test-internal" \
          "./run-all-internal.sh"

run_suite "GROMACS example tests (${PYRETIS_GROMACS_SCOPE} scope)" \
          "examples/test/test-gromacs" \
          "./run-all-gromacs.sh"

run_suite "PyVisA example tests" \
          "examples/test/test-pyvisa" \
          "./run-all-pyvisa.sh"

run_suite "CP2K example tests" \
          "examples/test/test-cp2k" \
          "./run-all-cp2k.sh"

run_suite "LAMMPS example tests" \
          "examples/test/test-lammps" \
          "./run-all-lammps.sh"
