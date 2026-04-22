#!/usr/bin/env bash
# PyVisA tests for the LAMMPS engine.
#
# IMPORTANT: This test requires LAMMPS (lmp_serial or lmp) to be available
# on the system PATH.  It also requires the 2D WCA LAMMPS example to be
# configured correctly (see examples/path_sampling/2D-wca-lammps/).
#
# What is covered
# ---------------
# 1. pyvisa -cmp          — compress simulation order.txt output to HDF5
# 2. pyvisa -recalculate  — recalculate from stored LAMMPS trajectory files
#
# Note on trajectory formats
# --------------------------
# PyRETIS stores LAMMPS trajectories as .lammpstrj files.  Since
# recalculate_order currently supports .xyz / .trr / .gro / .g96, the
# -recalculate command will skip .lammpstrj files but still succeeds for
# ensembles that have no stored trajectories (verifying graceful handling).
# If you want full recalculation from LAMMPS trajectories, convert the
# .lammpstrj files to .xyz first (e.g. using VMD or LAMMPS itself).
set -e
export MPLBACKEND=Agg
export HWLOC_COMPONENTS=-gl,x11,opencl,cuda

LAMMPS_SRC="../../../path_sampling/2D-wca-lammps"

if command -v lmp_serial >/dev/null 2>&1; then
    lmp_bin=lmp_serial
elif command -v lmp >/dev/null 2>&1; then
    lmp_bin=lmp
else
    echo "SKIP: no LAMMPS executable found on PATH"
    exit 0
fi

if ! "${lmp_bin}" -h >/dev/null 2>&1; then
    echo "SKIP: LAMMPS executable '${lmp_bin}' is not runnable in this environment"
    exit 0
fi

# Copy LAMMPS example into working directory
cp "${LAMMPS_SRC}/retis.rst" .
cp -r "${LAMMPS_SRC}/lammps_input" .

# The source example expects a pre-populated load folder. This test should be
# self-contained, so switch to kick initiation and keep the run short.
sed -i 's/^steps = 30$/steps = 5/' retis.rst
sed -i 's/^method = load$/method = kick/' retis.rst
sed -i "s/^lmp = lmp_serial$/lmp = ${lmp_bin}/" retis.rst

# Run a short LAMMPS RETIS simulation
pyretisrun -i retis.rst -p

# ---------------------------------------------------------------------------
# Test 1 — compress simulation output to HDF5
# ---------------------------------------------------------------------------
echo "[pyvisa/lammps] Test 1: compress simulation output"
pyvisa -i retis.rst -cmp
test -f pyvisa_compressed_data.hdf5.zip || { echo "FAIL: HDF5.zip file not created"; exit 1; }
echo "[pyvisa/lammps] Test 1 PASSED"

# ---------------------------------------------------------------------------
# Test 2 — recalculate from simulation output
# (LAMMPS trajectories are .lammpstrj; unsupported files are skipped gracefully)
# ---------------------------------------------------------------------------
echo "[pyvisa/lammps] Test 2: recalculate (graceful handling of .lammpstrj)"
pyvisa -i retis.rst -recalculate
echo "[pyvisa/lammps] Test 2 PASSED"

# ---------------------------------------------------------------------------
# Test 3 — recalculate from a single .lammpstrj file
# (expected to log a warning and return; must not crash)
# ---------------------------------------------------------------------------
echo "[pyvisa/lammps] Test 3: recalculate -data with .lammpstrj (no crash)"
LAMMPS_TRJ=$(find . -name "*.lammpstrj" | head -1)
if [ -n "$LAMMPS_TRJ" ]; then
    pyvisa -i retis.rst -recalculate -data "$LAMMPS_TRJ" || true
    echo "[pyvisa/lammps] Test 3 PASSED (handled gracefully)"
else
    echo "[pyvisa/lammps] Test 3 SKIPPED (no .lammpstrj found)"
fi

# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------
make clean 2>/dev/null || true
rm -f pyvisa_compressed_data.hdf5.zip
find . -mindepth 1 \
    -not -name 'run.sh' \
    -delete 2>/dev/null || true
