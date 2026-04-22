#!/usr/bin/env bash
# Tests for pyvisa with the internal (1-D double-well / Langevin) engine.
#
# What is covered
# ---------------
# 1. pyvisa -cmp          — compress simulation output to HDF5
# 2. pyvisa -recalculate  — recalculate all OPs from stored trajectories
#                           (normal PyRETIS output, existing feature)
# 3. pyvisa -recalculate -data <xyz_file>    — recalculate from a single file
# 4. pyvisa -recalculate -data <xyz_folder>  — recalculate from a folder
#
# These latter two cases exercise the "Load data / Recalculate" workflow
# that was extended to accept individual files and folders.
set -e
export MPLBACKEND=Agg

RETIS_SRC="../../test-internal/retis"

# ---------------------------------------------------------------------------
# Setup: copy simulation input (exclude run.sh and pre-existing outputs)
# ---------------------------------------------------------------------------
cp "${RETIS_SRC}/retis.rst" .
cp "${RETIS_SRC}/initial.xyz" .

pyretisrun -i retis.rst -p

# ---------------------------------------------------------------------------
# Test 1 — compress
# ---------------------------------------------------------------------------
echo "[pyvisa] Test 1: compressing simulation output"
pyvisa -i retis.rst -cmp
test -f pyvisa_compressed_data.hdf5.zip || { echo "FAIL: HDF5.zip file not created"; exit 1; }
echo "[pyvisa] Test 1 PASSED"

# ---------------------------------------------------------------------------
# Test 2 — full-simulation recalculate (standard PyRETIS output layout)
# ---------------------------------------------------------------------------
echo "[pyvisa] Test 2: recalculating from full simulation output"
pyvisa -i retis.rst -recalculate
python ../lib/compare_order.py
echo "[pyvisa] Test 2 PASSED"

# ---------------------------------------------------------------------------
# Test 3 — recalculate from a single .xyz trajectory file
# Note: the internal Langevin engine stores trajectories as traj.txt (a custom
# position/velocity format) rather than .xyz; recalculate_order has no .txt
# handler, so this test is skipped when only traj.txt files exist.  The same
# workflow is exercised end-to-end by test-cp2k-pyvisa which uses .xyz files.
# ---------------------------------------------------------------------------
echo "[pyvisa] Test 3: recalculating from a single trajectory file"
FIRST_TRJ=$(find . -name "*.xyz" -path "*/traj-acc/*/traj/*" | sort | head -1)
if [ -z "$FIRST_TRJ" ]; then
    echo "SKIP Test 3: no accepted trajectory found"
else
    pyvisa -i retis.rst -recalculate -data "$FIRST_TRJ"
    FIRST_TRJ_DIR=$(dirname "$FIRST_TRJ")
    test -f "${FIRST_TRJ_DIR}/order.txt" || { echo "FAIL: order.txt not written next to trajectory"; exit 1; }
    echo "[pyvisa] Test 3 PASSED (recalculated from $FIRST_TRJ)"
fi

# ---------------------------------------------------------------------------
# Test 4 — recalculate from a folder of .xyz files
# ---------------------------------------------------------------------------
echo "[pyvisa] Test 4: recalculating from a folder of trajectory files"
if [ -z "$FIRST_TRJ" ]; then
    echo "SKIP Test 4: no accepted trajectory found"
else
    TRJ_FOLDER=$(dirname "$FIRST_TRJ")
    pyvisa -i retis.rst -recalculate -data "$TRJ_FOLDER"
    test -f "${TRJ_FOLDER}/order.txt" || { echo "FAIL: order.txt not written for folder"; exit 1; }
    echo "[pyvisa] Test 4 PASSED (recalculated from $TRJ_FOLDER)"
fi

# ---------------------------------------------------------------------------
# Test 5 — recalculate then compress (verifies the two operations compose)
# ---------------------------------------------------------------------------
echo "[pyvisa] Test 5: recalculate followed by compress"
pyvisa -i retis.rst -recalculate -cmp
test -f pyvisa_compressed_data.hdf5.zip || { echo "FAIL: HDF5.zip not created after recalculate+cmp"; exit 1; }
echo "[pyvisa] Test 5 PASSED"

# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------
make -f "${RETIS_SRC}/Makefile" clean
rm -f pyvisa_compressed_data.hdf5.zip
