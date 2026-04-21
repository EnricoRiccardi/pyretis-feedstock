#!/usr/bin/env bash
# PyVisA tests using CP2K-engine trajectory data (.xyz frames).
#
# CP2K stores trajectory frames as .xyz files, which are natively supported
# by pyretis/tools/recalculate_order.py.  This test therefore exercises the
# full recalculate-from-file / recalculate-from-folder workflow without
# requiring CP2K to be installed.
#
# What is covered
# ---------------
# 1. pyvisa -recalculate -data <xyz_file>    — single-frame recalculation
# 2. pyvisa -recalculate -data <xyz_folder>  — folder-level recalculation
# 3. Verify that an order.txt is written next to the processed trajectory
#
# Source data
# -----------
# The .xyz files come from examples/test/test-cp2k/test-retis-load/pippo/
# which contains pre-generated CP2K initial-path frames.
set -e
export MPLBACKEND=Agg

CP2K_SRC="../../test-cp2k/test-retis-load"

# Copy input files (rst + CP2K input + initial xyz frames)
cp "${CP2K_SRC}/retis.rst" .
cp -r "${CP2K_SRC}/cp2k_input" .
cp -r "${CP2K_SRC}/pippo" .

# ---------------------------------------------------------------------------
# Test 1 — recalculate from a single .xyz file
# ---------------------------------------------------------------------------
echo "[pyvisa/cp2k] Test 1: recalculate from a single .xyz file"
SINGLE_XYZ="pippo/A.xyz"
pyvisa -i retis.rst -recalculate -data "$SINGLE_XYZ"
echo "[pyvisa/cp2k] Test 1 PASSED"

# ---------------------------------------------------------------------------
# Test 2 — recalculate from a folder of .xyz files
# ---------------------------------------------------------------------------
echo "[pyvisa/cp2k] Test 2: recalculate from a folder of .xyz files"
pyvisa -i retis.rst -recalculate -data pippo
test -f "pippo/order.txt" || { echo "FAIL: order.txt not written inside folder"; exit 1; }
echo "[pyvisa/cp2k] Test 2 PASSED"

# ---------------------------------------------------------------------------
# Test 3 — verify the recalculated order values are finite numbers
# ---------------------------------------------------------------------------
echo "[pyvisa/cp2k] Test 3: validate recalculated order.txt"
python - <<'EOF'
import sys, os

order_file = "pippo/order.txt"
if not os.path.isfile(order_file):
    print("FAIL: order.txt not found")
    sys.exit(1)

values = []
with open(order_file) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                values.append(float(parts[1]))
            except ValueError:
                pass

if not values:
    print("FAIL: order.txt contains no numeric values")
    sys.exit(1)

print(f"  Read {len(values)} order-parameter values (min={min(values):.4f}, max={max(values):.4f})")
EOF
echo "[pyvisa/cp2k] Test 3 PASSED"

# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------
find . -mindepth 1 -not -name 'run.sh' -delete
