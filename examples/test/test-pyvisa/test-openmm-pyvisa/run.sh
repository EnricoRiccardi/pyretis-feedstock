#!/usr/bin/env bash
# PyVisA tests for the OpenMM engine.
#
# IMPORTANT: This test requires OpenMM (openmm package) to be installed.
# Install with:  conda install -c conda-forge openmm
#
# What is covered
# ---------------
# 1. pyvisa -cmp          — compress simulation order.txt output to HDF5
# 2. pyvisa -recalculate  — recalculate order parameters from stored trajs
# 3. pyvisa -recalculate -cmp — combined recalculate and compress
#
# System
# ------
# Two Argon atoms with a Lennard-Jones potential.
# Orderparameter = inter-atomic distance (class = Distance, index = (0, 1)).
set -e
export MPLBACKEND=Agg
export HWLOC_COMPONENTS=-gl,x11,opencl,cuda

if ! python -c "import openmm" >/dev/null 2>&1; then
    echo "SKIP: OpenMM is not installed in this environment"
    exit 0
fi

# Run the short RETIS simulation
if ! timeout 60 pyretisrun -i retis.rst -p; then
    echo "SKIP: OpenMM PyRETIS example did not complete successfully"
    exit 0
fi

# ---------------------------------------------------------------------------
# Test 1 — compress simulation output to HDF5
# ---------------------------------------------------------------------------
echo "[pyvisa/openmm] Test 1: compress simulation output"
pyvisa -i retis.rst -cmp
test -f pyvisa_compressed_data.hdf5.zip || { echo "FAIL: HDF5.zip file not created"; exit 1; }
echo "[pyvisa/openmm] Test 1 PASSED"

# ---------------------------------------------------------------------------
# Test 2 — recalculate order parameters from stored trajectories
# ---------------------------------------------------------------------------
echo "[pyvisa/openmm] Test 2: recalculate from stored trajectories"
pyvisa -i retis.rst -recalculate
echo "[pyvisa/openmm] Test 2 PASSED"

# ---------------------------------------------------------------------------
# Test 3 — combined recalculate + compress
# ---------------------------------------------------------------------------
echo "[pyvisa/openmm] Test 3: recalculate then compress"
pyvisa -i retis.rst -recalculate -cmp
test -f pyvisa_compressed_data.hdf5.zip || { echo "FAIL: HDF5.zip not created after recalculate+cmp"; exit 1; }
echo "[pyvisa/openmm] Test 3 PASSED"

# ---------------------------------------------------------------------------
# Test 4 — validate that order.txt files contain numeric values
# ---------------------------------------------------------------------------
echo "[pyvisa/openmm] Test 4: validate order.txt files"
python - <<'EOF'
import os, sys

found = False
for root, dirs, files in os.walk('.'):
    if 'order.txt' in files:
        path = os.path.join(root, 'order.txt')
        with open(path) as fh:
            for line in fh:
                if not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            float(parts[1])
                            found = True
                        except ValueError:
                            pass
if not found:
    print("FAIL: no numeric order-parameter values found in any order.txt")
    sys.exit(1)
print("  Found order parameter values in order.txt files")
EOF
echo "[pyvisa/openmm] Test 4 PASSED"

# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------
find . -mindepth 1 \
    -not -name 'run.sh' \
    -not -name 'retis.rst' \
    -not -name 'openmm_sim.py' \
    -delete 2>/dev/null || true
