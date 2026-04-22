#!/usr/bin/env bash
# Run all OpenMM tests.
#
# IMPORTANT: This script requires OpenMM to be installed.
# Install with:  conda install -c conda-forge openmm
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== test-retis ==="
cd "$SCRIPT_DIR/test-retis"
bash run.sh
