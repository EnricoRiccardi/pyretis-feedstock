#!/usr/bin/env bash
set -e
python -m pytest -n auto --cov=pyretis test/
pycodestyle
pydocstyle --count pyretis
python devtools/run_linting.py -i pyretis
