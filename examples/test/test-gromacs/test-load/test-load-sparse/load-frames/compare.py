# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Run sparse-load frame checks via the shared GROMACS compare script."""
import importlib.util
import logging
import sys
from pathlib import Path


def _load_shared_compare():
    """Load the shared GROMACS compare module from the sibling gmx folder."""
    compare_path = Path(__file__).resolve().parents[3] / 'gmx' / 'compare.py'
    spec = importlib.util.spec_from_file_location(
        'pyretis_gmx_compare', compare_path
    )
    module = importlib.util.module_from_spec(spec)
    if spec.loader is None:
        raise ImportError(f'Could not load compare helper from {compare_path}')
    spec.loader.exec_module(module)
    return module


def main():
    """Run sparse-load checks using the shared GROMACS compare helper."""
    return _load_shared_compare().main_sparse_load(
        settings_file='retis-load-rc-run.rst',
        valid_ends=('R', 'L'),
    )


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,
                        format='%(levelname)s: %(message)s')
    sys.exit(main())
