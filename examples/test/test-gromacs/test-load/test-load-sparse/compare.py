# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Shared sparse-load validator for the GROMACS load tests."""
import argparse
import importlib.util
from pathlib import Path
import sys
import colorama


def _load_shared_compare():
    """Load the shared GROMACS comparison helper."""
    module_path = Path(__file__).resolve().parents[2] / 'gmx' / 'compare.py'
    spec = importlib.util.spec_from_file_location(
        'pyretis_gmx_compare', module_path
    )
    module = importlib.util.module_from_spec(spec)
    if spec.loader is None:
        raise ImportError(f'Could not load comparison helper: {module_path}')
    spec.loader.exec_module(module)
    return module


def main():
    """Run sparse-load checks for the current test directory."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--settings',
        default='retis-load-rc-run.rst',
        help='Settings file to validate.',
    )
    parser.add_argument(
        '--valid-ends',
        nargs='*',
        default=('R', 'L'),
        help='Accepted path end states.',
    )
    args = parser.parse_args()

    compare = _load_shared_compare()
    return compare.main_sparse_load(
        settings_file=args.settings,
        valid_ends=tuple(args.valid_ends),
    )


if __name__ == '__main__':
    colorama.init(autoreset=True)
    sys.exit(main())
