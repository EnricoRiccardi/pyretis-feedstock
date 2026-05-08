#!/usr/bin/env python3
"""Run pycodestyle on the PyRETIS codebase and summarise violations.

Usage
-----
    python check_style.py [paths...]

If no paths are given, scans the default set: pyretis/, test/, examples/test/,
examples/path_sampling/, and docs/_static/.

Exit code is the number of violations found (0 = clean).
"""
import os
import subprocess
import sys
from collections import defaultdict


DEFAULT_PATHS = [
    'pyretis',
    'test',
    'examples/test',
    'examples/path_sampling',
    'docs/_static',
]


def run_pycodestyle(paths):
    """Run pycodestyle on *paths* and return (violations, counts)."""
    cmd = [sys.executable, '-m', 'pycodestyle', '--statistics', '-qq'] + paths
    result = subprocess.run(  # nosec B603
        cmd, capture_output=True, text=True, cwd=os.getcwd(), check=False
    )
    lines = (result.stdout + result.stderr).splitlines()

    violations = []
    counts = defaultdict(int)
    for line in lines:
        line = line.strip()
        if not line:
            continue
        # Statistics lines look like: "5       E501 line too long"
        parts = line.split(None, 1)
        if len(parts) == 2 and parts[0].isdigit():
            counts[parts[1]] += int(parts[0])
        else:
            violations.append(line)

    return violations, dict(counts)


def main():
    """Run the style check script."""
    paths = sys.argv[1:] if len(sys.argv) > 1 else DEFAULT_PATHS
    # Filter to existing paths
    paths = [p for p in paths if os.path.exists(p)]
    if not paths:
        print("No paths found to scan.")
        sys.exit(0)

    print(f"Checking: {', '.join(paths)}")
    violations, counts = run_pycodestyle(paths)

    if violations:
        print(f"\n{'='*60}")
        print(f"  {len(violations)} violation(s) found:")
        print(f"{'='*60}")
        for v in violations:
            print(f"  {v}")
        if counts:
            print("\n  Summary by code:")
            for code, n in sorted(counts.items()):
                print(f"    {n:4d}  {code}")
    else:
        print("  No violations found.")

    sys.exit(len(violations))


if __name__ == '__main__':
    main()
