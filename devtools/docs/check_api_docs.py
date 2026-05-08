# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Check that source modules are included in the API docs.

Typical usage is: python check_api_docs.py ../../docs/api/ ../../pyretis

"""
import pathlib
import sys


def read_api_rst(filepath):
    """Yield modules documented in the given file."""
    with filepath.open('r') as infile:
        for lines in infile:
            if lines.find('automodule') != -1:
                yield lines.strip().split()[-1]


def read_api_rst_files(dirpath):
    """Read all .rst files in the given directory."""
    rst_files = [
        i for i in dirpath.rglob('*rst')
    ]
    modules = set([])
    for filepath in rst_files:
        for mod in read_api_rst(filepath):
            modules.add(mod)
    return modules


def find_source_modules(source_dir):
    """Return source modules in the given directory."""
    source_files = set([])
    filenames = source_dir.rglob('*.py')
    for filename in filenames:
        # Figure out the file placement relative to the source
        filename = filename.relative_to(source_dir)
        # Get the stem filename
        filename_e = filename.stem
        if filename_e == '__init__':
            continue
        # Remove the .py suffix by replacing it with an empty string
        filename = filename.with_suffix("")
        mod = 'pyretis.{}'.format('.'.join(filename.parts))
        source_files.add(mod)
    return source_files


def main(docdir, source_dir):
    """Run the check."""
    docdir = pathlib.Path(docdir)
    source_dir = pathlib.Path(source_dir)
    modules = read_api_rst_files(docdir)
    source = find_source_modules(source_dir)
    missing = [i for i in source if i not in modules]
    if missing:
        print('Missing the following modules in the API documentation:')
        print()
        for i in sorted(missing):
            print(i)
        print()
        return 1
    print('No modules missing.')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1], sys.argv[2]))
