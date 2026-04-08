# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Update copyright info for PyRETIS."""
import pathlib
import re
import subprocess
import sys
import tempfile
from os import environ
from datetime import date

NEW_YEAR = str(date.today().year)
COPYRIGHT = re.compile('copyright', re.IGNORECASE)
# Assume that we won't do PyRETIS after 2099 ;-)
YEAR = re.compile('20[0-9][0-9]')
SUFFIX_TO_IGNORE = {'pyc', 'coverage', 'swp'}
# Set up environment for subprocess:
ENV = {}
for key in {'SYSTEMROOT', 'PATH'}:
    val = environ.get(key)
    if val is not None:
        ENV[key] = val
ENV['LANGUAGE'] = 'C'
ENV['LANG'] = 'C'
ENV['LC_ALL'] = 'C'


def is_file_in_repo(filename):
    """Test if file is part of the git repo."""
    is_in_repo = False
    try:
        out = subprocess.Popen(
            ['git', 'ls-files', '--error-unmatch', filename],
            stdin=subprocess.PIPE,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            env=ENV,
            shell=False,
        )
        out.communicate()
        is_in_repo = out.returncode == 0
    except OSError:
        is_in_repo = False
    return is_in_repo


def look_for_files(rootdir):
    """Locate files recursively in the given root directory."""
    filepaths = [x for x in rootdir.rglob('*') if x.is_file()]
    for filepath in filepaths:
        ignore = (filepath.suffix in SUFFIX_TO_IGNORE)
        if not ignore:
            is_in_repo = is_file_in_repo(filepath)
            if is_in_repo:
                yield filepath


def replace_copyright(filename):
    """Replace copyright information in the given input file."""
    replace = False
    print('Processing file: {}'.format(filename))
    with tempfile.NamedTemporaryFile('r+') as tempf:
        with filename.open('r') as inputfile:
            try:
                for lines in inputfile:
                    copy = [i for i in COPYRIGHT.findall(lines)]
                    year = [i for i in YEAR.findall(lines)]
                    if copy and year:
                        print('\t- Found copyright & year')
                        replace = True
                        # Replace the years found with the new year:
                        for i in year:
                            lines = lines.replace(i, NEW_YEAR)
                    tempf.write(lines)
            except UnicodeDecodeError:
                print('\t- Skipping file due to decoding error...')
        tempf.flush()
        if replace:
            # Move to start of file again
            tempf.seek(0)
            with filename.open('w') as outputfile:
                outputfile.write(tempf.read())
    return replace


def main(rootdir):
    """Search and replace copyright information in all files found."""
    replaced = []
    rootdir = pathlib.Path(rootdir)
    for filename in look_for_files(rootdir):
        if replace_copyright(filename):
            replaced.append(filename)
    if replaced:
        print('The following files were updated:')
        for i in replaced:
            print('\t{}'.format(i))
    else:
        print('No files were updated.')


if __name__ == '__main__':
    main(sys.argv[1])
