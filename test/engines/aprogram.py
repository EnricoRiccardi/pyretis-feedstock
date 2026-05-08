#! /usr/bin/env python3
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This module will pretend to be an external program."""
import sys


if __name__ == '__main__':
    # pylint: disable=invalid-name
    print('This is a program for testing external commands.', file=sys.stdout)
    args = sys.argv
    if len(args) > 1:
        print('ERROR: Program got arguments:', file=sys.stderr, end='\n')
        argstring = ' '.join(args)
        print(argstring, file=sys.stderr, end='\n')
        sys.exit(1)
    else:
        sys.exit(0)
