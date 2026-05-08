# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Methods that might be useful for testing.

This module defines generic methods for testing.

"""
import os


def search_for_files(rootdir, match=None):
    """Find files by walking the given directory.

    Parameters
    ----------
    rootdir : string
        The path where we will search from.
    match : string, optional
        If given, the method will only return files that
        are equal to the given match.

    Return
    ------
    out : list of strings
        The paths of the found files.

    """
    files = []
    for root, _, filei in os.walk(rootdir):
        for i in filei:
            if match is None:
                files.append(os.path.join(root, i))
            elif i == match:
                files.append(os.path.join(root, i))
    return files


def clean_dir(dirname):
    """Remove ALL files in the given directory."""
    for files in os.listdir(dirname):
        filename = os.path.join(dirname, files)
        if os.path.isfile(filename):
            os.remove(filename)
