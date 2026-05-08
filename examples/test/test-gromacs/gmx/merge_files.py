# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Merge the output from two PyRETIS simulations.

Here we merge the output of two consecutive simulations (i.e. one
is the continuation of the other). Note that we make some specific
assumtions on how the output is written and that we do not attempt
to merge trajectory archive files.

"""
import os
import sys
from pyretis.inout.common import make_dirs
from pyretis.inout.settings import parse_settings_file
from pyretis.core.pathensemble import generate_ensemble_name


def merge_files(file1, file2, outputfile):
    """Merge two files to a single file.

    This method is used to merge PyRETIS output files into a
    single file. Here we assume that the data for the first
    cycle in the second file is identical to the data for the
    last cycle in the first file, and we write this cycle only
    once to the merged output file.

    Parameters
    ----------
    file1 : string
        Path to the first file.
    file2 : string
        Path to the second file.
    outputfile : string
        The path to the merged file we will create here.

    """
    with open(outputfile, 'w') as output:
        # Write full contents of file1:
        with open(file1, 'r') as infile1:
            output.write(infile1.read())
        # Skip data for the first cycle of file2:
        skip = 1
        with open(file2, 'r') as infile2:
            for lines in infile2:
                # The pathensemble.txt files are different from the rest,
                # and we only need to skip the first line.
                if os.path.basename(file2) == 'pathensemble.txt':
                    skip -= 1
                elif lines.startswith('# Cycle'):
                    # For the other files we need to skip to the next
                    # cycle output:
                    skip -= 1
                if skip < 0:
                    output.write(lines)


def main(indir1, indir2, outdir):
    """Merge files in the two given directories."""
    settings = parse_settings_file(os.path.join(indir1, 'retis.rst'))
    inter = settings['simulation']['interfaces']
    for i, _ in enumerate(inter):
        ensemble_dir = generate_ensemble_name(i)
        output_ensemble_dir = os.path.join(outdir, ensemble_dir)
        make_dirs(output_ensemble_dir)
        for filei in {'energy.txt', 'order.txt', 'pathensemble.txt'}:
            file1 = os.path.join(indir1, ensemble_dir, filei)
            file2 = os.path.join(indir2, ensemble_dir, filei)
            merged = os.path.join(output_ensemble_dir, filei)
            merge_files(file1, file2, merged)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1], sys.argv[2], sys.argv[3]))
