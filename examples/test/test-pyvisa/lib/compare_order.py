# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Simple script to compare the outcome of two pandas dataframe files."""
import os
import sys


def main():
    """Compare all order.txt and order.txt_000 in each ensemble.

    Returns
    -------
    out : int
        0 if the comparison was successful, 1 otherwise.

    """
    compared = 0
    for ensembles in os.listdir('.'):
        if ensembles.isdigit():
            cur = ensembles + '/order.txt'
            bak = ensembles + '/order.txt_000'
            if not os.path.isfile(cur) or not os.path.isfile(bak):
                continue
            with open(cur) as one, open(bak) as two:
                for line1, line2 in zip(one, two):
                    parts1, parts2 = line1.split(), line2.split()
                    if len(parts1) < 2 or len(parts2) < 2:
                        continue
                    if parts1[1] != parts2[1]:
                        if round(float(parts1[1]), 2) != \
                                round(float(parts2[1]), 2):
                            print('File differs')
                            return 1
            compared += 1
    if compared == 0:
        print('No ensemble pairs to compare (no accepted trajectories?)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
