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
    for ensembles in os.listdir('.'):
        if ensembles.isdigit():
            with open(ensembles + '/order.txt') as one:
                with open(ensembles + '/order.txt_000') as two:
                    for line1, line2 in zip(one, two):
                        if line1.split()[1] != line2.split()[1]:
                            if round(float(line1.split()[1]), 2) != \
                                    round(float(line2.split()[1]), 2):
                                print('File differs')
                                return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
