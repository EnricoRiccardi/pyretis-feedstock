# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Simple script to compare the outcome of two pandas dataframe files."""
import sys
import pandas as pd


def main():
    """Compare two pandas files.

    Returns
    -------
    out : int
        0 if the comparison was successful, 1 otherwise.

    """
    df1 = pd.read_hdf('pyvisa_compressed_data.hdf5', key='data')
    df2 = pd.read_hdf('results/pyvisa_compressed_data.hdf5', key='data')

    # Compare pandas DataFrames _data_

    # Remove NaN, since NaN != NaN
    for ens1, ens2 in zip(df1, df2):
        for trj1, trj2 in zip(ens1, ens2):
            ens1[trj1].frames.notna()
            ens2[trj2].frames.notna()

            # Now compare the remaining
            if not ens1[trj1].frames.equals(ens1[trj2].frames):
                print(ens1[trj1].frames)
                print('The files contain different data')
                return 1

    # Compare pandas DataFrames _settings_
    df1 = pd.read_hdf('pyvisa_compressed_data.hdf5', key='infos')
    df2 = pd.read_hdf('results/pyvisa_compressed_data.hdf5', key='infos')
    if not df1.equals(df2):
        print('The files contain different settings')
        return 1

    return 0


if __name__ == '__main__':
    sys.exit(main())
