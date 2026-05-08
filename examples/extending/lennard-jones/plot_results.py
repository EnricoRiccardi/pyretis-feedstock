# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Plot results from time tests of the Lennard-Jones potential.

Here we assume that the results are available in files named
``timings.txt``.
"""
import numpy as np
from matplotlib import pylab as plt

RESULTS = {'c': 'c/timings.txt',
           'fortran': 'fortran/timings.txt',
           'fortran (pointer)': 'fortran-pointer/timings.txt',
           'numpy': 'timings-numpy.txt',
           'pure python': 'timings-python.txt'}

MARKERS = ('o', 'v', '^', 's', '*', 'p', 'h', 'x', '+')


if __name__ == '__main__':
    # pylint: disable=invalid-name
    data = {}
    for key in RESULTS:
        data[key] = np.loadtxt(RESULTS[key])

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    for i, key in enumerate(data):
        result = data[key]
        npart, avg, std = result[:, 0], result[:, 2], result[:, 3]
        ax1.errorbar(npart, avg, yerr=std, label=key, lw=3, markersize=9,
                     marker=MARKERS[i], alpha=0.9)
    ax1.legend(ncol=2, loc='center right')
    ax1.set_xlabel('System size')
    ax1.set_ylabel('Time')
    fig.tight_layout()

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    norm_avg = data['c'][:, 2]

    for i, key in enumerate(data):
        result = data[key]
        npart, avg = result[:, 0], result[:, 2]
        nmin = min(len(norm_avg), len(avg))
        avg_n = avg[:nmin] / norm_avg[:nmin]
        ax2.plot(npart, avg_n, lw=3, label=key, markersize=9, alpha=0.9,
                 marker=MARKERS[i])
    ax2.legend(ncol=2, loc='best', bbox_to_anchor=(0.3, 0.7))
    ax2.set_yscale('log')
    ax2.set_xlabel('System size')
    ax2.set_ylabel('Time relative to python+c')
    fig2.tight_layout()
    plt.show()
