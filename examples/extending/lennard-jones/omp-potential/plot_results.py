"""Plot results from time tests of the Lennard-Jones potential.

Here we assume that the results are available in files named
``timeings.txt``.
"""
import numpy as np
from matplotlib import pylab as plt


data = {'Serial': np.loadtxt('timings-plain.txt'),
        'OpenMP, 2 threads': np.loadtxt('timings-OMP_NUM_THREADS_2.txt')}

markers = ['o', 's', '^', 'v']

fig = plt.figure(figsize=(12, 5))
ax1 = fig.add_subplot(121)
for i, key in enumerate(data):
    result = data[key]
    npart, avg, std = result[:, 0], result[:, 2], result[:, 3]
    ax1.errorbar(npart, avg, yerr=std, label=key, lw=3, markersize=9,
                 marker=markers[i], alpha=0.9)
ax1.legend(loc='best')
ax1.set_xlabel('System size')
ax1.set_ylabel('Time')

ax2 = fig.add_subplot(122)
norm_avg = data['OpenMP, 2 threads'][:, 2]

for i, key in enumerate(data):
    result = data[key]
    npart, avg = result[:, 0], result[:, 2]
    nmin = min(len(norm_avg), len(avg))
    avg_n = avg[:nmin] / norm_avg[:nmin]
    ax2.plot(npart, avg_n, lw=3, label=key, markersize=9, alpha=0.9,
             marker=markers[i])
ax2.legend(loc='best', bbox_to_anchor=(0.3, 0.7))
ax2.set_xlabel('System size')
ax2.set_ylabel('Time relative to OpenMP with 2 threads')
fig.tight_layout()
plt.show()
