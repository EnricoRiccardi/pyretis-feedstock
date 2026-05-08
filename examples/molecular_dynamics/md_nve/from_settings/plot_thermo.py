# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Plot results using matplotlib."""
import numpy as np
from matplotlib import pyplot as plt


def make_plot():
    """Create the plot and show it."""
    plt.style.use('seaborn-poster')
    raw_data = np.loadtxt('thermo.txt')
    data = {}
    for i, key in enumerate(('step', 'temp', 'pot', 'kin', 'tot', 'press')):
        data[key] = raw_data[:, i]

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.grid()
    ax2.grid()
    ax1.plot(data['step'], data['pot'], lw=3, alpha=0.8,
             label='Potential')
    ax1.plot(data['step'], data['kin'], lw=3,
             alpha=0.8, label='Kinetic')
    ax1.plot(data['step'], data['tot'], lw=3,
             alpha=0.8, label='Total')
    ax1.set_xlabel('Step')
    ax1.set_ylabel('Energy per particle')
    ax1.legend()
    ax1.set_xlim(data['step'][0], data['step'][-1])
    linep, = ax2.plot(data['step'], data['press'], lw=3,
                      alpha=0.8, label='Pressure')
    ax3 = ax2.twinx()

    color_cycle = ax2._get_lines.prop_cycler
    color = next(color_cycle)['color']
    linet, = ax3.plot(data['step'], data['temp'], lw=3,
                      alpha=0.8, label='Temperature', color=color)
    ax2.set_ylabel('Pressure')
    ax3.set_ylabel('Temperature')
    ax2.legend([linep, linet], ['Pressure', 'Temperature'], loc='center')
    ax2.set_xlim(data['step'][0], data['step'][-1])
    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    make_plot()
