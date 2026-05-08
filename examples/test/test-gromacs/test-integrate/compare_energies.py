# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Compare the energy output for GROMACS and PyRETIS."""
import sys
import colorama
from matplotlib import pyplot as plt
import numpy as np
from pyretis.inout.formats.gromacs import (
    read_xvg_file,
)
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


plt.style.use('seaborn-v0_8-poster')


def _show_or_close():
    """Show figures only when the backend supports interactive display."""
    if 'agg' in plt.get_backend().lower():
        plt.close('all')
        return
    plt.show()


def main(energy_file, xvg_file, plot=False):
    """Perform the test."""
    logger.info(f'Reading energy file: {energy_file}')
    energy = np.loadtxt(energy_file)
    logger.info(f'Reading xvg file: {xvg_file}')
    energy_xvg = read_xvg_file(xvg_file)

    mse_ok = obtain_mses(energy, energy_xvg)

    if plot:
        logger.info('\nPlotting for comparison')
        plot_comparison(energy, energy_xvg, energy_file, xvg_file)

    if not mse_ok:
        logger.error('\nComparison failed!')
        sys.exit(1)


def obtain_mses(energy, energy_xvg, tol=1.0e-5):
    """Obtain some mean squared errors."""
    pairs = ((1, 'potential'), (2, 'kinetic en.'), (3, 'total energy'),
             (4, 'temperature'))
    for pair in pairs:
        if not len(energy[:, pair[0]]) == len(energy_xvg[pair[1]]):
            return False
        mse = ((energy[:, pair[0]] - energy_xvg[pair[1]])**2).mean(axis=0)
        level = 'info'
        tol_ok = True
        if tol:
            tol_ok = abs(mse) < tol
            if not tol_ok:
                level = 'error'
        if level == 'error':
            logger.error(f'MSE {pair[1]}: {mse}')
        else:
            logger.info(f'MSE {pair[1]}: {mse}')
        if not tol_ok:
            return False
    return True


def plot_comparison(energy, energy_xvg, energy_file, xvg_file):
    """Just plot some properties for the paths."""
    fig1 = plt.figure()
    ax11 = fig1.add_subplot(221)
    ax12 = fig1.add_subplot(222)
    ax21 = fig1.add_subplot(223)
    ax22 = fig1.add_subplot(224)
    plabel = f'PyRETIS ({energy_file})'
    glabel = f'gmx ({xvg_file})'
    ax11.plot(energy[:, 1], lw=4, ls='-', marker='o', label=plabel)
    ax11.plot(energy_xvg['potential'], lw=2, ls='--', marker='^',
              label=glabel)
    ax11.legend()
    ax11.set_ylabel('Potential')

    ax12.plot(energy[:, 2], lw=4, ls='-', marker='o')
    ax12.plot(energy_xvg['kinetic en.'], lw=2, ls='--', marker='^')
    ax12.set_ylabel('Kinetic')

    ax21.plot(energy[:, 3], lw=4, ls='-', marker='o')
    ax21.plot(energy_xvg['total energy'], lw=2, ls='--', marker='^')
    ax21.set_ylabel('Total')

    ax22.plot(energy[:, 4], lw=4, ls='-', marker='o')
    ax22.plot(energy_xvg['temperature'], lw=2, ls='--', marker='^')
    ax22.set_ylabel('Temperature')

    fig1.tight_layout()
    _show_or_close()


if __name__ == '__main__':
    colorama.init(autoreset=True)
    PLOT = len(sys.argv) > 3
    main(sys.argv[1], sys.argv[2], plot=PLOT)
