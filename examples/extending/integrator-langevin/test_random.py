# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""The Langevin integrator implemented in FORTRAN."""
import logging
import numpy as np
from matplotlib import pyplot as plt
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())
try:
    from vvintegrator import vvintegrator
except ImportError:
    MSG = ('Could not import external FORTRAN library.'
           '\nPlease compile with "make"!')
    logger.critical(MSG)
    raise ImportError(MSG)


plt.style.use('seaborn-colorblind')
SEED = 1000


def test_rangaussian(sigma, numbers=10000, bins=100):
    """Test the gaussian generator."""
    rndf = np.array([vvintegrator.rangaussian(sigma) for _ in range(numbers)])
    rnd = np.random.normal(loc=0.0, scale=sigma, size=numbers)

    mini = min(rndf.min(), rnd.min())
    maxi = max(rndf.max(), rnd.max())

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    histf, edgesf = np.histogram(rndf, bins=bins, range=(mini, maxi))
    midf = 0.5 * (edgesf[1:] + edgesf[:-1])
    hist, edges = np.histogram(rnd, bins=bins, range=(mini, maxi))
    mid = 0.5 * (edges[1:] + edges[:-1])
    ax1.plot(midf, histf, 'o-', label='FORTRAN', lw=2, alpha=0.7)
    ax1.plot(mid, hist, 's-', label='numpy', lw=2, alpha=0.7)
    ax1.set_title(f'Numbers = {numbers}, bins = {bins}')
    ax1.legend()
    plt.show()


def test_gssbivar(s12os11, sqrts11, sqrtsos11, numbers=10000, bins=100):
    """Test the gaussian bivariate distribution."""
    rndf = []
    for _ in range(numbers):
        rndf.append(vvintegrator.gssbivar(s12os11, sqrts11, sqrtsos11))
    rndf = np.array(rndf)
    histf, xedges, yedges = np.histogram2d(rndf[:, 0], rndf[:, 1],
                                           bins=bins, range=None)
    xvarf, yvarf = np.meshgrid(xedges, yedges)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.pcolormesh(xvarf, yvarf, histf)
    plt.show()


def main():
    """Run the test."""
    size = vvintegrator.get_seed_size()
    seeds = np.array([SEED + i for i in range(size)], dtype='i4')
    vvintegrator.seed_random_generator(seeds)
    for sig in [1.0, 0.5, 2.0, 10.0]:
        test_rangaussian(sig, numbers=10**6)
    test_gssbivar(1.0, 2.0, 1.0)


if __name__ == '__main__':
    main()
