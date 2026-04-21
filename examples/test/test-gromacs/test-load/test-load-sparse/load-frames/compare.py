# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Result comparison for GROMACS RETIS load-frames simulation.

We compare a RETIS simulation of 5 steps — with initial paths loaded
from individual gromacs frame files — against stored reference results.

Reproducibility depends on:
  - Fixed TIS seed in retis-load-rc.rst (seed = 0, rgen = rgen-borg)
  - Initial paths loaded from committed frame files in pippo/
  - The GROMACS integrator is deterministic given the loaded velocities

To generate reference results, run the test once with GROMACS available:
    ./run.sh
    mkdir -p results/000 results/001 results/002 results/003
    for d in 000 001 002 003; do
        cp $d/pathensemble.txt $d/order.txt $d/energy.txt results/$d/
    done
"""
import logging
import os
import sys
import colorama
from pyretis.inout.settings import parse_settings_file
from pyretis.core.pathensemble import generate_ensemble_name
from pyretis.testing.simulation_comparison import (
    compare_path_ensemble_data,
    compare_numerical_data,
)

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

RESULTS = 'results'
INPUT_FILE = 'retis-load-rc.rst'


def _section(msg):
    logger.info('\n%s', msg)
    logger.info('=' * len(msg))


def compare_files(settings):
    """Compare pathensemble.txt, order.txt and energy.txt for every ensemble.

    Parameters
    ----------
    settings : dict
        Parsed simulation settings.

    Returns
    -------
    status : int
        0 if all comparisons pass, positive integer otherwise.
    """
    if not os.path.isdir(RESULTS):
        logger.error(
            "Reference directory '%s' does not exist. "
            "Run the test once with GROMACS to generate reference results.",
            RESULTS,
        )
        return 1

    inter = settings['simulation']['interfaces']
    files = ('pathensemble.txt', 'order.txt', 'energy.txt')
    status = 0
    for i in range(len(inter)):
        ens_dir = generate_ensemble_name(i)
        _section(f'Ensemble {ens_dir}')
        for fname in files:
            ref = os.path.join(RESULTS, ens_dir, fname)
            new = os.path.join(ens_dir, fname)
            logger.info('  Comparing %s ...', fname)
            if fname == 'pathensemble.txt':
                equal, msg = compare_path_ensemble_data(new, ref)
            else:
                equal, msg = compare_numerical_data(new, ref)
            if equal:
                logger.info('    -> OK')
            else:
                logger.error('    -> MISMATCH: %s', msg)
                status += 1
    return status


def main():
    """Run the comparison; return combined exit code."""
    settings = parse_settings_file(INPUT_FILE)

    _section('Comparing output files against reference data')
    status = compare_files(settings)

    if status == 0:
        logger.info('\nAll checks passed.')
    else:
        logger.error('\n%d check(s) failed.', status)
    return status


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,
                        format='%(levelname)s: %(message)s')
    colorama.init(autoreset=True)
    sys.exit(main())
