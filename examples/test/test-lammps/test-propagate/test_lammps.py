# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""A test using the LAMMPS engine."""
import os
import sys
import colorama
import numpy as np
from matplotlib import pyplot as plt
from pyretis.core.path import Path
from pyretis.inout.common import make_dirs
from pyretis.engines.lammps import LAMMPSEngine
from pyretis.core.random_gen import create_random_generator
from pyretis.testing.helpers import clean_dir
from pyretis.testing.systemhelp import create_system_ext
from plotting import plot_compare, plot_xy
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


HERE = os.path.abspath(os.path.dirname(__file__))
STEPS = 100
SUBCYCLES = 2


def run_forward():
    """Use the LAMMPS engine to propagate a MD simulation forward in time."""
    logger.info('\nTesting LAMMPS engine: propagating forward & backward.')
    engine = LAMMPSEngine('lmp_serial', 'lammps_input', SUBCYCLES)
    # Create a dummy system:
    system = create_system_ext(pos=('system.data', 0))
    exe_dir = os.path.join(HERE, 'run-forward-backward')
    make_dirs(exe_dir)
    clean_dir(exe_dir)
    engine.exe_dir = exe_dir
    pathf = Path(rgen=None, maxlen=STEPS)
    logger.info('-> Propagating forward....')
    ensemble = {
        'system': system,
        'order_function': None,
        'interfaces': [-1, 0.0, 2.0]
    }
    engine.propagate(path=pathf, ensemble=ensemble, reverse=False)
    # Propagate from the last point, but backward:
    logger.info('-> Setting system to last point in the forward path.')
    last = pathf.phasepoints[-1]
    pathb = Path(rgen=None, maxlen=pathf.length)
    logger.info('-> Propagating backward....')
    ensemble = {
        'system': last,
        'order_function': None,
        'interfaces': [-1, 0.0, 3.0]
    }
    engine.propagate(path=pathb, ensemble=ensemble, reverse=True)
    logger.info('-> Comparing forward & backward:')
    order_f = np.array([i.order for i in pathf.phasepoints])
    order_b = np.array([i.order for i in reversed(pathb.phasepoints)])
    data_sets = [
        [
            (range(len(order_f)), order_f[:, 0], 'Forward.'),
            (range(len(order_b)), order_b[:, 0], 'Backward.'),
        ],
        [
            (range(len(order_f)), order_f[:, 1], 'Forward.'),
            (range(len(order_b)), -1*order_b[:, 1], 'Backward.'),
        ]
    ]
    plot_compare(data_sets, ['OP1', 'OP2'])
    data_sets_xy = [
        (order_f[:, 0], order_b[:, 0], 'OP1, forward', 'OP1, backward'),
        (order_f[:, 1], -1*order_b[:, 1], 'OP2, forward', 'OP2, backward'),
    ]
    plot_xy(data_sets_xy)


def run_backward():
    """Use the LAMMPS engine to propagate a MD simulation backward in time."""
    logger.info('\nTesting LAMMPS engine: propagating backward & forward.')
    engine = LAMMPSEngine('lmp_serial', 'lammps_input', SUBCYCLES)
    # Create a dummy system:
    system = create_system_ext(pos=('system.data', 0))
    exe_dir = os.path.join(HERE, 'run-backward-forward')
    make_dirs(exe_dir)
    clean_dir(exe_dir)
    engine.exe_dir = exe_dir
    rgen = create_random_generator({'seed': 2})
    ensemble = {
        'system': system,
        'rgen': rgen,
    }
    vel_settings = {'aimless': True}
    engine.modify_velocities(ensemble, vel_settings)
    pathb = Path(rgen=None, maxlen=STEPS)
    logger.info('-> Propagating backward....')
    ensemble = {
        'system': system,
        'order_function': None,
        'interfaces': [1.2, 5.0, 12.0],
    }
    engine.propagate(path=pathb, ensemble=ensemble, reverse=True)
    # Propagate from the last point, but forward:
    logger.info('-> Setting system to last point in the backward path.')
    last = pathb.phasepoints[-1]
    pathf = Path(rgen=None, maxlen=pathb.length)
    logger.info('-> Propagating forward....')
    ensemble = {
        'system': last,
        'order_function': None,
        'interfaces': [1.0, 5.0, 12.0],
    }
    engine.propagate(path=pathf, ensemble=ensemble, reverse=False)
    logger.info('-> Comparing backward & forward:')
    order_f = np.array([i.order for i in reversed(pathf.phasepoints)])
    order_b = np.array([i.order for i in pathb.phasepoints])
    data_sets = [
        [
            (range(len(order_f)), order_f[:, 0], 'Forward.'),
            (range(len(order_b)), order_b[:, 0], 'Backward.'),
        ],
        [
            (range(len(order_f)), -1*order_f[:, 1], 'Forward.'),
            (range(len(order_b)), order_b[:, 1], 'Backward.'),
        ]
    ]
    plot_compare(data_sets, ['OP1', 'OP2'])
    data_sets_xy = [
        (order_f[:, 0], order_b[:, 0], 'OP1, forward', 'OP1, backward'),
        (-1*order_f[:, 1], order_b[:, 1], 'OP2, forward', 'OP2, backward'),
    ]
    plot_xy(data_sets_xy)


def main(plot=False):
    """Run the comparisons."""
    run_forward()
    run_backward()
    if plot:
        plt.show()


if __name__ == '__main__':
    colorama.init(autoreset=True)
    main(plot=len(sys.argv) >= 2)
