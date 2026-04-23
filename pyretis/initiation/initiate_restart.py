# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This file contains functions used for initiation of paths.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

initiate_restart (:py:func:`.initiate_restart`)
    A method that will get the initial path from the output from
    a previous simulation.
"""
import logging
import os
from pyretis.core.pathensemble import generate_ensemble_name
from pyretis.core.random_gen import create_random_generator
from pyretis.inout.restart import read_restart_file
from pyretis.inout.screen import REFERENCE
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = ['initiate_restart']


def initiate_restart(simulation, settings, cycle):
    """Initialise paths by loading restart data.

    Parameters
    ----------
    simulation : object like :py:class:`.Simulation`
        The simulation we are setting up.
    settings : dictionary
        A dictionary with settings for the initiation.
    cycle : integer
        The simulation cycles we are starting at.

    """
    for idx, ensemble in enumerate(simulation.ensembles):
        path_ensemble = ensemble['path_ensemble']
        name = path_ensemble.ensemble_name
        logger.log(REFERENCE, '\n--- Restarting path ensemble: %s ---', name)
        restart_file = os.path.join(
            generate_ensemble_name(path_ensemble.ensemble_number),
            'ensemble.restart')

        if not os.path.isfile(restart_file):
            raise FileNotFoundError(
                f'Restart file not found: {restart_file}')
        restart_info = read_restart_file(restart_file)
        ensemble['engine'].load_restart_info(restart_info['engine'])
        ensemble['engine'].exe_dir = path_ensemble.directory['generate']
        ensemble['system'].load_restart_info(restart_info['system'])
        if 'rgen' in restart_info:
            ensemble['rgen'] = create_random_generator(restart_info['rgen'])
            logger.info('Restored per-ensemble rgen from %s', restart_file)
        ensemble['order_function'].load_restart_info(
            restart_info.get('order_function', []))

        path_ensemble.load_restart_info(restart_info['path_ensemble'],
                                        cycle)
        # This allows ensemble renumbering
        if settings['simulation']['task'] == 'retis':
            path_ensemble.ensemble_number = idx

        # The Force field is not part of the restart, just explicitly
        # set it:
        for point in path_ensemble.last_path.phasepoints:
            point.forcefield = ensemble['system'].forcefield
        path = path_ensemble.last_path

        yield True, path, path.status, path_ensemble
