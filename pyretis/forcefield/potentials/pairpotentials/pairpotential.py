# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This module defines some helper functions for pair potentials.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mixing_parameters (:py:func:`.mixing_parameters`)
    Definition of mixing rules which are used in some force fields
    for generating cross interactions.

generate_pair_interactions (:py:func:`generate_pair_interactions`)
    Function to generate pair parameters from atom parameters.
"""
import itertools
import logging
import numpy as np
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = ['mixing_parameters', 'generate_pair_interactions']


def _check_pair_parameters(parameters):
    """Check that the required pair parameters are given.

    If the parameters are not given, we will just set them to default
    values.

    Parameters
    ----------
    parameters : dict
        The settings for the potential.

    Note
    ----
    This function will not return anything, but it will update the input
    `parameters` with default values.

    """
    for pair in parameters:
        for key in ('epsilon', 'sigma'):
            if key not in parameters[pair]:
                msg = f'{key} for {pair} not given. Set to 0.0'
                logger.warning(msg)
                parameters[pair][key] = 0.0
        if 'rcut' not in parameters[pair]:
            factor = parameters[pair].get('factor', 0.0)
            rcut = factor * parameters[pair]['sigma']
            parameters[pair]['rcut'] = rcut
            msg = ('rcut for {} not given. Using factor to set it to: '
                   '{} * {} = {}')
            msg = msg.format(pair, parameters[pair]['sigma'], factor, rcut)
            logger.info(msg)


def generate_pair_interactions(parameters, mixing):
    """Generate pair parameters from atom parameters.

    The parameters are given as a dictionary where the keys are
    either just integers -- which defines atom parameters -- or tuples
    which define pair interactions.

    Parameters
    ----------
    parameters : dict
        This dict contain the atom parameters.
    mixing : string
        Determines how we should mix pair interactions.

    """
    _check_pair_parameters(parameters)
    atoms = []
    pair_param = {}
    for key in parameters:
        if key == 'mixing' or isinstance(key, tuple):
            continue
        atoms.append(key)
    for (atmi, atmj) in itertools.product(atoms, atoms):
        pari = parameters[atmi]
        parj = parameters[atmj]
        if (atmi, atmj) in pair_param:
            continue
        if (atmi, atmj) in parameters:
            pair_param[atmi, atmj] = dict(parameters[atmi, atmj])
            pair_param[atmj, atmi] = pair_param[atmi, atmj]
            continue
        if (atmj, atmi) in parameters:
            pair_param[atmj, atmi] = dict(parameters[atmj, atmi])
            pair_param[atmi, atmj] = pair_param[atmj, atmi]
            continue
        if atmi == atmj:
            eps_ij = pari['epsilon']
            sig_ij = pari['sigma']
            rcut_ij = pari['rcut']
            pair_param[atmi, atmi] = {'epsilon': eps_ij, 'sigma': sig_ij,
                                      'rcut': rcut_ij}
        else:
            eps_ij, sig_ij, rcut_ij = mixing_parameters([pari['epsilon'],
                                                         parj['epsilon']],
                                                        [pari['sigma'],
                                                         parj['sigma']],
                                                        [pari['rcut'],
                                                         parj['rcut']],
                                                        mixing)
            pair_param[atmi, atmj] = {'epsilon': eps_ij, 'sigma': sig_ij,
                                      'rcut': rcut_ij}
            pair_param[atmj, atmi] = pair_param[atmi, atmj]
    #     if (atmj, atmi) not in pair_param:
    #         pair_param[atmj, atmi] = pair_param[atmi, atmj]
    return pair_param


def mixing_parameters(epsilon, sigma, rcut, mixing='geometric'):
    r"""Define the so-called mixing rules.

    These mixing rules are used for some force fields when generating
    cross interactions. The known mixing rules are:

    1. Geometric:

       * .. math::

            \epsilon_{ij} = \sqrt{\epsilon_{i} \times \epsilon_{j}}

       * .. math::

            \sigma_{ij} = \sqrt{\sigma_{i} \times \sigma_{j}}

       * .. math::

            r_{\text{c},ij} = \sqrt{r_{\text{c},i} \times r_{\text{c},j}}

    2. Arithmetic:

       * .. math::

            \epsilon_{ij} = \sqrt{\epsilon_{i} \times \epsilon_{j}}

       * .. math::

            \sigma_{ij} = \frac{\sigma_{i} \times \sigma_{j}}{2}

       * .. math::

            r_{\text{c},ij} = \frac{r_{\text{c},i} \times r_{\text{c},j}}{2}

    3. Sixthpower

       * .. math::

            \epsilon_{ij} = 2 \sqrt{\epsilon_{i} \times \epsilon_{j}}
            \frac{\sigma_i^3 \times \sigma_j^3}{\sigma_i^6 + \sigma_j^6}

       * .. math::

            \sigma_{ij} = \left( \frac{\sigma_{i}^6 \times
            \sigma_{j}^6}{2} \right)^{1/6}

       * .. math::

            r_{\text{c},ij} = \left(\frac{r_{\text{c},i}^6 \times
            r_{\text{c},j}^6}{2}\right)^{1/6}


    Parameters
    ----------
    epsilon : list of float
        Lennard-Jones epsilon parameter for a particle of type `i` and 'j'.
    sigma : list of float
        Lennard-Jones sigma parameter for a particle of type `i` and 'j'.
    rcut : list of float
        Lennard-Jones cut-off value for a particle of type `i` and 'j'.
    mixing :  string, optional
        Represents what kind of mixing that should be done.

    Returns
    -------
    out[0] : float
        The mixed ``epsilon_ij`` parameter.
    out[1] : float
        The mixed ``sigma_ij`` parameter.
    out[2] : float
        The mixed ``rcut_ij`` parameter.

    """
    if mixing == 'geometric':
        epsilon_ij = np.sqrt(epsilon[0] * epsilon[1])
        sigma_ij = np.sqrt(sigma[0] * sigma[1])
        rcut_ij = np.sqrt(rcut[0] * rcut[1])
    elif mixing == 'arithmetic':
        epsilon_ij = np.sqrt(epsilon[0] * epsilon[1])
        sigma_ij = 0.5 * (sigma[0] + sigma[1])
        rcut_ij = 0.5 * (rcut[0] + rcut[1])
    elif mixing == 'sixthpower':
        si3 = sigma[0]**3
        si6 = si3**2
        sj3 = sigma[1]**3
        sj6 = sj3**2
        avgs6 = 0.5 * (si6 + sj6)
        epsilon_ij = np.sqrt(epsilon[0] * epsilon[1]) * si3 * sj3 / avgs6
        sigma_ij = avgs6**(1.0 / 6.0)
        rcut_ij = (0.5 * (rcut[0]**6 + rcut[1]**6))**(1.0 / 6.0)
    else:
        epsilon_ij = 0.0
        sigma_ij = 0.0
        rcut_ij = 0.0
        logger.warning('No mixing requested. Cross interactions set to 0')
    return epsilon_ij, sigma_ij, rcut_ij
