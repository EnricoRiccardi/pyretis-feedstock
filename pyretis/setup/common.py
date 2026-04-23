# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This module defines common methods for the settings handling.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

create_external (:py:func:`.create_external`)
    Method to create objects from settings.

check_settings (:py:func:`.check_settings`)
    Check that required simulation settings are actually given.

create_engine (:py:func:`.create_engine`)
    Method to create an engine from settings.

create_orderparameter (:py:func:`.create_orderparameter`)
    Method to create order parameters from settings.

create_potential (:py:func:`.create_potential`)
    Method to create a potential from settings.
"""
import functools
import logging
import os
import numpy as np
from pyretis.core.common import initiate_instance, import_from
from pyretis.engines import engine_factory
from pyretis.orderparameter import order_factory
from pyretis.orderparameter.orderparameter import CompositeOrderParameter
from pyretis.forcefield.factory import potential_factory
logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = ['create_external', 'check_settings', 'create_orderparameter',
           'create_engine', 'create_potential']


def _normalize_order_output(order, klass_name):
    """Return order-parameter data as a flat list for one frame."""
    msg = (
        f'Order parameter "{klass_name}" must return a flat list of values '
        'for a single frame. PyRETIS handles multiple frames as a list of '
        'lists outside calculate().'
    )
    if isinstance(order, np.ndarray):
        if order.ndim == 0:
            return [order.item()]
        if order.ndim == 1:
            return order.tolist()
        raise ValueError(msg)
    if np.isscalar(order):
        return [order.item() if hasattr(order, 'item') else order]
    try:
        values = list(order)
    except TypeError as exc:
        raise ValueError(msg) from exc
    if any(isinstance(value, (list, tuple, np.ndarray)) for value in values):
        raise ValueError(msg)
    return values


def _wrap_orderparameter(instance):
    """Normalize output from order parameters and collective variables."""
    if getattr(instance, '_pyretis_order_output_wrapped', False):
        return instance
    calculate = instance.calculate

    @functools.wraps(calculate)
    def wrapped(*args, **kwargs):
        order = calculate(*args, **kwargs)
        return _normalize_order_output(order, instance.__class__.__name__)

    instance.calculate = wrapped
    instance._pyretis_order_output_wrapped = True
    return instance


def check_settings(settings, required):
    """Check that required simulation settings are actually given.

    This method will look for required settings in the given
    `settings`. If one or more keys from the given `required` list of
    strings are not found, this method will return False. Otherwise,
    it will return True. Typically, an exception should be raised if
    False is returned, this is handled outside the method in case
    someone wants to add some magic handling of missing settings.

    Parameters
    ----------
    settings : dict
        This dict contains the given settings
    required : list of strings
        This list contains the settings that are required and which
        we will check the presence of.

    Returns
    -------
    result : boolean
        True if all required settings are present, False otherwise.
    not_found : list of strings
        There are the required settings we did not find.

    """
    result = True
    not_found = []
    for setting in required:
        if setting not in settings:
            result = False
            not_found.append(setting)
    return result, not_found


def create_external(settings, key, factory, required_methods,
                    key_settings=None):
    """Create external objects from settings.

    This method will handle the creation of objects from settings. The
    requested objects can be PyRETIS internals or defined in external
    modules.

    Parameters
    ----------
    settings : dict
        This dictionary contains the settings for the simulation.
    key : string
        The setting we are creating for.
    factory : callable
        A method to call that can handle the creation of internal
        objects for us.
    required_methods : list of strings
        The methods we need to have if creating an object from external
        files.
    key_settings : dict, optional
        This dictionary contains the settings for the specific key we
        are processing. If this is not given, we will try to obtain
        these settings by `settings[key]`. The reason why we make it
        possible to pass these as settings is in case we are processing
        a key which does not give a simple setting, but a list of settings.
        It that case `settings[key]` will give a list to process. That list
        is iterated somewhere else and `key_settings` can then be used to
        process these elements.

    Returns
    -------
    out : object
        This object represents the class we are requesting here.

    """
    if key_settings is None:
        try:
            key_settings = settings[key]
        except KeyError:
            logger.debug('No "%s" setting found. Skipping set-up', key)
            return None
    module = key_settings.get('module', None)
    klass = None
    try:
        klass = key_settings['class']
    except KeyError:
        logger.debug('No "class" setting for "%s" specified. Skipping set-up',
                     key)
        return None
    if module is None:
        return factory(key_settings)
    # Here we assume we are to load from a file. Before we import
    # we need to check that the path is ok or if we should include
    # the 'exe_path' from settings.
    # 1) Check if we can find the module:
    if os.path.isfile(module):
        obj = import_from(module, klass)
    else:
        if 'exe_path' in settings['simulation']:
            module = os.path.join(settings['simulation']['exe_path'],
                                  module)
            obj = import_from(module, klass)
        else:
            msg = f'Could not find module "{module}" for {key}!'
            raise ValueError(msg)
    # run some checks:
    for function in required_methods:
        objfunc = getattr(obj, function, None)
        if not objfunc:
            msg = f'Could not find method {klass}.{function}'
            logger.critical(msg)
            raise ValueError(msg)
        if not callable(objfunc):
            msg = f'Method {klass}.{function} is not callable!'
            logger.critical(msg)
            raise ValueError(msg)
    instance = initiate_instance(obj, key_settings)
    if key in ('orderparameter', 'collective-variable'):
        return _wrap_orderparameter(instance)
    return instance


def create_orderparameter(settings):
    """Create order parameters from settings.

    Parameters
    ----------
    settings : dict
        This dictionary contains the settings for the simulation.

    Returns
    -------
    out : object like :py:class:`.OrderParameter`
        This object represents the order parameter.

    """
    main_order = create_external(
        settings,
        'orderparameter',
        order_factory,
        ['calculate'],
    )
    if main_order is None:
        logger.info('No order parameter created')
        return None
    main_order = _wrap_orderparameter(main_order)
    logger.progress('Created main order parameter:\n%s', main_order)

    extra_cv = []
    order_settings = settings.get('collective-variable', [])
    for order_setting in order_settings:
        order = create_external(
            settings,
            'collective-variable',
            order_factory,
            ['calculate'],
            key_settings=order_setting
        )
        order = _wrap_orderparameter(order)
        logger.progress('Created additional collective variable:\n%s', order)
        extra_cv.append(order)
    if not extra_cv:
        return main_order
    all_order = [main_order] + extra_cv
    order = CompositeOrderParameter(order_parameters=all_order)
    order = _wrap_orderparameter(order)
    logger.progress('Composite order parameter:\n%s', order)
    return order


def create_engine(settings):
    """Create an engine from settings.

    Parameters
    ----------
    settings : dict
        This dictionary contains the settings for the simulation.

    Returns
    -------
    out : object like :py:class:`.EngineBase`
        This object represents the engine.

    """
    engine = create_external(settings, 'engine', engine_factory,
                             ['integration_step'])
    if not engine:
        raise ValueError('Could not create engine from settings!')
    return engine


def create_potential(settings, key_settings):
    """Create a potential from settings.

    Parameters
    ----------
    settings : dict
        This dictionary contains the settings for the simulation.
    key_settings : dict
        Settings for the potential we are creating.

    Returns
    -------
    out : object like :py:class:`.PotentialFunction`
        The object representing the potential function.

    """
    return create_external(settings, 'potential', potential_factory,
                           ['force', 'potential', 'potential_and_force'],
                           key_settings=key_settings)
