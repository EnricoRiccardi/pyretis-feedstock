# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Definition of order parameters.

This package defines order parameters for use with PyRETIS.


Package structure
-----------------

Modules
~~~~~~~

orderparameter.py (:py:mod:`pyretis.orderparameter.orderparameter`)
    Defines the base class for order parameters and some simple
    example order parameters.

orderangle.py (:py:mod:`pyretis.orderparameter.orderangle`)
    Defines a class for an angle order parameter.

orderdihedral.py (:py:mod:`pyretis.orderparameter.orderdihedral`)
    Defines a class for a dihedral angle order parameter.

Important methods defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

order_factory (:py:func:`.order_factory`)
    A method to create order parameters from settings.

"""
from pyretis.core.common import generic_factory
from .orderparameter import (
    OrderParameter,
    Position,
    Velocity,
    Distance,
    Distancevel,
    PositionVelocity,
    DistanceVelocity,
    Permeability,
    PermeabilityMinusOffset,
    CompositeOrderParameter,
    expand_order_names,
    normalize_order_output,
    wrap_orderparameter,
)
from .orderangle import Angle
from .orderdihedral import Dihedral


def order_factory(settings):
    """Create order parameters according to the given settings.

    This function is included as a convenient way of setting up and
    selecting the order parameter.

    Parameters
    ----------
    settings : dict
        This defines how we set up and select the order parameter.

    Returns
    -------
    out : object like :py:class:`.OrderParameter`
        An object representing the order parameter.

    """
    factory_map = {
        'orderparameter': {
            'cls': OrderParameter
        },
        'position': {
            'cls': Position
        },
        'velocity': {
            'cls': Velocity
        },
        'distance': {
            'cls': Distance
        },
        'distancevel': {
            'cls': Distancevel
        },
        'positionvelocity': {
            'cls': PositionVelocity
        },
        'distancevelocity': {
            'cls': DistanceVelocity
        },
        'angle': {
            'cls': Angle
        },
        'dihedral': {
            'cls': Dihedral
        },
        'permeability': {
            'cls': Permeability
        },
        'permeabilityminusoffset': {
            'cls': PermeabilityMinusOffset
        },
    }
    return generic_factory(settings, factory_map, name='orderparameter')
