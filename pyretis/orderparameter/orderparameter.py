# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""This file contains classes to represent order parameters.

The order parameters are assumed to all be completely determined
by the system properties and they will all return at least one
value - the order parameter itself. The order parameters can also
return several order parameters which can be used for further analysis.

Important classes defined here
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OrderParameter (:py:class:`.OrderParameter`)
    Base class for the order parameters.

Position (:py:class:`.Position`)
    An order parameter equal to the position of a particle.

Velocity (:py:class:`.Velocity`)
    An order parameter equal to the velocity of a particle.

Distance (:py:class:`.Distance`)
    A class for a particle-particle distance order parameter.

Distancevel (:py:class:`.Distancevel`)
    A class for the rate of change in a particle-particle distance
    order parameter.

CompositeOrderParameter (:py:class:`.CompositeOrderParameter`)
    A class for an order parameter which is made up of several order
    parameters, i.e. of several objects like
    :py:class:`.OrderParameter`.

PositionVelocity (:py:class:`.PositionVelocity`)
    An order parameter which is equal to the composition of
    :py:class:`.Position` and
    :py:class:`.Velocity`.

DistanceVelocity (:py:class:`.DistanceVelocity`)
    An order parameter which is equal to the composition of
    :py:class:`.Distance` and :py:class:`.Distancevel`.

"""

from abc import abstractmethod
import functools
import logging
import numpy as np

logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())


__all__ = [
    "OrderParameter",
    "Position",
    "Permeability",
    "PermeabilityMinusOffset",
    "Velocity",
    "Distance",
    "Distancevel",
    "DistanceVelocity",
    "PositionVelocity",
    "CompositeOrderParameter",
    "expand_order_names",
    "normalize_order_output",
    "wrap_orderparameter",
]


def normalize_order_output(order, klass_name="unknown"):
    """Normalise the return value of ``calculate()`` to a flat list.

    This function accepts the raw return value from any ``calculate``
    implementation — a scalar, a 0-D or 1-D :py:class:`numpy.ndarray`,
    a list, or a tuple — and returns a plain Python list whose elements
    are all plain Python scalars (no nested lists or arrays).

    Parameters
    ----------
    order : scalar, list, tuple, or numpy.ndarray
        The raw return value from a ``calculate`` method.
    klass_name : str, optional
        Name of the order-parameter class; used only in error messages.

    Returns
    -------
    out : list of floats
        A flat list of scalar values.

    Raises
    ------
    ValueError
        If *order* is a >1-D array or contains nested sequences.

    """
    msg = (
        f'Order parameter "{klass_name}" must return a flat list of values '
        "for a single frame. PyRETIS handles multiple frames as a list of "
        "lists outside calculate()."
    )
    if isinstance(order, np.ndarray):
        if order.ndim == 0:
            return [order.item()]
        if order.ndim == 1:
            return order.tolist()
        raise ValueError(msg)
    if np.isscalar(order):
        return [order.item() if hasattr(order, "item") else order]
    try:
        values = list(order)
    except TypeError as exc:
        raise ValueError(msg) from exc
    if any(isinstance(value, (list, tuple, np.ndarray)) for value in values):
        raise ValueError(msg)
    return values


def expand_order_names(name, count, default_prefix="op", start_index=1):
    """Expand an order-parameter ``name`` into ``count`` per-element labels.

    The naming convention is shared by the main order parameter and any
    collective variable: each ``calculate()`` returns a list with
    ``count`` entries and this function produces one label per entry.

    Parameters
    ----------
    name : str, list/tuple of str, or None
        The user-provided name. When ``None`` the labels are always
        ``"<default_prefix>_<i>"``. When ``name`` is a string and
        ``count > 1`` the labels are ``"<name>_1", "<name>_2", ...``;
        when ``count == 1`` the bare ``name`` is used. When ``name`` is
        a list or tuple its length must equal ``count``.
    count : int
        Number of labels to produce (the number of values returned by
        ``calculate()`` for this order parameter or collective variable).
    default_prefix : str, optional
        Prefix used to build labels when ``name`` is ``None``. Use
        ``"op"`` for the main order parameter and ``"cv"`` for any
        collective variable.
    start_index : int, optional
        Starting index used when generating ``"<prefix>_<i>"`` labels.
        Useful for combining several blocks into one continuous index
        space.

    Returns
    -------
    out : list of str
        A list of length ``count``.

    Raises
    ------
    ValueError
        If ``count`` is not a positive integer, if ``name`` is a
        sequence whose length differs from ``count``, or if ``name`` is
        of an unsupported type.

    """
    if not isinstance(count, int) or count < 1:
        raise ValueError(
            f"count must be a positive integer (got {count!r})."
        )
    if name is None:
        return [
            f"{default_prefix}_{i}"
            for i in range(start_index, start_index + count)
        ]
    if isinstance(name, str):
        if count == 1:
            return [name]
        return [
            f"{name}_{i}"
            for i in range(start_index, start_index + count)
        ]
    if isinstance(name, (list, tuple)):
        if len(name) != count:
            raise ValueError(
                "Length of name list "
                f"({len(name)}) does not match the number of values "
                f"returned by the order parameter ({count})."
            )
        if not all(isinstance(item, str) for item in name):
            raise ValueError(
                "All entries in a name list must be strings."
            )
        return list(name)
    raise ValueError(
        "Order-parameter name must be a string, a list/tuple of "
        f"strings, or None (got {type(name).__name__})."
    )


def wrap_orderparameter(instance):
    """Wrap an OP instance so ``calculate()`` always returns a list.

    The wrapper is idempotent: calling this function twice on the same
    instance is safe and has no additional effect.

    Parameters
    ----------
    instance : object with a ``calculate`` method
        Any order-parameter object (internal or external).

    Returns
    -------
    instance : same object, mutated in place.

    """
    if getattr(instance, "_pyretis_order_output_wrapped", False):
        return instance
    calculate = instance.calculate

    @functools.wraps(calculate)
    def wrapped(*args, **kwargs):
        order = calculate(*args, **kwargs)
        return normalize_order_output(order, instance.__class__.__name__)

    instance.calculate = wrapped
    instance._pyretis_order_output_wrapped = True
    return instance


class OrderParameter:
    """Base class for order parameters.

    This class represents an order parameter and other collective
    variables. The order parameter is assumed to be a function
    that can uniquely be determined by the system object and its
    attributes.

    Attributes
    ----------
    description : string
        This is a short description of the order parameter.
    velocity_dependent : boolean
        This flag indicates whether or not the order parameter
        depends on the velocity direction. If so, we need to
        recalculate the order parameter when reversing trajectories.

    """

    def __init__(self, description="Generic order parameter", velocity=False):
        """Initialise the OrderParameter object.

        Parameters
        ----------
        description : string
            Short description of the order parameter.

        """
        self.description = description
        self.velocity_dependent = velocity
        if self.velocity_dependent:
            logger.debug(
                'Order parameter "%s" was marked as velocity dependent.',
                self.description,
            )

    @abstractmethod
    def calculate(self, system):
        """Calculate the order parameter (or collective variable).

        The same interface is used both for the main order parameter
        and for additional collective variables: each call returns
        the per-frame value(s) of the quantity. Implementations may
        return a **single scalar**, a :py:class:`numpy.ndarray`, a
        tuple, or a list. When loaded through
        :py:func:`pyretis.setup.common.create_orderparameter` or
        wrapped via :py:func:`.wrap_orderparameter`, the return value
        is normalised to a plain Python list of scalars by
        :py:func:`.normalize_order_output`. Downstream code therefore
        always sees ``system.order`` and ``phasepoint.order`` as a
        list whose first element is the progress coordinate used by
        path-sampling moves.

        Parameters
        ----------
        system : object like :py:class:`.System`
            This object contains the information needed to calculate
            the order parameter.

        Returns
        -------
        out : scalar, list, tuple, or numpy.ndarray
            The order-parameter (or collective-variable) value(s).
            After normalisation the first element is used as the
            progress coordinate in path-sampling simulations and any
            additional elements are stored alongside it.

        """
        return

    def __str__(self):
        """Return a simple string representation of the order parameter."""
        msg = [
            f'Order parameter: "{self.__class__.__name__}"',
            f"{self.description}",
        ]
        if self.velocity_dependent:
            msg.append("This order parameter is velocity dependent.")
        return "\n".join(msg)

    def load_restart_info(self, info):
        """Load the order parameter restart info."""

    def restart_info(self):
        """Save any mutable parameters for the restart."""


class Position(OrderParameter):
    """A positional order parameter.

    This class defines a very simple order parameter which is just
    the position of a given particle.

    Attributes
    ----------
    index : integer
        This is the index of the atom which will be used, i.e.
        ``system.particles.pos[index]`` will be used.
    dim : integer
        This is the dimension of the coordinate to use.
        0, 1 or 2 for 'x', 'y' or 'z'.
    periodic : boolean
        This determines if periodic boundaries should be applied to
        the position or not.

    """

    def __init__(self, index, dim="x", periodic=False, description=None):
        """Initialise the order parameter.

        Parameters
        ----------
        index : int
            This is the index of the atom we will use the position of.
        dim : string
            This select what dimension we should consider,
            it should equal 'x', 'y' or 'z'.
        periodic : boolean, optional
            This determines if periodic boundary conditions should be
            applied to the position.

        """
        if description is None:
            description = f"Position of particle {index} (dim: {dim})"
        super().__init__(description=description, velocity=False)
        self.periodic = periodic
        self.index = index
        self.dim = {"x": 0, "y": 1, "z": 2}.get(dim, None)
        if self.dim is None:
            msg = f"Unknown dimension {dim} requested"
            logger.critical(msg)
            raise ValueError(msg)

    def calculate(self, system):
        """Calculate the position order parameter.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The object containing the positions.

        Returns
        -------
        out : list of floats
            The position order parameter.

        """
        particles = system.particles
        pos = particles.pos[self.index]
        lamb = pos[self.dim]
        if self.periodic:
            lamb = system.box.pbc_coordinate_dim(lamb, self.dim)
        return [lamb]


class Permeability(Position):
    """A positional order parameter for calculating parmeability.

    This class defines a very simple order parameter which is just
    the position of a given particle, but it allows for the mirror move.

    Attributes
    ----------
    index : integer
        This is the index of the atom which will be used, i.e.
        ``system.particles.pos[index]`` will be used.
    dim : integer
        This is the dimension of the coordinate to user:
        0, 1 or 2 for 'x', 'y' or 'z'.
    periodic : boolean
        This determines if periodic boundaries should be applied to
        the position or not.
    offset : float
        The offset to apply before returning, l(x) = pos(x+offset).
        This allows for effectively moving the periodic boundary position.
        Should be defined in the same space as this order parameter
        (relative or not).
    mirror_pos : float
        The position of the mirror plane around which we mirror when the
        mirror function is called. Should be defined in the same space as this
        order parameter (relative or not).
    relative : boolean
        If we should map the position to be relative to the box-vector.
        Defaults to True as that is the only correct option for NPT.
        For NVT sampling it can be set to False.

    """

    def __init__(
        self, index, dim="z", offset=0.0, mirror_pos=0.0, relative=True
    ):
        """Initialise the order parameter.

        Parameters
        ----------
        index : int
            This is the index of the atom we will use the position of.
        dim : string
            This select what dimension we should consider,
            it should equal 'x', 'y' or 'z'.
        offset : float,
            The offset to add to the position, used to change the
            effective location of the periodic boundary.
        mirror_pos : float
            The value of the position to mirror around. Allowing for the mirror
            move.
        relative : boolean
            This determines if the position is returned as a relative value
            to the box size.

        """
        description = f"Permeability position of particle {index} (dim: {dim})"
        super().__init__(
            index=index, dim=dim, periodic=True, description=description
        )

        self.offset = offset
        self.mirror_pos = mirror_pos
        self.relative = relative
        self._mirror = False
        if relative and abs(offset) >= 1:
            raise ValueError(
                "Mapping to relative space, but offset is not "
                f"between -1 and 1, offset was : {offset}"
            )
        if relative and abs(mirror_pos) >= 1:
            raise ValueError(
                "Mapping to relative space, but mirror_pos is"
                "not between -1 and 1, mirror_pos was : "
                f"{mirror_pos}"
            )

    def calculate_s(self, system):
        """Calculate a value that we the alter in the `calculate` function.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The object containing the positions.

        Returns
        -------
        out : list of floats
            The position order parameter.

        """
        return super().calculate(system)

    def calculate(self, system):
        """Calculate the `compute_s` order parameter and alters it.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The object containing the positions.

        Returns
        -------
        out : list of floats
            The position order parameter.

        """
        # Just get the position and unwrap it
        values = self.calculate_s(system)
        value = values.pop(0)

        # Now map to relative space or not
        if self.relative and system.box is not None:
            value /= system.box.length[self.dim]
            box = 1
            l_box = system.box.low[self.dim] / system.box.length[self.dim]
        elif system.box is not None:
            box = system.box.length[self.dim]
            l_box = system.box.low[self.dim]
        else:
            box = None
            l_box = 0
        # Map to mirrored space if required
        if self._mirror:
            value = 2 * self.mirror_pos - value
        # Apply lowest_box
        value -= l_box

        # Add offset
        value += self.offset

        # Wrap box
        if box is not None:
            value = np.mod(value, box)
        # Reapply lowest box if  not in relative space

        if not self.relative:
            value += l_box
        return [value, self.index, (-1) ** self._mirror] + values

    def mirror(self):
        """Swap this object around the mirror plane."""
        # Make mirror of the OP
        self._mirror = not self._mirror

    def restart_info(self):
        """Return the mutable attributes for the restart."""
        info = {"index": self.index, "mirror": self._mirror}

        return info

    def load_restart_info(self, info):
        """Load the mutable attributes from restart."""
        self.index = info["index"]
        self._mirror = info["mirror"]


class PermeabilityMinusOffset(Permeability):
    """A positional order parameter for calculating parmeability - offset."""

    def calculate(self, system):
        """Calculate the permeability op and subtract the offset."""
        # Just get the permeability and unwrap it
        values = super().calculate(system)
        return [values[0] - self.offset] + values[1:]


class Velocity(OrderParameter):
    """Initialise the order parameter.

    This class defines a very simple order parameter which is just
    the velocity of a given particle.

    Attributes
    ----------
    index : integer
        This is the index of the atom which will be used, i.e.
        ``system.particles.vel[index]`` will be used.
    dim : integer
        This is the dimension of the coordinate to use.
        0, 1 or 2 for 'x', 'y' or 'z'.

    """

    def __init__(self, index, dim="x"):
        """Initialise the order parameter.

        Parameters
        ----------
        index : int
            This is the index of the atom we will use the velocity of.
        dim : string
            This select what dimension we should consider,
            it should equal 'x', 'y' or 'z'.

        """
        txt = f"Velocity of particle {index} (dim: {dim})"
        super().__init__(description=txt, velocity=True)
        self.index = index
        self.dim = {"x": 0, "y": 1, "z": 2}.get(dim, None)
        if self.dim is None:
            logger.critical("Unknown dimension %s requested", dim)
            raise ValueError

    def calculate(self, system):
        """Calculate the velocity order parameter.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The object containing the velocities.

        Returns
        -------
        out : list of floats
            The velocity order parameter.

        """
        return [system.particles.vel[self.index][self.dim]]


def _verify_pair(index):
    """Check that the given index contains a pair."""
    try:
        if len(index) != 2:
            msg = (
                "Wrong number of atoms for pair definition. "
                f"Expected 2 got {len(index)}"
            )
            logger.error(msg)
            raise ValueError(msg)
    except TypeError as err:
        msg = "Atom pair should be defined as a tuple/list of integers."
        logger.error(msg)
        raise TypeError(msg) from err


class Distance(OrderParameter):
    """A distance order parameter.

    This class defines a very simple order parameter which is just
    the scalar distance between two particles.

    Attributes
    ----------
    index : tuple of integers
        These are the indices used for the two particles.
        `system.particles.pos[index[0]]` and
        `system.particles.pos[index[1]]` will be used.
    periodic : boolean
        This determines if periodic boundaries should be applied to
        the distance or not.

    """

    def __init__(self, index, periodic=True):
        """Initialise order parameter.

        Parameters
        ----------
        index : tuple of ints
            This is the indices of the atom we will use the position of.
        periodic : boolean, optional
            This determines if periodic boundary conditions should be
            applied to the position.

        """
        _verify_pair(index)
        pbc = "Periodic" if periodic else "Non-periodic"
        txt = f"{pbc} distance, particles {index[0]} and {index[1]}"
        super().__init__(description=txt, velocity=False)
        self.periodic = periodic
        self.index = index

    def calculate(self, system):
        """Calculate the order parameter.

        Here, the order parameter is just the distance between two
        particles.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The object containing the positions and box used for the
            calculation.

        Returns
        -------
        out : list of floats
            The distance order parameter.

        """
        particles = system.particles
        delta = particles.pos[self.index[1]] - particles.pos[self.index[0]]
        if self.periodic:
            delta = system.box.pbc_dist_coordinate(delta)
        lamb = np.sqrt(np.dot(delta, delta))
        return [lamb]


class Distancevel(OrderParameter):
    """A rate of change of the distance order parameter.

    This class defines a very simple order parameter which is just
    the time derivative of the scalar distance between two particles.

    Attributes
    ----------
    index : tuple of integers
        These are the indices used for the two particles.
        `system.particles.pos[index[0]]` and
        `system.particles.pos[index[1]]` will be used.
    periodic : boolean
        This determines if periodic boundaries should be applied to
        the distance or not.

    """

    def __init__(self, index, periodic=True):
        """Initialise the order parameter.

        Parameters
        ----------
        index : tuple of ints
            This is the indices of the atom we will use the position of.
        periodic : boolean, optional
            This determines if periodic boundary conditions should be
            applied to the position.

        """
        _verify_pair(index)
        pbc = "Periodic" if periodic else "Non-periodic"
        txt = (
            f"{pbc} rate-of-change-distance, "
            f"particles {index[0]} and {index[1]}"
        )
        super().__init__(description=txt, velocity=True)
        self.periodic = periodic
        self.index = index

    def calculate(self, system):
        """Calculate the order parameter.

        Here, the order parameter is just the distance between two
        particles.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The object containing the positions and box used for the
            calculation.

        Returns
        -------
        out : list of floats
            The rate-of-change of the distance order parameter.

        """
        particles = system.particles
        delta = particles.pos[self.index[1]] - particles.pos[self.index[0]]
        if self.periodic:
            delta = system.box.pbc_dist_coordinate(delta)
        lamb = np.sqrt(np.dot(delta, delta))
        # Add the velocity as an additional collective variable:
        delta_v = particles.vel[self.index[1]] - particles.vel[self.index[0]]
        cv1 = np.dot(delta, delta_v) / lamb
        return [cv1]


class CompositeOrderParameter(OrderParameter):
    """A composite order parameter.

    This class represents a composite order parameter. It does not
    actually calculate order parameters itself, but it has references
    to several objects like :py:class:`.OrderParameter` which it can
    use to obtain the order parameters. Note that the first one of
    these objects will be interpreted as the main progress coordinate
    in path sampling simulations.

    Attributes
    ----------
    extra : list of objects like :py:class:`OrderParameter`
        This is a list of order parameters to calculate.

    """

    def __init__(self, order_parameters=None):
        """Set up the composite order parameter.

        Parameters
        ----------
        order_parameters : list of objects like :py:class:`.OrderParameter`
            A list of order parameters we can add.

        """
        super().__init__(
            description="Combined order parameter", velocity=False
        )
        self.order_parameters = []
        if order_parameters is not None:
            for order_function in order_parameters:
                self.add_orderparameter(order_function)

    def calculate(self, system):
        """Calculate the main order parameter and return it.

        This is defined as a method just to ensure that at least this
        method will be defined in the different order parameters.

        Parameters
        ----------
        system : object like :py:class:`.System`
            This object contains the information needed to calculate
            the order parameter.

        Returns
        -------
        out : list of floats
            The order parameter(s). The first order parameter returned
            is assumed to be the progress coordinate for path sampling
            simulations.

        """
        all_order = []
        for order_function in self.order_parameters:
            all_order.extend(order_function.calculate(system))
        return all_order

    def mirror(self):
        """Mirrors all of the functions that allow it."""
        order_p = self.order_parameters[0]
        op_mirror_func = getattr(order_p, "mirror", None)
        if op_mirror_func is not None:
            op_mirror_func()
        else:
            msg = "Attempting a mirror move, but orderparameter of \n class"
            msg += f" '{type(order_p).__name__}' does not have the function"
            msg += " 'mirror()'.\n"
            msg += "Please use an OP of type 'Permeability' or implement your"
            msg += " own mirror() function"
            logger.warning(msg)
        # This is safe as compound OPs should always have more than 1 OP.
        for order_function in self.order_parameters[1:]:
            mirror_func = getattr(order_function, "mirror", None)
            if mirror_func is not None:
                mirror_func()

    def restart_info(self):
        """Return the mutable attributes for restart."""
        return [op.restart_info() for op in self.order_parameters]

    def load_restart_info(self, info):
        """Load the mutable attributes for restart."""
        for i, op_info in enumerate(info):
            self.order_parameters[i].load_restart_info(op_info)

    def add_orderparameter(self, order_function):
        """Add an extra order parameter to calculate.

        Parameters
        ----------
        order_function : object like :py:class:`.OrderParameter`
            An object we can use to calculate the order parameter.

        Returns
        -------
        out : boolean
            Return True if we added the function, False otherwise.

        """
        # We check that the ``calculate`` method is present and callable.
        for func in ("calculate",):
            objfunc = getattr(order_function, func, None)
            name = order_function.__class__.__name__
            if not objfunc:
                msg = f'Missing method "{func}" in order parameter {name}'
                logger.error(msg)
                raise ValueError(msg)
            if not callable(objfunc):
                msg = f'"{func}" in order parameter {name} is not callable!'
                raise ValueError(msg)
        self.velocity_dependent |= order_function.velocity_dependent
        if self.velocity_dependent:
            logger.debug(
                'Order parameter "%s" was marked as velocity dependent.',
                self.description,
            )
        self.order_parameters.append(order_function)
        return True

    @property
    def index(self):
        """Get only the index that is tracked by the orderparameter."""
        return self.order_parameters[0].index

    @index.setter
    def index(self, var):
        """Set only the index that is tracked by the orderparameter."""
        self.order_parameters[0].index = var

    def __str__(self):
        """Return a simple string representation of the order parameter."""
        txt = ["Order parameter, combination of:"]
        for i, order in enumerate(self.order_parameters):
            txt.append(f"{i}: {str(order)}")
        msg = "\n".join(txt)
        return msg


class PositionVelocity(CompositeOrderParameter):
    """An order parameter equal to the position & velocity of a given atom."""

    def __init__(self, index, dim="x", periodic=False):
        """Initialise the order parameter.

        Parameters
        ----------
        index : int
            This is the index of the atom we will use the position
            and velocity of.
        dim : string
            This select what dimension we should consider,
            it should equal 'x', 'y' or 'z'.
        periodic : boolean, optional
            This determines if periodic boundary conditions should be
            applied to the position.

        """
        position = Position(index, dim=dim, periodic=periodic)
        velocity = Velocity(index, dim=dim)
        orderparameters = [position, velocity]
        super().__init__(order_parameters=orderparameters)


class DistanceVelocity(CompositeOrderParameter):
    """An order parameter equal to a distance and its rate of change."""

    def __init__(self, index, periodic=True):
        """Initialise the order parameter.

        Parameters
        ----------
        index : tuple of integers
            These are the indices used for the two particles.
            `system.particles.pos[index[0]]` and
            `system.particles.pos[index[1]]` will be used.
        periodic : boolean, optional
            This determines if periodic boundary conditions should be
            applied to the position.

        """
        position = Distance(index, periodic=periodic)
        velocity = Distancevel(index, periodic=periodic)
        orderparameters = [position, velocity]
        super().__init__(order_parameters=orderparameters)
