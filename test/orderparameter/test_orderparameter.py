# -*- coding: utf-8 -*-
# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the order parameter classes from pyretis.orderparameter."""

import logging
import pytest
import numpy as np
from pyretis.orderparameter import order_factory
from pyretis.orderparameter import (
    OrderParameter,
    Position,
    Velocity,
    PositionVelocity,
    Distance,
    Distancevel,
    DistanceVelocity,
    CompositeOrderParameter,
    Angle,
    Dihedral,
    Permeability,
    PermeabilityMinusOffset,
    expand_order_names,
    normalize_order_output,
    wrap_orderparameter,
)
from pyretis.core import System, create_box, Particles

logging.disable(logging.CRITICAL)


class SimpleOrder(OrderParameter):
    """An order parameter which equals the system temperature."""

    def __init__(self):
        super().__init__(description="Simple order parameter")

    def calculate(self, system):
        return [system.temperature["set"]]


class SimpleOrderTemp(OrderParameter):
    """An order parameter which equals the system temperature."""

    def __init__(self):
        super().__init__(description="Simple order parameter", velocity=True)

    def calculate(self, system):
        return [system.temperature["set"] * -1]


class SimpleOrderFaulty:  # pylint: disable=too-few-public-methods
    """An order parameter which is faulty - missing calculate."""

    def calculate_all(self, system):  # pylint: disable=no-self-use
        """Just return the set temperature."""
        return [system.temperature["set"]]


class SimpleOrderFaulty2:  # pylint: disable=too-few-public-methods
    """An order parameter which is faulty - missing callable calculate."""

    calculate = 100


class TestOrderGeneric:
    """Test that we can create a class and define some parameters."""

    def test_simple_order(self):
        """Test that we can create a very simple order parameter."""
        system = System(temperature=123.0, units="lj", box=None)
        order = SimpleOrder()
        correct = [123.0, "lj"]
        val = order.calculate(system)
        assert val[0] == pytest.approx(correct[0])
        assert len(val) == 1

        vals = order.calculate(system)
        assert vals[0] == pytest.approx(correct[0])
        assert len(vals) == 1


def create_system(ndim, npart, periodic=False):
    """Create a simple system for testing."""
    if periodic:
        box = create_box(low=[0] * ndim, high=[1] * ndim)
    else:
        box = create_box(periodic=[False] * ndim)
    system = System(box=box)
    system.particles = Particles(system.get_dim())
    for _ in range(npart):
        pos = np.random.random(box.dim)
        vel = np.random.random(box.dim)
        system.add_particle(name="Ar", pos=pos, vel=vel)
    return system, box


class TestOrderPosition:
    """Run the tests for the Position class."""

    def _check_order_parameter(self, orderp, correct, system, idim, ndim):
        """Verify the order parameter."""
        for orderpi, correcti in zip(orderp, correct):
            if idim > ndim - 1:
                with pytest.raises(IndexError):
                    orderpi.calculate(system)
            else:
                for i, j in zip(orderpi.calculate(system), correcti):
                    assert i == pytest.approx(j)

    @staticmethod
    def _get_order(index, xdim, periodic=False):
        """Return the initiated order parameters."""
        orderp = [
            Position(index, dim=xdim, periodic=periodic),
            Velocity(index, dim=xdim),
            PositionVelocity(index, dim=xdim, periodic=periodic),
        ]
        return orderp

    def _get_correct(self, system, index, idim, periodic=False):
        """Return the correct order parameters."""
        correct = [
            self._correct_order1(system, index, idim, periodic=periodic),
            self._correct_order2(system, index, idim),
        ]
        correct.append(correct[0] + correct[1])
        return correct

    def test_without_pbc(self):
        """Test the order parameters without periodic boundaries."""
        periodic = False
        for npart, index in zip((1, 10), (0, 2)):
            for ndim in [1, 2, 3]:
                system, _ = create_system(ndim, npart, periodic=periodic)
                for idim, xdim in enumerate(("x", "y", "z")):
                    orderp = self._get_order(index, xdim, periodic=periodic)
                    correct = self._get_correct(
                        system, index, idim, periodic=periodic
                    )
                    self._check_order_parameter(
                        orderp, correct, system, idim, ndim
                    )

    def test_with_pbc(self):
        """Test the order parameters with periodic boundaries."""
        periodic = True
        for npart, index in zip((1, 10), (0, -1)):
            for ndim in [1, 2, 3]:
                # Just test for some displacements:
                for disp in [0.0, 1.5, -1.5, 100.0, -100.0]:
                    system, _ = create_system(ndim, npart, periodic=periodic)
                    system.particles.pos += (
                        np.ones_like(system.particles.pos) * disp
                    )
                    for idim, xdim in enumerate(("x", "y", "z")):
                        orderp = self._get_order(
                            index, xdim, periodic=periodic
                        )
                        correct = self._get_correct(
                            system, index, idim, periodic=periodic
                        )
                        self._check_order_parameter(
                            orderp, correct, system, idim, ndim
                        )

    @staticmethod
    def _correct_order1(system, index, idim, periodic=False):
        """The correct position order parameter."""
        try:
            pos = system.particles.pos[index][idim]
            if periodic:
                box = system.box
                return [box.pbc_coordinate_dim(pos, idim)]
            return [pos]
        except IndexError:
            return [None]

    @staticmethod
    def _correct_order2(system, index, idim):
        """The correct velocity order parameter."""
        try:
            return [system.particles.vel[index][idim]]
        except IndexError:
            return [None]

    def test_init_fail(self):
        """Check that the initiation fails if we supply strange input."""
        with pytest.raises(ValueError):
            Position(0, dim="a")
        with pytest.raises(ValueError):
            Velocity(0, dim="pingu")
        with pytest.raises(ValueError):
            PositionVelocity(123, dim="chonky")


class TestOrderDistance:
    """Run the tests for the Distance class."""

    def _check_order_parameter(self, orderp, correct, system):
        """Verify the order parameter."""
        for orderpi, correcti in zip(orderp, correct):
            for i, j in zip(orderpi.calculate(system), correcti):
                assert i == pytest.approx(j)

    @staticmethod
    def _get_order(index, periodic=False):
        """Return the initiated order parameters."""
        orderp = [
            Distance(index, periodic=periodic),
            Distancevel(index, periodic=periodic),
            DistanceVelocity(index, periodic=periodic),
        ]
        return orderp

    def _get_correct(self, system, index, periodic=False):
        """Return the correct order parameters."""
        correct = [
            self._correct_order1(system, index, periodic=periodic),
            self._correct_order2(system, index, periodic=periodic),
        ]
        correct.append(correct[0] + correct[1])
        return correct

    def test_two_particles(self):
        """Test the distance order parameter without pbc."""
        # Test for a one-particle system:
        index = (0, 1)
        for ndim in [1, 2, 3]:
            system, _ = create_system(ndim, 2, periodic=False)
            orderp = self._get_order(index, periodic=False)
            correct = self._get_correct(system, index, periodic=False)
            self._check_order_parameter(orderp, correct, system)

    def test_two_particles_pbs(self):
        """Test the distance order parameter with pbc."""
        # Test for a one-particle system:
        index = (0, 1)
        for disp in [0.0, 1.5, -1.5, 100.0, -100.0]:
            for ndim in [1, 2, 3]:
                system, box = create_system(ndim, 2, periodic=True)
                for i in index:
                    system.particles.pos[i] = (
                        np.random.random(box.dim) + np.ones(box.dim) * disp
                    )
                orderp = self._get_order(index, periodic=True)
                correct = self._get_correct(system, index, periodic=True)
                self._check_order_parameter(orderp, correct, system)

    @staticmethod
    def _correct_order1(system, index, periodic=False):
        """The correct position order parameter."""
        try:
            i, j = index
            delta = system.particles.pos[j] - system.particles.pos[i]
            if periodic:
                box = system.box
                delta = box.pbc_dist_coordinate(delta)
            return [np.sqrt(np.dot(delta, delta))]
        except IndexError:
            return [None]

    @staticmethod
    def _correct_order2(system, index, periodic=False):
        """The correct velocity order parameter."""
        try:
            i, j = index
            delta = system.particles.pos[j] - system.particles.pos[i]
            delta_v = system.particles.vel[j] - system.particles.vel[i]
            if periodic:
                box = system.box
                delta = box.pbc_dist_coordinate(delta)
            return [np.dot(delta, delta_v) / np.sqrt(np.dot(delta, delta))]
        except IndexError:
            return [None]

    def test_init_fail(self):
        """Check that the initiation fails if we supply strange input."""
        inputs = [0, [0], (0,), (0, 1, 2)]
        errors = [TypeError, ValueError, ValueError, ValueError]
        klasses = (Distance, Distancevel, DistanceVelocity)
        for cls in klasses:
            for i, j in zip(inputs, errors):
                with pytest.raises(j):
                    cls(i)


def water_molecule(box):
    """Return a simple system with a single water molecule."""
    system = System(box=box)
    system.particles = Particles(system.get_dim())
    system.add_particle(name="O", pos=np.array([0.230, 0.628, 0.113]))
    system.add_particle(name="H", pos=np.array([0.137, 0.626, 0.150]))
    system.add_particle(name="H", pos=np.array([0.231, 0.589, 0.021]))
    return system


def triangle():
    """Return a 2D system with particles in a tright riangle."""
    box = create_box(periodic=[False, False])
    system = System(box=box)
    system.particles = Particles(system.get_dim())
    system.add_particle(name="X", pos=np.array([0.0, 0.0]))
    system.add_particle(name="X", pos=np.array([1.0, 0.0]))
    system.add_particle(name="X", pos=np.array([0.0, 1.0]))
    angles = [
        ((1, 0, 2), 90.0),
        ((0, 1, 2), 45.0),
        ((0, 2, 1), 45.0),
    ]
    return system, angles


class TestOrderAngle:
    """Run the tests for the Angle class."""

    def test_without_pbc(self):
        """Test the angle order parameter without pbc."""
        orderp = Angle((1, 0, 2), periodic=False)
        box = create_box(periodic=[False] * 3)
        # Test angle for the SPC water geometry.
        system = water_molecule(box)
        angle = orderp.calculate(system)[0]
        assert np.degrees(angle) == pytest.approx(109.984398, abs=1e-3)

    def test_with_pbc(self):
        """Test the angle order parameter with pbc."""
        orderp = Angle((1, 0, 2), periodic=True)
        # Test angle for the SPC water geometry.
        box = create_box(periodic=[True, True, True], cell=[1.0, 1.0, 1.0])
        system = water_molecule(box)
        angle = orderp.calculate(system)[0]
        assert np.degrees(angle) == pytest.approx(109.984398, abs=1e-3)

    def test_triangle(self):
        """Test the angle order parameter for a 2D case."""
        system, angles = triangle()
        for idx, correct in angles:
            orderp = Angle(idx, periodic=False)
            angle = orderp.calculate(system)[0]
            assert np.degrees(angle) == pytest.approx(correct)

    def test_initiate_fail(self):
        """Test that we fail if we give incorrect number of indices."""
        with pytest.raises(TypeError):
            Angle(0, periodic=False)
        with pytest.raises(ValueError):
            Angle((0,), periodic=False)
        with pytest.raises(ValueError):
            Angle((0, 1), periodic=False)
        with pytest.raises(ValueError):
            Angle((0, 1, 2, 3), periodic=False)

    def test_special_cases(self):
        """Test the angle order parameter for some special cases:

        1. The angle between (1, 0, 0) and (0, 1, 0) -> pi/2.
        2. The angle between (1, 0, 0) and (1, 0, 0) -> 0.
        3. The angle between (1, 0, 0) and (-1, 0, 0) -> -pi.

        """
        orderp = Angle((0, 1, 2), periodic=False)
        test_cases = [
            {
                "angle": 0.5 * np.pi,
                "pos": np.array([
                    [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]
                ]),
            },
            {
                "angle": 0.0,
                "pos": np.array([
                    [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]
                ]),
            },
            {
                "angle": np.pi,
                "pos": np.array([
                    [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [-1.0, 0.0, 0.0]
                ]),
            },
        ]
        system, _ = create_system(3, 3, periodic=False)
        for case in test_cases:
            system.particles.pos = case["pos"]
            angle = orderp.calculate(system)[0]
            assert angle == pytest.approx(case["angle"])


class TestOrderDihedral:
    """Run the tests for the Dihedral class."""

    test_cases = [
        {
            "angle": 180.0,
            "pos": np.array([
                [0.0, 1.0, 0.0], [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0], [1.0, -1.0, 0.0],
            ]),
        },
        {
            "angle": 0.0,
            "pos": np.array([
                [0.0, 1.0, 0.0], [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0], [1.0, 1.0, 0.0],
            ]),
        },
        {
            "angle": -90.0,
            "pos": np.array([
                [0.0, 0.0, 1.0], [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0], [1.0, 1.0, 0.0],
            ]),
        },
        {
            "angle": 90.0,
            "pos": np.array([
                [0.0, 0.0, -1.0], [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0], [1.0, 1.0, 0.0],
            ]),
        },
        {
            "angle": -60.0127,
            "pos": np.array(
                [
                    [0.354, -2.210, -7.248],
                    [-0.290, -2.221, -6.483],
                    [0.472, -2.265, -5.191],
                    [1.036, -3.090, -5.164],
                ]
            ),
        },
        {
            "angle": 60.0319,
            "pos": np.array(
                [
                    [0.354, -2.210, -7.248],
                    [-0.290, -2.221, -6.483],
                    [0.472, -2.265, -5.191],
                    [1.058, -1.458, -5.122],
                ]
            ),
        },
        {
            "angle": 0.0,
            "pos": np.array(
                [
                    [1.499, -0.043, 0.000],
                    [2.055, 1.361, 0.000],
                    [3.481, 1.470, 0.000],
                    [3.898, 0.528, 0.000],
                ]
            ),
        },
        {
            "angle": -59.365971,
            "pos": np.array(
                [
                    [0.039, -0.028, 0.000],
                    [1.499, -0.043, 0.000],
                    [1.956, -0.866, -1.217],
                    [1.571, -1.903, -1.181],
                ]
            ),
        },
        {
            "angle": 60.833130,
            "pos": np.array(
                [
                    [0.039, -0.028, 0.000],
                    [1.499, -0.043, 0.000],
                    [1.956, -0.866, -1.217],
                    [1.610, -0.425, -2.172],
                ]
            ),
        },
        {
            "angle": -62.290916,
            "pos": np.array(
                [
                    [-0.543, -0.938, 0.000],
                    [0.039, -0.028, 0.000],
                    [1.499, -0.043, 0.000],
                    [1.847, -0.534, 0.928],
                ]
            ),
        },
    ]

    def test_without_pbc(self):
        """Test the angle order parameter without pbc."""
        orderp = Dihedral((3, 2, 1, 0), periodic=False)
        system, _ = create_system(3, 4, periodic=False)
        # Test some pre-defined cases:
        for case in self.test_cases:
            system.particles.pos = case["pos"]
            angle = orderp.calculate(system)[0]
            angle_deg = np.degrees(angle)  # pylint: disable=no-member
            assert angle_deg == pytest.approx(case["angle"], abs=1e-4)

    def test_with_pbc(self):
        """Test the angle order parameter with pbc."""
        orderp = Dihedral((3, 2, 1, 0), periodic=True)
        system, _ = create_system(3, 4, periodic=True)
        # Define a new box for this test:
        box = create_box(periodic=[True, True, True], cell=[8.0, 8.0, 8.0])
        system.box = box
        # Test same cases, just displaced.
        for case in self.test_cases:
            displace = np.ones_like(case["pos"]) * 9
            system.particles.pos = case["pos"] + displace
            angle = orderp.calculate(system)[0]
            angle_deg = np.degrees(angle)  # pylint: disable=no-member
            assert angle_deg == pytest.approx(case["angle"], abs=1e-4)

    def test_order(self):
        """Test if we get the same angle if we reverse indices."""
        order1 = Dihedral((0, 1, 2, 3), periodic=False)
        order2 = Dihedral((3, 2, 1, 0), periodic=False)
        system, _ = create_system(3, 4, periodic=False)
        for _ in range(3):
            system.particles.pos = np.random.rand(4, 3)
            angle1 = order1.calculate(system)[0]
            angle2 = order2.calculate(system)[0]
            assert angle1 == pytest.approx(angle2)

    def test_initiate_fail(self):
        """Test that we fail if we give incorrect number of indices."""
        with pytest.raises(TypeError):
            Dihedral(0, periodic=False)
        with pytest.raises(ValueError):
            Dihedral((0,), periodic=False)
        with pytest.raises(ValueError):
            Dihedral((0, 1), periodic=False)
        with pytest.raises(ValueError):
            Dihedral((0, 1, 2), periodic=False)
        with pytest.raises(ValueError):
            Dihedral((0, 1, 2, "tre"), periodic=False)


class TestOrderFactory:
    """Test the order factory."""

    def test_factory(self):
        """Test that we can create order parameters with the factory."""
        test_cases = [
            {
                "setting": {"class": "orderparameter"},
                "class": OrderParameter,
            },
            {
                "setting": {"class": "OrderPARAMetEr"},
                "class": OrderParameter,
            },
            {
                "setting": {"class": "position", "index": 0},
                "class": Position,
            },
            {
                "setting": {"class": "velocity", "index": 0},
                "class": Velocity,
            },
            {
                "setting": {"class": "distance", "index": (0, 1)},
                "class": Distance,
            },
            {
                "setting": {"class": "distancevel", "index": (0, 1)},
                "class": Distancevel,
            },
            {
                "setting": {"class": "PositionVelocity", "index": 0},
                "class": PositionVelocity,
            },
            {
                "setting": {"class": "distancevelocity", "index": (0, 1)},
                "class": DistanceVelocity,
            },
            {
                "setting": {"class": "angle", "index": (0, 1, 2)},
                "class": Angle,
            },
            {
                "setting": {"class": "dihedral", "index": (0, 1, 2, 3)},
                "class": Dihedral,
            },
        ]
        for case in test_cases:
            orderp = order_factory(case["setting"])
            assert isinstance(orderp, case["class"])


class TestCollection:
    """Test that we can create collections of order parameters."""

    def test_init(self):
        """Test creation of an object."""
        orderp = CompositeOrderParameter()
        cv1 = Distance((0, 1), periodic=False)
        orderp.add_orderparameter(cv1)
        cv2 = Distance((1, 2), periodic=False)
        orderp.add_orderparameter(cv2)
        cv3 = Distance((0, 2), periodic=False)
        orderp.add_orderparameter(cv3)
        orderp2 = CompositeOrderParameter(order_parameters=[cv1, cv2, cv3])
        for i, j in zip(orderp.order_parameters, orderp2.order_parameters):
            assert i is j
        system, _ = create_system(3, 3, periodic=False)
        system.particles.pos = np.array(
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.1], [0.0, 0.0, 0.3]]
        )
        order = orderp.calculate(system)
        correct = [0.1, 0.2, 0.3]
        assert len(order) == len(correct)
        for i, j in zip(order, correct):
            assert i == pytest.approx(j)

    def test_faulty_input(self):
        """Test if we supply faulty input."""
        orderp = CompositeOrderParameter()

        cv1 = SimpleOrder()
        orderp.add_orderparameter(cv1)
        assert cv1 is orderp.order_parameters[0]
        # Try to add some faulty order parameters.
        cv2 = SimpleOrderFaulty()
        with pytest.raises(ValueError):
            orderp.add_orderparameter(cv2)
        cv3 = SimpleOrderFaulty2()
        with pytest.raises(ValueError):
            orderp.add_orderparameter(cv3)

    def test_velocity_dependence(self):
        """Test that the combined is marked as velocity dependent."""
        orderp = CompositeOrderParameter()

        cv1 = SimpleOrder()
        orderp.add_orderparameter(cv1)
        assert not orderp.velocity_dependent

        cv2 = SimpleOrder()
        orderp.add_orderparameter(cv2)
        assert not orderp.velocity_dependent

        cv3 = SimpleOrderTemp()
        orderp.add_orderparameter(cv3)
        assert orderp.velocity_dependent

        cv4 = SimpleOrder()
        orderp.add_orderparameter(cv4)
        assert orderp.velocity_dependent


class TestPermeability:
    def setup_method(self):
        self.system, self.box = create_system(2, 2, periodic=True)
        self.op = Permeability(index=0, dim="x", relative=False)
        self.x = self.system.particles.pos[0, 0]

    def test_non_relative_offset(self):
        with pytest.raises(ValueError, match="offset"):
            Permeability(index=0, dim="x", offset=-1.01)
        # Test that we do pss if not relative
        _ = Permeability(index=0, dim="x", offset=-1.01, relative=False)

    def test_non_relative_mirror(self):
        with pytest.raises(ValueError, match="mirror_pos"):
            Permeability(index=0, dim="x", mirror_pos=1.01)
        # Test that we do pss if not relative
        _ = Permeability(index=0, dim="x", mirror_pos=1.01, relative=False)

    def test_x_calculation(self):
        # Test that this is just the position
        assert self.op.calculate(self.system) == [self.x, 0, 1]

    def test_x_wrap(self):
        # set x to a value, that +ofsett will be wrapped
        op1 = Permeability(index=0, dim="x")
        op2 = Permeability(index=0, dim="x", offset=0.5)
        x = 0.8
        self.system.particles.pos[0, 0] = x
        assert op1.calculate(self.system)[0] == x
        assert op2.calculate(self.system)[0] == pytest.approx(0.3)

    def test_box_min_under_0(self):
        # set new box
        self.system.box.length[0] = 2
        self.system.box.low[0] = -1
        assert self.op.calculate(self.system)[0] == pytest.approx(self.x)
        op2 = Permeability(index=0, dim="x", relative=True)
        # x should be 0.5 + x/2 (adding a box lenght to the left)
        assert op2.calculate(self.system)[0] == pytest.approx(0.5 + self.x / 2)

    def test_broken_box(self):
        # Break box
        self.system.box = None
        # Make an OP that would normally break the box

        op2 = Permeability(index=0, dim="x", relative=False, offset=-12)
        self.op.periodic = False
        op2.periodic = False

        assert self.op.calculate(self.system)[0] == pytest.approx(self.x)
        assert op2.calculate(self.system)[0] == pytest.approx(self.x - 12)

    def test_mirrored_function(self):
        assert not self.op._mirror
        assert self.op.calculate(self.system)[0] == pytest.approx(self.x)
        assert self.op.calculate(self.system)[2] == pytest.approx(1)
        self.op.mirror()
        assert self.op._mirror
        assert self.op.calculate(self.system)[0] == pytest.approx(1 - self.x)
        assert self.op.calculate(self.system)[2] == pytest.approx(-1)

    def test_mirror_composite(self):
        orderp = CompositeOrderParameter()
        assert not self.op._mirror
        orderp.add_orderparameter(self.op)
        op2 = Permeability(index=0, dim="x", relative=True)
        op2.mirror()
        assert op2._mirror
        orderp.add_orderparameter(op2)
        orderp.mirror()
        assert self.op._mirror
        assert not op2._mirror

    def test_mirror_composite_warning(self, caplog):
        orderp = CompositeOrderParameter()
        op1 = Position(index=0)
        orderp.add_orderparameter(op1)
        assert not self.op._mirror
        orderp.add_orderparameter(self.op)
        op2 = Permeability(index=0, dim="x", relative=True)
        op2.mirror()
        assert op2._mirror
        orderp.add_orderparameter(op2)
        ln = "pyretis.orderparameter.orderparameter"
        logging.disable(logging.INFO)

        with caplog.at_level(logging.INFO, logger=ln):
            orderp.mirror()
        logging.disable(logging.CRITICAL)
        assert "Attempting a mirror move, but" in caplog.text
        assert self.op._mirror
        assert not op2._mirror

    def test_index_get_set_composite(self):
        orderp = CompositeOrderParameter()
        op1 = Position(index=1)
        orderp.add_orderparameter(op1)
        orderp.add_orderparameter(self.op)
        assert orderp.index == 1
        orderp.index = 12
        assert op1.index == 12
        assert self.op.index == 0

    def test_permeabilityminusoffset(self):
        # Tak e offset that is bigger than the box (normally be wrapped)
        offset = 10
        opmin = PermeabilityMinusOffset(index=0, dim="x", relative=False)
        # Check that with an offset of 0, this is equal
        assert self.op.calculate(self.system) == opmin.calculate(self.system)
        self.op.offset = offset
        opmin.offset = offset
        assert (
            self.op.calculate(self.system)[0]
            == opmin.calculate(self.system)[0] + offset
        )

    def test_restart_cycle(self):
        # See if all the info is there
        info = self.op.restart_info()
        assert info == {"index": 0, "mirror": False}
        # See if mirror is given properly
        self.op.mirror()
        info = self.op.restart_info()
        assert info == {"index": 0, "mirror": True}
        # See if index is given properly
        self.op.index = "bla"
        info = self.op.restart_info()
        assert info == {"index": "bla", "mirror": True}
        # See if the restart is done properly
        info = {"index": "blub", "mirror": "Fish"}
        self.op.load_restart_info(info)
        assert self.op.index == "blub"
        assert self.op._mirror == "Fish"

    def test_restart_cycle_composite(self):
        orderp = CompositeOrderParameter()
        op1 = Position(index=1)
        orderp.add_orderparameter(op1)
        orderp.add_orderparameter(self.op)
        info = orderp.restart_info()
        ref_info = [None, {"index": 0, "mirror": False}]
        assert info == ref_info
        # See if we can set properly as well
        info = [None, {"index": 42, "mirror": "towel"}]
        orderp.load_restart_info(info)
        assert self.op.index == 42
        assert self.op._mirror == "towel"


# ---------------------------------------------------------------------------
# Tests for normalize_order_output and wrap_orderparameter
# ---------------------------------------------------------------------------


class TestNormalizeOrderOutput:
    """Test normalize_order_output for all supported input types."""

    def test_python_int(self):
        assert normalize_order_output(3) == [3]

    def test_python_float(self):
        assert normalize_order_output(1.5) == [1.5]

    def test_python_list_single(self):
        assert normalize_order_output([0.5]) == [0.5]

    def test_python_list_multi(self):
        assert normalize_order_output([0.1, 0.2]) == [0.1, 0.2]

    def test_python_tuple(self):
        assert normalize_order_output((0.1, 0.2)) == [0.1, 0.2]

    def test_numpy_scalar(self):
        val = np.float64(2.5)
        result = normalize_order_output(val)
        assert result == [2.5]
        assert isinstance(result, list)

    def test_numpy_0d_array(self):
        arr = np.array(3.14)
        result = normalize_order_output(arr)
        assert result == pytest.approx([3.14])
        assert isinstance(result, list)

    def test_numpy_1d_array(self):
        arr = np.array([0.1, 0.2, 0.3])
        result = normalize_order_output(arr)
        assert result == pytest.approx([0.1, 0.2, 0.3])
        assert isinstance(result, list)

    def test_numpy_2d_array_raises(self):
        arr = np.array([[0.1, 0.2], [0.3, 0.4]])
        with pytest.raises(ValueError):
            normalize_order_output(arr)

    def test_nested_list_raises(self):
        with pytest.raises(ValueError):
            normalize_order_output([[0.1, 0.2]])

    def test_non_iterable_raises(self):
        with pytest.raises((ValueError, TypeError)):
            normalize_order_output(object())


class _ScalarOP(OrderParameter):
    """Helper: returns a bare float."""

    def __init__(self):
        super().__init__(description="scalar")

    def calculate(self, system):
        return system.temperature["set"]


class _ArrayOP(OrderParameter):
    """Helper: returns a 1-D numpy array."""

    def __init__(self):
        super().__init__(description="array")

    def calculate(self, system):
        return np.array([system.temperature["set"], 0.0])


class _AlreadyWrapped(OrderParameter):
    """Helper: returns a proper list already."""

    def __init__(self):
        super().__init__(description="list")

    def calculate(self, system):
        return [system.temperature["set"]]


class TestWrapOrderparameter:
    """Test wrap_orderparameter for all return-type variants."""

    def setup_method(self):
        self.system = System(temperature=42.0, units="lj", box=None)

    def test_scalar_normalised_to_list(self):
        op = wrap_orderparameter(_ScalarOP())
        result = op.calculate(self.system)
        assert isinstance(result, list)
        assert result == [42.0]

    def test_numpy_array_normalised_to_list(self):
        op = wrap_orderparameter(_ArrayOP())
        result = op.calculate(self.system)
        assert isinstance(result, list)
        assert result == pytest.approx([42.0, 0.0])

    def test_list_return_unchanged(self):
        op = wrap_orderparameter(_AlreadyWrapped())
        result = op.calculate(self.system)
        assert isinstance(result, list)
        assert result == [42.0]

    def test_idempotent(self):
        """Wrapping twice must not double-wrap."""
        op = _ScalarOP()
        op1 = wrap_orderparameter(op)
        op2 = wrap_orderparameter(op)
        assert op1 is op2
        result = op2.calculate(self.system)
        assert result == [42.0]

    def test_composite_with_scalar_subop(self):
        """CompositeOrderParameter handles sub-OPs that return scalars."""
        sub = wrap_orderparameter(_ScalarOP())
        composite = CompositeOrderParameter(order_parameters=[sub, sub])
        result = composite.calculate(self.system)
        assert isinstance(result, list)
        assert result == [42.0, 42.0]


class TestExpandOrderNames:
    """Tests for the expand_order_names helper."""

    def test_default_single(self):
        """A single value with no name uses the indexed default."""
        assert expand_order_names(None, 1) == ["op_1"]

    def test_default_multi(self):
        """Default labels for a multi-element OP."""
        assert expand_order_names(None, 5) == [
            "op_1", "op_2", "op_3", "op_4", "op_5",
        ]

    def test_default_cv_prefix(self):
        """The default prefix can be customised (used for CVs)."""
        assert expand_order_names(None, 3, default_prefix="cv") == [
            "cv_1", "cv_2", "cv_3",
        ]

    def test_string_single(self):
        """A single value with a string name uses the bare name."""
        assert expand_order_names("pippo", 1) == ["pippo"]

    def test_string_multi(self):
        """A string name expands with index suffixes for multi values."""
        assert expand_order_names("pippo", 3) == [
            "pippo_1", "pippo_2", "pippo_3",
        ]

    def test_list_match(self):
        """A list of names is used as-is when its length matches."""
        labels = ["a", "b", "c"]
        assert expand_order_names(labels, 3) == labels

    def test_list_mismatch_raises(self):
        """A list whose length differs from count raises ValueError."""
        with pytest.raises(ValueError, match="does not match"):
            expand_order_names(["a", "b"], 3)

    def test_list_with_non_string_raises(self):
        """A list containing non-strings raises ValueError."""
        with pytest.raises(ValueError, match="strings"):
            expand_order_names(["a", 2, "c"], 3)

    def test_invalid_type_raises(self):
        """An unsupported name type raises ValueError."""
        with pytest.raises(ValueError, match="must be a string"):
            expand_order_names(42, 1)

    def test_invalid_count_raises(self):
        """A non-positive count raises ValueError."""
        with pytest.raises(ValueError, match="positive integer"):
            expand_order_names(None, 0)

    def test_start_index(self):
        """start_index shifts the indexed labels."""
        assert expand_order_names(None, 2, start_index=4) == [
            "op_4", "op_5",
        ]
        assert expand_order_names("foo", 2, start_index=3) == [
            "foo_3", "foo_4",
        ]
