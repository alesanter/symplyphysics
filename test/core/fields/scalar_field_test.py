from collections import namedtuple
from pytest import fixture, raises
from sympy import atan, cos, pi, sin, sqrt, symbols, simplify, Expr
from sympy.vector import express
from symplyphysics.core.test_decorators import unsupported_usage
from symplyphysics.core.points.cylinder_point import CylinderPoint
from symplyphysics.core.points.sphere_point import SpherePoint
from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem, coordinates_rotate, coordinates_transform
from symplyphysics.core.points.cartesian_point import CartesianPoint
from symplyphysics.core.points.point import Point
from symplyphysics.core.fields.scalar_field import ScalarField

Args = namedtuple("Args", ["C"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    C = CoordinateSystem()
    return Args(C=C)


# Test ScalarField constructor


def test_basic_field() -> None:

    def field_function(p: CartesianPoint) -> Expr:
        return p.z * p.y

    field = ScalarField(field_function)
    field_point = CartesianPoint(1, 2, 3)
    assert field(field_point) == 6
    assert len(field.basis) == 3
    assert callable(field.field_function)


def test_empty_field() -> None:
    field = ScalarField()
    field_point = Point(1, 2, 3)
    assert field(field_point) == 0
    # Scalar field components size is always 1
    assert field.field_function == 0


def test_4d_point_field() -> None:
    field = ScalarField(lambda p: p.coordinate(3))
    field_point = Point(1, 2, 3, 4)
    assert field(field_point) == 4


@unsupported_usage
def test_invalid_lambda_field() -> None:
    field = ScalarField(lambda p: p.y + "string")
    field_point = CartesianPoint(1, 2, 3)
    # cannot add integer and string in field lambda
    with raises(TypeError):
        field(field_point)


def test_coord_system_field(test_args: Args) -> None:
    field = ScalarField(lambda p: p.y * p.z, test_args.C)
    field_point = CartesianPoint(1, 2, 3)
    assert field(field_point) == 6
    assert field.basis == [
        test_args.C.coord_system.x, test_args.C.coord_system.y, test_args.C.coord_system.z
    ]
    assert field.coordinate_system == test_args.C


def test_wrong_coord_system_field(test_args: Args) -> None:
    field = ScalarField(lambda p: p.y * p.z, test_args.C)
    field_point_sphere = SpherePoint(1, 2, 3)
    with raises(ValueError):
        field(field_point_sphere)
    field_point_cylinder = CylinderPoint(1, 2, 3)
    with raises(ValueError):
        field(field_point_cylinder)


def test_wrong_lambda_field() -> None:
    field = ScalarField(lambda p: p.y * p.z)
    field_point = Point(1, 2, 3)
    with raises(AttributeError):
        field(field_point)


# Test ScalarField.from_expression()


def test_basic_vector_to_field_conversion(test_args: Args) -> None:
    field = ScalarField.from_expression(test_args.C.coord_system.x + test_args.C.coord_system.y,
        test_args.C)
    field_point = Point(1, 2, 3)
    assert field(field_point) == 3
    assert field.basis == [
        test_args.C.coord_system.x, test_args.C.coord_system.y, test_args.C.coord_system.z
    ]
    assert field.coordinate_system == test_args.C


# Coordinate system base vectors (C.i, C.j) are not being processed by ScalarField,
# so as a result of applying field they are kept untouched, like all other free variables in SymPy expression.
def test_dimensional_vector_to_field_conversion(test_args: Args) -> None:
    field = ScalarField.from_expression(
        test_args.C.coord_system.x * test_args.C.coord_system.i +
        test_args.C.coord_system.y * test_args.C.coord_system.j, test_args.C)
    field_point = Point(1, 2, 3)
    assert field(field_point) == test_args.C.coord_system.i + 2 * test_args.C.coord_system.j


def test_empty_vector_to_field_conversion() -> None:
    field = ScalarField.from_expression(0)
    # applying empty field to a point results in zero value
    field_point = Point(1, 2, 3)
    assert field(field_point) == 0
    assert len(field.basis) == 3


def test_only_integer_vector_to_field_conversion() -> None:
    field = ScalarField.from_expression(1)
    field_point = Point(1, 2, 3)
    assert field(field_point) == 1
    assert len(field.basis) == 3


# Only scalars from requested coordinate system are being applied
def test_external_symbols_vector_to_field(test_args: Args) -> None:
    phi = symbols("phi")
    field = ScalarField.from_expression(test_args.C.coord_system.x + 2 * phi, test_args.C)
    assert callable(field.field_function)
    field_point = Point(1, 2, 3)
    assert field(field_point) == 1 + 2 * phi


def test_custom_names_vector_to_field_conversion() -> None:
    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    r, theta, z = C1.coord_system.base_scalars()
    field = ScalarField.from_expression(r + 2 * theta, C1)
    field_point = CylinderPoint(1, 2, 3)
    assert field(field_point) == 5
    assert field.basis == [r, theta, z]
    assert field.coordinate_system == C1


def test_rotate_coordinates_vector_to_field_conversion(test_args: Args) -> None:
    sympy_vector_field = test_args.C.coord_system.x + test_args.C.coord_system.y
    field = ScalarField.from_expression(sympy_vector_field, test_args.C)
    field_point = Point(1, 2, 3)
    assert field(field_point) == 3
    theta = symbols("theta")
    B = coordinates_rotate(test_args.C, theta, test_args.C.coord_system.k)
    x, y, z = B.coord_system.base_scalars()
    transformed_vector = express(sympy_vector_field, B.coord_system, variables=True)
    result_transformed_field = ScalarField.from_expression(transformed_vector, B)
    assert transformed_vector == x * sin(theta) + x * cos(theta) - y * sin(theta) + y * cos(theta)
    assert result_transformed_field(
        field_point) == sin(theta) + cos(theta) - 2 * sin(theta) + 2 * cos(theta)
    assert result_transformed_field.basis == [x, y, z]
    assert result_transformed_field.coordinate_system == B


# when we express SymPy Vector to another coordinate system and base scalars (C.x, C.y) are
# not transformed, we have same base scalars as before.
def test_rotate_coordinates_without_variables_vector_to_field_conversion(test_args: Args) -> None:
    theta = symbols("theta")
    B = coordinates_rotate(test_args.C, theta, test_args.C.coord_system.k)
    transformed_field = express(test_args.C.coord_system.x + test_args.C.coord_system.y,
        B.coord_system)
    field = ScalarField.from_expression(transformed_field, B)
    field_point = Point(1, 2, 3)
    assert field(field_point) == test_args.C.coord_system.x + test_args.C.coord_system.y


# Test field.apply()


# Result is a function that returns a scalar value at any point of the trajectory.
# Input field has X * Y, so as we are moving along X-axis or Y-axis of the trajectory,
# resulting value grows proportionally to X * Y.
def test_basic_field_apply(test_args: Args) -> None:
    field = ScalarField(lambda p: p.y * p.x)
    # represents surface
    trajectory = [test_args.C.coord_system.x, test_args.C.coord_system.y]
    trajectory_value = field.apply(trajectory)
    assert trajectory_value == test_args.C.coord_system.y * test_args.C.coord_system.x


# Coordinate system is not necessary to apply field.
def test_parametrized_no_coord_system_field_apply() -> None:
    field = ScalarField(lambda p: -p.y)
    parameter = symbols("parameter")
    # represents y = x trajectory
    trajectory = [parameter, parameter]
    trajectory_value = field.apply(trajectory)
    assert trajectory_value == -1 * parameter


def test_parametrized_field_apply(test_args: Args) -> None:
    result_field = ScalarField.from_expression(
        -test_args.C.coord_system.y + 2 * test_args.C.coord_system.x, test_args.C)
    parameter = symbols("parameter")
    # represents y = x trajectory
    trajectory = [parameter, parameter]
    trajectory_value = result_field.apply(trajectory)
    assert trajectory_value == parameter


def test_sympy_field_apply(test_args: Args) -> None:
    result_field = ScalarField.from_expression(
        -test_args.C.coord_system.y + test_args.C.coord_system.x, test_args.C)
    trajectory = [test_args.C.coord_system.x, test_args.C.coord_system.y]
    trajectory_value = result_field.apply(trajectory)
    assert trajectory_value == -test_args.C.coord_system.y + test_args.C.coord_system.x


def test_uncallable_field_apply(test_args: Args) -> None:
    result_field = ScalarField(1)
    trajectory = [test_args.C.coord_system.x, test_args.C.coord_system.y]
    trajectory_value = result_field.apply(trajectory)
    assert trajectory_value == 1


# ScalarField is not rebased automatically and should be rebased to the same coordinate
# system as in trajectory with 'field_rebase'.
def test_different_coord_systems_field_apply(test_args: Args) -> None:
    result_field = ScalarField(lambda p: p.y * p.x, test_args.C)
    C1 = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    r, theta, _ = C1.coord_system.base_scalars()
    trajectory = [r, theta]
    trajectory_value = result_field.apply(trajectory)
    assert trajectory_value == theta * r


# While ScalarField should contain information about coordinate system and
# can be rebased to new coordinate system with ScalarField invariance, it is not
# necessary to do so. Instead user can define a trajectory, that contains this point, transform
# trajectory to new coordinate system and apply field to it.
# ScalarField invariant will hold.
def test_invariant_transformed_trajectory_field_apply(test_args: Args) -> None:
    field = ScalarField(lambda p: p.x**2 + 2 * p.y**2)
    point = [1, 2]
    trajectory = [test_args.C.coord_system.x, test_args.C.coord_system.y + 5]
    trajectory_value = field.apply(trajectory)
    assert trajectory_value == test_args.C.coord_system.x**2 + 2 * (test_args.C.coord_system.y +
        5)**2
    assert trajectory_value.subs({
        test_args.C.coord_system.x: point[0],
        test_args.C.coord_system.y: point[1]
    }) == 99

    B = coordinates_rotate(test_args.C, pi / 4, test_args.C.coord_system.k)
    x, y, _ = B.coord_system.base_scalars()
    transformed_trajectory = [
        express(trajectory[0], B.coord_system, variables=True),
        express(trajectory[1], B.coord_system, variables=True)
    ]
    transformed_trajectory_value = field.apply(transformed_trajectory)

    p1 = test_args.C.coord_system.origin.locate_new(
        "p1", point[0] * test_args.C.coord_system.i + point[1] * test_args.C.coord_system.j)
    p1_coordinates = p1.express_coordinates(test_args.C.coord_system)
    assert p1_coordinates[0] == point[0]
    assert p1_coordinates[1] == point[1]

    p1_coordinates_in_b = p1.express_coordinates(B.coord_system)
    assert p1_coordinates_in_b[0] != point[0]

    assert transformed_trajectory_value.subs({
        x: p1_coordinates_in_b[0],
        y: p1_coordinates_in_b[1]
    }) == 99


# Test field.apply_to_basis()


def test_basic_field_apply_to_basis(test_args: Args) -> None:
    field = ScalarField(lambda p: p.y * p.x, test_args.C)
    field_space = field.apply_to_basis()
    assert field_space == test_args.C.coord_system.y * test_args.C.coord_system.x


# Test field_rebase()


def test_basic_field_rebase(test_args: Args) -> None:
    field = ScalarField(lambda p: p.x + p.y, test_args.C)
    point = [1, 2]
    point_value = field.apply(point)
    assert point_value == 3
    assert field.coordinate_system == test_args.C

    # B is located at [1, 2] origin instead of [0, 0] of test_args.C
    Bi = test_args.C.coord_system.locate_new(
        "B", test_args.C.coord_system.i + 2 * test_args.C.coord_system.j)
    B = CoordinateSystem(test_args.C.coord_system_type, Bi)
    field_rebased = field.rebase(B)
    x, y, z = B.coord_system.base_scalars()
    assert field_rebased.basis == [x, y, z]
    assert field_rebased.coordinate_system == B
    # Original field is not changed
    assert field.basis == [
        test_args.C.coord_system.x, test_args.C.coord_system.y, test_args.C.coord_system.z
    ]
    assert field.coordinate_system == test_args.C

    transformed_point_value = field_rebased.apply(point)
    assert transformed_point_value != point_value
    assert transformed_point_value == 6


# ScalarField invariant does not hold, when applied to some fixed point in space. Use
# 'field_rebase' to let ScalarField know about new coordinate system.
def test_invariant_field_rebase_and_apply(test_args: Args) -> None:
    field = ScalarField(lambda p: p.x**2 + 2 * p.y**2, test_args.C)
    point = [1, 2]
    p1 = test_args.C.coord_system.origin.locate_new(
        "p1", point[0] * test_args.C.coord_system.i + point[1] * test_args.C.coord_system.j)
    p1_coordinates = p1.express_coordinates(test_args.C.coord_system)
    assert p1_coordinates[0] == point[0]
    assert p1_coordinates[1] == point[1]

    point_value = field.apply(point)
    assert point_value == 9

    B = coordinates_rotate(test_args.C, pi / 4, test_args.C.coord_system.k)
    p1_coordinates_in_b = p1.express_coordinates(B.coord_system)
    assert p1_coordinates_in_b[0] != point[0]

    transformed_point = [p1_coordinates_in_b[0], p1_coordinates_in_b[1]]
    transformed_point_value = field.apply(transformed_point)
    # invariant does not hold if field is not rebased to new coordinate system
    assert transformed_point_value != point_value

    field_rebased = field.rebase(B)
    transformed_point_value = field_rebased.apply(transformed_point)
    assert transformed_point_value == point_value


# Field is not rebased if no original coordinate system was set.
def test_no_coord_system_field_rebase(test_args: Args) -> None:
    field = ScalarField(lambda p: p.x + p.y)
    point = [1, 2]
    point_value = field.apply(point)
    assert point_value == 3
    Bi = test_args.C.coord_system.locate_new(
        "B", test_args.C.coord_system.i + 2 * test_args.C.coord_system.j)
    B = CoordinateSystem(test_args.C.coord_system_type, Bi)
    with raises(ValueError):
        field.rebase(B)


# Test non-cartesian coordinate systems


def test_cylindrical_field_create(test_args: Args) -> None:
    field = ScalarField(lambda p: p.x + p.y, test_args.C)
    point = [1, 2]
    point_value = field.apply(point)
    assert point_value == 3

    B = coordinates_transform(test_args.C, CoordinateSystem.System.CYLINDRICAL)
    field_rebased = field.rebase(B)
    assert field_rebased.coordinate_system == B

    # point should have r = sqrt(5) in polar coordinates
    # theta angle is atan(2/1)
    point_polar = [sqrt(5), atan(2)]
    point_polar_value = field_rebased.apply(point_polar)
    assert simplify(point_polar_value) == 3
