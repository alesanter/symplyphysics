from collections import namedtuple
from pytest import fixture, raises
from sympy import Eq
from symplyphysics import assert_equal_vectors, units, errors, Quantity, QuantityVector
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.experimental.vectors import VectorSymbol, VectorNorm as norm
from symplyphysics.core.experimental.solvers import solve_for_scalar, solve_for_vector, vector_equals, apply
from symplyphysics.laws.dynamics.vector import acceleration_from_force as newton_second_law

Args = namedtuple("Args", ["m", "a", "f"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(1 * units.kilogram)
    a = QuantityVector([Quantity(3 * units.meter / units.second**2)])
    f = QuantityVector([Quantity(3 * units.newton)])
    return Args(m=m, a=a, f=f)


def test_basic_force(test_args: Args) -> None:
    result = newton_second_law.calculate_force(test_args.m, test_args.a)
    assert_equal_vectors(
        result,
        QuantityVector([3 * units.newton]),
    )


def test_basic_acceleration(test_args: Args) -> None:
    result = newton_second_law.calculate_acceleration(test_args.m, test_args.f)
    assert_equal_vectors(
        result,
        QuantityVector([3 * units.meter / units.second**2]),
    )


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        newton_second_law.calculate_force(mb, test_args.a)
    with raises(TypeError):
        newton_second_law.calculate_force(100, test_args.a)
    with raises(errors.UnitsError):
        newton_second_law.calculate_acceleration(mb, test_args.f)
    with raises(TypeError):
        newton_second_law.calculate_acceleration(100, test_args.f)


def test_bad_acceleration(test_args: Args) -> None:
    ab = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        newton_second_law.calculate_force(test_args.m, ab)
    with raises(TypeError):
        newton_second_law.calculate_force(test_args.m, 100)


def test_bad_force(test_args: Args) -> None:
    fb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        newton_second_law.calculate_acceleration(test_args.m, fb)
    with raises(TypeError):
        newton_second_law.calculate_acceleration(test_args.m, 100)


# TODO: Remove once the law is fixed to account for the new vector functionality.
def test_vector_symbols() -> None:
    force = VectorSymbol("F", units.force)
    acceleration = VectorSymbol("a", units.acceleration)
    mass = newton_second_law.mass

    force_law = Eq(force, acceleration * mass)

    acceleration_law = solve_for_vector(force_law, acceleration)
    assert vector_equals(acceleration_law.lhs, acceleration)
    assert vector_equals(acceleration_law.rhs, force / mass)

    mass_eqn = apply(force_law, norm)
    mass_law = solve_for_scalar(mass_eqn, mass)[0]
    assert expr_equals(mass_law.lhs, mass)
    assert expr_equals(mass_law.rhs, norm(force) / norm(acceleration))
