from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.nuclear.buckling import geometric_buckling_for_uniform_cylinder as buckling


@fixture(name="test_args")
def test_args_fixture():
    cylinder_radius = Quantity(1.8 * units.meter)
    cylinder_height = Quantity(4.2 * units.meter)
    Args = namedtuple("Args", ["R", "H"])
    return Args(R=cylinder_radius, H=cylinder_height)


def test_basic_geometric_buckling(test_args):
    result = buckling.calculate_geometric_buckling_squared(test_args.R, test_args.H)
    assert_equal(result, 2.3447 / units.meter**2)


def test_bad_cylinder_radius(test_args):
    Rb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        buckling.calculate_geometric_buckling_squared(Rb, test_args.H)
    with raises(TypeError):
        buckling.calculate_geometric_buckling_squared(100, test_args.H)


def test_bad_cylinder_height(test_args):
    Hb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        buckling.calculate_geometric_buckling_squared(test_args.R, Hb)
    with raises(TypeError):
        buckling.calculate_geometric_buckling_squared(test_args.R, 100)
