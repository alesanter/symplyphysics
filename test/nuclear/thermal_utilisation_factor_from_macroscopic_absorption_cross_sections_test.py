from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.core.symbols.probability import Probability
from symplyphysics.laws.nuclear import thermal_utilisation_factor_from_macroscopic_absorption_cross_sections as utilisation_factor

Args = namedtuple("Args", ["Saf", "Sat"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    macro_abs_fuel_cross_section = Quantity(0.2028 / units.centimeter)
    macro_abs_total_cross_section = Quantity(0.2356 / units.centimeter)
    return Args(Saf=macro_abs_fuel_cross_section, Sat=macro_abs_total_cross_section)


def test_basic_utilisation_factor(test_args: Args) -> None:
    result = utilisation_factor.calculate_utilisation_factor(test_args.Saf, test_args.Sat)
    assert isinstance(result, Probability)
    assert_equal(result, 0.861)


def test_bad_macroscopic_cross_section(test_args: Args) -> None:
    Sb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        utilisation_factor.calculate_utilisation_factor(Sb, test_args.Sat)
    with raises(TypeError):
        utilisation_factor.calculate_utilisation_factor(100, test_args.Sat)
    with raises(errors.UnitsError):
        utilisation_factor.calculate_utilisation_factor(test_args.Saf, Sb)
    with raises(TypeError):
        utilisation_factor.calculate_utilisation_factor(test_args.Saf, 100)

    Sb = Quantity(test_args.Sat.scale_factor + 1)
    with raises(TypeError):
        utilisation_factor.calculate_utilisation_factor(Sb, test_args.Sat)
