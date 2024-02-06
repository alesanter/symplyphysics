from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity import current_is_proportional_to_voltage as ohms_law

# Description
## Assert we have 3 Volts applied to 2-Ohm resistor.
## According to Ohm's Law we should have 3/2 = 1.5Amps current flowing through this resistor.


@fixture(name="test_args")
def test_args_fixture():
    V = Quantity(3 * units.volt)
    R = Quantity(2 * units.ohm)
    Args = namedtuple("Args", ["V", "R"])
    return Args(V=V, R=R)


def test_basic_current(test_args):
    result = ohms_law.calculate_current(test_args.V, test_args.R)
    assert_equal(result, 1.5 * units.ampere)


def test_bad_voltage(test_args):
    Vb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        ohms_law.calculate_current(Vb, test_args.R)
    with raises(TypeError):
        ohms_law.calculate_current(100, test_args.R)


def test_bad_resistance(test_args):
    Rb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        ohms_law.calculate_current(test_args.V, Rb)
    with raises(TypeError):
        ohms_law.calculate_current(test_args.V, 100)
