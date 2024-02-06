from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, errors, units, Quantity, prefixes)
from symplyphysics.laws.electricity import energy_accumulated_in_inductor_from_inductance_and_current as inductor_law

# Description
## Assert we have 150mH inductor with 0.5A current flowing through it.
## According to law we should have amount of energy accumulated in this inductor equals to 0.150 * 0.5**2 / 2 = 0.01875 Joules.


@fixture(name="test_args")
def test_args_fixture():
    I = Quantity(150 * prefixes.milli * units.henry)
    C = Quantity(0.5 * units.ampere)
    Args = namedtuple("Args", ["I", "C"])
    return Args(I=I, C=C)


def test_basic_energy(test_args):
    result = inductor_law.calculate_accumulated_energy(test_args.I, test_args.C)
    assert_equal(result, 0.01875 * units.joule)


def test_bad_inductance(test_args):
    Ib = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        inductor_law.calculate_accumulated_energy(Ib, test_args.C)
    with raises(TypeError):
        inductor_law.calculate_accumulated_energy(100, test_args.C)


def test_bad_current(test_args):
    Cb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        inductor_law.calculate_accumulated_energy(test_args.I, Cb)
    with raises(TypeError):
        inductor_law.calculate_accumulated_energy(test_args.I, 100)
