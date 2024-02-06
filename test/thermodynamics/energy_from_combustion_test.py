from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
    prefixes,
)
from symplyphysics.laws.thermodynamics import energy_from_combustion as amount_energy

# Test example from https://easyfizika.ru/zadachi/termodinamika/skolko-tepla-vydelitsya-pri-sgoranii-2-kg-benzina/
# When burning 2 kg of gasoline, 92 MJ energy should be released


@fixture(name="test_args")
def test_args_fixture():
    k_q = Quantity(46 * prefixes.mega * units.joule / units.kilogram)
    m = Quantity(2 * units.kilogram)
    Args = namedtuple("Args", ["k_q", "m"])
    return Args(k_q=k_q, m=m)


def test_basic_amount(test_args):
    result = amount_energy.calculate_amount_energy(test_args.k_q, test_args.m)
    assert_equal(result, 92 * prefixes.mega * units.joules)


def test_bad_specific_heat_combustion(test_args):
    k_qb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        amount_energy.calculate_amount_energy(k_qb, test_args.m)
    with raises(TypeError):
        amount_energy.calculate_amount_energy(100, test_args.m)


def test_bad_body_mass(test_args):
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        amount_energy.calculate_amount_energy(test_args.k_q, mb)
    with raises(TypeError):
        amount_energy.calculate_amount_energy(test_args.k_q, 100)
