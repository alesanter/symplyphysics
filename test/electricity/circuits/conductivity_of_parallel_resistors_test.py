from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity.circuits import conductivity_of_parallel_resistors as parallel_resistor

# Description
## Assert we have two resistors with 1/2 Siemens and 1/4 Siemens conductance.
## Accordind to calculator (https://www.chipdip.ru/calc/parallel-resistors) resulting conductance should be 0.75 Siemens.


@fixture(name="test_args")
def test_args_fixture():
    S1 = Quantity(1 / 2 * units.siemens)
    S2 = Quantity(1 / 4 * units.siemens)
    Args = namedtuple("Args", ["S1", "S2"])
    return Args(S1=S1, S2=S2)


def test_basic_conductivity(test_args):
    result = parallel_resistor.calculate_parallel_conductance([test_args.S1, test_args.S2])
    assert_equal(result, 0.75 * units.siemens)


def test_three_resistors_array(test_args):
    S3 = Quantity(1 * units.siemens)
    result = parallel_resistor.calculate_parallel_conductance([test_args.S1, test_args.S2, S3])
    assert_equal(result, 1.75 * units.siemens)


def test_bad_conductivity(test_args):
    Sb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        parallel_resistor.calculate_parallel_conductance([Sb, test_args.S2])
    with raises(TypeError):
        parallel_resistor.calculate_parallel_conductance([100, test_args.S2])
    with raises(errors.UnitsError):
        parallel_resistor.calculate_parallel_conductance([test_args.S1, Sb])
    with raises(TypeError):
        parallel_resistor.calculate_parallel_conductance([test_args.S1, 100])
    with raises(errors.UnitsError):
        parallel_resistor.calculate_parallel_conductance([Sb, Sb])
    with raises(TypeError):
        parallel_resistor.calculate_parallel_conductance([100, 100])
    with raises(TypeError):
        parallel_resistor.calculate_parallel_conductance(test_args.S1)
