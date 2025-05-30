from pytest import raises
from sympy import Derivative, cos, pi, Symbol as SymSymbol
from symplyphysics import (units, Quantity, SI, dimensionless)
from symplyphysics.core.symbols.quantities import scale_factor

# Test Quantity constructor


def test_basic_quantity() -> None:
    a = Quantity(10 * units.meter / units.second)
    expr = a
    q = Quantity(expr)
    assert q.scale_factor == 10
    assert SI.get_dimension_system().equivalent_dims(q.dimension, units.speed)


def test_prefixed_quantity() -> None:
    a = Quantity(10 * units.kilometer)
    expr = a
    q = Quantity(expr)
    assert q.scale_factor == 10000
    assert SI.get_dimension_system().equivalent_dims(q.dimension, units.length)


def test_dimensionless_quantity() -> None:
    a = Quantity(10)
    expr = a
    q = Quantity(expr)
    assert q.scale_factor == 10
    assert SI.get_dimension_system().equivalent_dims(q.dimension, dimensionless)


def test_sum_quantity() -> None:
    a = Quantity(10 * units.meter / units.second)
    b = Quantity(5 * units.meter / units.second)
    expr = a + b
    q = Quantity(expr)
    assert q.scale_factor == 15
    assert SI.get_dimension_system().equivalent_dims(q.dimension, units.speed)

    expr = b + a
    q = Quantity(expr)
    assert q.scale_factor == 15
    assert SI.get_dimension_system().equivalent_dims(q.dimension, units.speed)


def test_invalid_sum_quantity() -> None:
    a = Quantity(10 * units.meter / units.second)
    b = Quantity(5 * units.meter)
    expr = a + b
    with raises(ValueError):
        Quantity(expr)
    expr = b + a
    with raises(ValueError):
        Quantity(expr)


def test_mul_quantity() -> None:
    a = Quantity(10 * units.meter / units.second)
    b = Quantity(5 * units.second)
    expr = a * b
    q = Quantity(expr)
    assert q.scale_factor == 50
    assert SI.get_dimension_system().equivalent_dims(q.dimension, units.length)

    expr = b * a
    q = Quantity(expr)
    assert q.scale_factor == 50
    assert SI.get_dimension_system().equivalent_dims(q.dimension, units.length)


def test_pow_quantity() -> None:
    a = Quantity(10 * units.meter / units.second)
    b = Quantity(2)
    expr = a**b
    q = Quantity(expr)
    assert q.scale_factor == 100
    assert SI.get_dimension_system().equivalent_dims(q.dimension, units.speed**2)


def test_invalid_pow_quantity() -> None:
    a = Quantity(10 * units.meter / units.second)
    b = Quantity(2)
    expr = b**a
    with raises(ValueError):
        Quantity(expr)


def test_invalid_derivative_quantity() -> None:
    a = Quantity(10 * units.second)
    expr = Derivative(a**2, a)
    with raises(ValueError):
        Quantity(expr)


def test_function_quantity() -> None:
    a = Quantity(pi)
    expr = cos(a)
    # Make sure it is unevaluated
    assert expr != -1
    q = Quantity(expr)
    assert q.scale_factor == -1
    assert SI.get_dimension_system().equivalent_dims(q.dimension, dimensionless)


def test_invalid_function_quantity() -> None:
    a = Quantity(pi * units.second)
    expr = cos(a)
    # Make sure it is unevaluated
    assert expr != -1
    with raises(ValueError):
        Quantity(expr)


def test_dimensionless_expr_function_quantity() -> None:
    a = Quantity(pi * units.second)
    b = Quantity(1 * units.second)
    expr = cos(a / b)
    # Make sure it is unevaluated
    assert expr != -1
    q = Quantity(expr)
    assert q.scale_factor == -1
    assert SI.get_dimension_system().equivalent_dims(q.dimension, dimensionless)


def test_scale_factor() -> None:
    a_quantity = Quantity(1.0 * units.ampere)
    a_float = 1.0
    assert scale_factor(a_quantity) == 1.0
    assert scale_factor(a_float) == 1.0


def test_invalid_free_symbol_quantity() -> None:
    a = Quantity(10 * units.meter / units.second)
    b = SymSymbol("b")
    expr = a * b
    with raises(ValueError):
        Quantity(expr)
    expr = b * a
    with raises(ValueError):
        Quantity(expr)
