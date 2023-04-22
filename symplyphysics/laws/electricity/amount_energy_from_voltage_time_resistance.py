from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input_symbols, validate_output_symbol)

# Description
# The amount of energy released by a conductor with a current is directly proportional
# to the square of the applied voltage, the time of the current and inversely proportional
# to the resistance of the conductor. This is the differential form of the law

# Amount of energy Q = U**2 * t / R
# where:
# U - voltage to conductor
# t - time of current action
# R - resistance of conductor
# The resultant energy is the energy generated by the passage of current through the conductor.

amount_energy = Symbol("amount_energy", units.energy)
voltage = Symbol("voltage", units.voltage)
time = Symbol("time", units.time)
resistance = Symbol("resistance", units.impedance)

law = Eq(amount_energy, (voltage**2 * time) / resistance)


def print() -> str:
    return print_expression(law)


@validate_input_symbols(voltage_=voltage, time_=time, resistance_=resistance)
@validate_output_symbol(amount_energy)
def calculate_amount_energy(voltage_: Quantity, time_: Quantity, resistance_: Quantity) -> Quantity:
    result_energy_expr = solve(law, amount_energy, dict=True)[0][amount_energy]
    result_expr = result_energy_expr.subs({voltage: voltage_, time: time_, resistance: resistance_})
    return expr_to_quantity(result_expr)
