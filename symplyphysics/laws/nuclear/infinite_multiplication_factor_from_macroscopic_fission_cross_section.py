from sympy import (Eq, solve, S)
from symplyphysics import (
    units,
    expr_to_quantity,
    Quantity,
    Symbol,
    print_expression,
    Dimensionless,
    convert_to,
    validate_input_symbols,
)

# Description
## Infinite multiplication factor: k_infinite = v * Σf / Σa
## Where:
## v - average number of neutrons produced per fission.
## Σf - overall macroscopic fission cross-section.
##   See [macroscopic cross-section](./macroscopic_cross_section_from_free_mean_path.py) implementation.
## Σa - overall macroscopic absorption cross-section.
## k_infinite - infinite multiplication factor.
##   See [infinite multiplication factor](./infinite_multiplication_factor.py)

neutrons_per_fission = Symbol("neutrons_per_fission", Dimensionless)
macroscopic_fission_cross_section = Symbol("macroscopic_fission_cross_section", 1 / units.length)
macroscopic_absorption_cross_section = Symbol("macroscopic_absorption_cross_section",
    1 / units.length)
infinite_multiplication_factor = Symbol("infinite_multiplication_factor", Dimensionless)

law = Eq(
    infinite_multiplication_factor,
    neutrons_per_fission * macroscopic_fission_cross_section / macroscopic_absorption_cross_section)


def print() -> str:
    return print_expression(law)


@validate_input_symbols(neutrons_per_fission_=neutrons_per_fission,
    macroscopic_fission_cross_section_=macroscopic_fission_cross_section,
    macroscopic_absorption_cross_section_=macroscopic_absorption_cross_section)
def calculate_multiplication_factor(neutrons_per_fission_: float,
    macroscopic_fission_cross_section_: Quantity,
    macroscopic_absorption_cross_section_: Quantity) -> float:

    result_factor_expr = solve(law, infinite_multiplication_factor,
        dict=True)[0][infinite_multiplication_factor]
    result_expr = result_factor_expr.subs({
        neutrons_per_fission: neutrons_per_fission_,
        macroscopic_fission_cross_section: macroscopic_fission_cross_section_,
        macroscopic_absorption_cross_section: macroscopic_absorption_cross_section_
    })
    result_factor = expr_to_quantity(result_expr)
    return convert_to(result_factor, S.One).evalf()
