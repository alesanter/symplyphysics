from functools import reduce
import operator
from typing import Any, Self
from sympy import Expr, Basic, flatten


class InvSumArray(Expr):
    """
    Represents unevaluated Sum of inverted quantities over array.

    """

    def __new__(cls, *args: Any) -> Self:
        # We remove all levels of nested tuples, because we want to sum everything we got.
        flat_args = tuple(flatten(args))
        obj = Expr.__new__(cls, *flat_args)
        return obj

    def _eval_nseries(self, x: Any, n: Any, logx: Any, cdir: Any) -> Any:
        pass

    def doit(self, **_hints: Any) -> Basic:
        return reduce(operator.add, map(lambda x: 1 / x, self.args))
