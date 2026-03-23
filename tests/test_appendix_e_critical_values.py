from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from mcb import constrained_critical_value, unconstrained_edwards_hsu_critical_value

APPENDIX_E_PATH = Path(__file__).resolve().parent / "data" / "appendix_e_tables.json"
BALANCED_LAMBDA = 1.0 / math.sqrt(2.0)
TABLE_FUNCTIONS = {
    "E.2": constrained_critical_value,
    "E.3": unconstrained_edwards_hsu_critical_value,
    "E.5": constrained_critical_value,
    "E.6": unconstrained_edwards_hsu_critical_value,
}
TABLE_TOLERANCES = {
    "E.2": 0.0015,
    "E.3": 0.0015,
    "E.5": 0.005,
    "E.6": 0.005,
}


def _load_appendix_e_tables() -> dict[str, dict]:
    return json.loads(APPENDIX_E_PATH.read_text())


def _to_float_nu(value: int | str) -> float:
    return math.inf if value == "inf" else float(value)


def _build_cases() -> list[object]:
    cases = []
    for table_id in ("E.2", "E.3", "E.5", "E.6"):
        table = _load_appendix_e_tables()[table_id]
        alpha = float(table["alpha"])
        for nu, row in zip(table["nu_values"], table["matrix"]):
            nu_value = _to_float_nu(nu)
            nu_label = "inf" if nu == "inf" else str(int(nu))
            for k, expected in zip(table["k_values"], row):
                cases.append(
                    pytest.param(
                        table_id,
                        alpha,
                        nu_value,
                        int(k),
                        float(expected),
                        id=f"{table_id}-k{k}-nu{nu_label}",
                    )
                )
    return cases


@pytest.mark.parametrize(("table_id", "alpha", "nu", "k", "expected"), _build_cases())
def test_critical_values_match_appendix_e(
    table_id: str,
    alpha: float,
    nu: float,
    k: int,
    expected: float,
) -> None:
    func = TABLE_FUNCTIONS[table_id]
    tolerance = TABLE_TOLERANCES[table_id]
    lambdas = (BALANCED_LAMBDA,) * (k - 1)
    computed = float(func(lambdas, nu, alpha))
    rounded_match = round(computed, 3) == round(expected, 3)
    within_tolerance = abs(computed - expected) <= tolerance

    assert rounded_match or within_tolerance, (
        f"{table_id} k={k} nu={nu!r}: expected {expected:.3f}, got {computed:.6f} (tol={tolerance})"
    )
