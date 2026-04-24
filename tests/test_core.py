from __future__ import annotations

import math

import pytest
from scipy.stats import t

import mcb.core as core
from mcb import (
    MCBResult,
    constrained_critical_value,
    constrained_mcb,
    pooled_mse_from_group_stats,
    pooled_mse_from_samples,
    print_result,
    summarize_from_means_ns_and_mse,
    summarize_from_samples,
    unconstrained_edwards_hsu_critical_value,
    unconstrained_mcb_edwards_hsu,
)


def test_pooled_mse_from_samples_matches_formula() -> None:
    samples = {
        "A": [1.0, 3.0, 5.0],
        "B": [2.0, 4.0, 6.0],
    }

    pooled_mse = pooled_mse_from_samples(samples)

    assert pooled_mse == pytest.approx(4.0)


def test_pooled_mse_from_group_stats_matches_raw_samples() -> None:
    pooled_mse = pooled_mse_from_group_stats(ns=[3, 3], variances=[4.0, 4.0])

    assert pooled_mse == pytest.approx(4.0)


@pytest.mark.parametrize(
    ("func", "kwargs", "message"),
    [
        (
            pooled_mse_from_samples,
            {"samples": {"A": [1.0, 2.0]}},
            "Need at least two variants to compute pooled MSE",
        ),
        (
            pooled_mse_from_samples,
            {"samples": {"A": [1.0, 2.0], "B": [3.0]}},
            "Each variant must have at least two observations",
        ),
        (
            pooled_mse_from_group_stats,
            {"ns": [3], "variances": [1.0]},
            "Need at least two variants to compute pooled MSE",
        ),
        (
            pooled_mse_from_group_stats,
            {"ns": [3, 3], "variances": [1.0]},
            "ns and variances must have the same length",
        ),
        (
            pooled_mse_from_group_stats,
            {"ns": [1, 3], "variances": [1.0, 1.0]},
            "Each variant must have at least two observations",
        ),
        (
            pooled_mse_from_group_stats,
            {"ns": [3, 3], "variances": [1.0, -1.0]},
            "variances must be nonnegative",
        ),
    ],
)
def test_pooled_mse_helpers_validate_inputs(func, kwargs: dict, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        func(**kwargs)


def test_summarize_from_samples_computes_pooled_summary() -> None:
    summary = summarize_from_samples(
        {
            "A": [1.0, 3.0, 5.0],
            "B": [2.0, 4.0, 6.0],
        }
    )

    assert summary.labels == ("A", "B")
    assert summary.means == pytest.approx((3.0, 4.0))
    assert summary.ns == (3, 3)
    assert summary.df == 4
    assert summary.mse == pytest.approx(pooled_mse_from_samples({"A": [1.0, 3.0, 5.0], "B": [2.0, 4.0, 6.0]}))


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"means": [1.0], "ns": [3], "mse": 1.0}, "Need at least two variants"),
        ({"means": [1.0, 2.0], "ns": [3], "mse": 1.0}, "means and ns must have the same length"),
        ({"means": [1.0, 2.0], "ns": [1, 3], "mse": 1.0}, "Each variant must have at least two observations"),
        ({"means": [1.0, 2.0], "ns": [3, 3], "mse": 0.0}, "mse must be positive"),
        (
            {"means": [1.0, 2.0], "ns": [3, 3], "mse": 1.0, "labels": ["A"]},
            "labels must have the same length as means",
        ),
    ],
)
def test_summarize_from_means_ns_and_mse_validates_inputs(kwargs: dict, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        summarize_from_means_ns_and_mse(**kwargs)


def test_constrained_mcb_two_variant_case_matches_expected_bounds() -> None:
    summary = summarize_from_means_ns_and_mse(
        means=[1.0, 0.0],
        ns=[10, 10],
        mse=1.0,
        labels=["A", "B"],
    )

    result = constrained_mcb(summary, alpha=0.05)
    critical_value = constrained_critical_value((1.0 / math.sqrt(2.0),), float(summary.df), 0.05)
    standard_error = math.sqrt(1.0 / 10.0 + 1.0 / 10.0)

    assert critical_value == pytest.approx(t.ppf(0.95, df=18))
    assert result.estimates == pytest.approx((1.0, -1.0))
    assert result.lower == pytest.approx((0.0, -1.0 - critical_value * standard_error))
    assert result.upper == pytest.approx((1.0 + critical_value * standard_error, 0.0))
    assert result.candidate_best_set == ("A",)


def test_unconstrained_mcb_two_variant_case_matches_expected_bounds() -> None:
    summary = summarize_from_means_ns_and_mse(
        means=[1.0, 0.0],
        ns=[10, 10],
        mse=1.0,
        labels=["A", "B"],
    )

    result = unconstrained_mcb_edwards_hsu(summary, alpha=0.05)
    critical_value = unconstrained_edwards_hsu_critical_value((1.0 / math.sqrt(2.0),), float(summary.df), 0.05)
    standard_error = math.sqrt(1.0 / 10.0 + 1.0 / 10.0)

    assert critical_value == pytest.approx(t.ppf(0.975, df=18))
    assert result.estimates == pytest.approx((0.0, -1.0))
    assert result.lower == pytest.approx((0.0, -1.0 - critical_value * standard_error))
    assert result.upper == pytest.approx((0.0, -max(1.0 - critical_value * standard_error, 0.0)))
    assert result.candidate_best_set == ("A",)


def test_constrained_sat_example_matches_hsu_finite_df_values() -> None:
    summary = summarize_from_means_ns_and_mse(
        means=[619.0, 629.0, 575.0],
        ns=[103, 31, 122],
        mse=82.52**2,
        labels=["Computer science", "Engineering", "Other"],
    )

    result = constrained_mcb(summary, alpha=0.10)

    assert tuple(result.critical_values.values()) == pytest.approx(
        (1.605, 1.505, 1.612),
        abs=0.002,
    )
    assert result.lower == pytest.approx((-35.44, -17.14, -78.98), abs=0.03)
    assert result.upper == pytest.approx((17.14, 35.44, 0.0), abs=0.03)
    assert result.candidate_best_set == ("Computer science", "Engineering")


def test_effective_nu_uses_high_df_asymptotic_switch() -> None:
    just_below_switch = core.ASYMPTOTIC_NU_SWITCH - 1.0

    assert core._effective_nu(just_below_switch) == pytest.approx(just_below_switch)
    assert math.isinf(core._effective_nu(core.ASYMPTOTIC_NU_SWITCH))
    assert math.isinf(core._effective_nu(1.0e12))
    assert math.isinf(core._effective_nu(math.inf))


def test_critical_values_use_asymptotic_value_for_very_large_finite_df() -> None:
    lambdas = (1.0 / math.sqrt(2.0),) * 2

    for func in (constrained_critical_value, unconstrained_edwards_hsu_critical_value):
        finite_value = func(lambdas, 1.0e12, 0.05)
        asymptotic_value = func(lambdas, math.inf, 0.05)

        assert finite_value == pytest.approx(asymptotic_value)


def test_print_result_formats_output(capsys: pytest.CaptureFixture[str]) -> None:
    result = MCBResult(
        method="Example method",
        parameter="theta",
        alpha=0.05,
        labels=("A",),
        estimates=(0.1,),
        lower=(-0.2,),
        upper=(0.3,),
        critical_values={"A": 2.0},
        candidate_best_set=("A",),
    )

    print_result(result, digits=2)
    captured = capsys.readouterr()

    assert "Example method" in captured.out
    assert "Candidate-best set" in captured.out
    assert "A: 2.00" in captured.out
