from __future__ import annotations

import math
import re

import numpy as np
import pytest

import mcb.core as core
from mcb import (
    constrained_mcb_correct_useful_probability,
    unconstrained_mcb_edwards_hsu_correct_useful_probability,
)

# Hsu (1996) Appendix C.4 alloy-steel example.
BOOK = {"k": 8, "delta_star": 5000.0, "sigma": 3000.0, "alpha": 0.05}
PROB_FUNCS = [
    constrained_mcb_correct_useful_probability,
    unconstrained_mcb_edwards_hsu_correct_useful_probability,
]


# Reference values from scipy.integrate.quad, independent of the repo quadrature;
# they match Figure C.1 by eye.
@pytest.mark.parametrize(
    ("n", "expected"),
    [(2, 0.0391), (3, 0.1218), (4, 0.3640), (5, 0.7104), (6, 0.9088), (7, 0.9477)],
)
def test_constrained_probability_matches_independent_quadrature(n, expected):
    assert constrained_mcb_correct_useful_probability(n=n, **BOOK) == pytest.approx(expected, abs=5e-4)


@pytest.mark.parametrize(
    ("n", "expected"),
    [(4, 0.1218), (5, 0.3515), (6, 0.6728), (7, 0.8844), (8, 0.9434)],
)
def test_edwards_hsu_probability_matches_independent_quadrature(n, expected):
    assert unconstrained_mcb_edwards_hsu_correct_useful_probability(n=n, **BOOK) == pytest.approx(expected, abs=5e-4)


@pytest.mark.parametrize("func", PROB_FUNCS)
def test_probability_approaches_one_minus_alpha(func):
    assert func(n=10_000, **BOOK) == pytest.approx(0.95, abs=1e-4)


@pytest.mark.parametrize("func", PROB_FUNCS)
def test_probability_stays_at_ceiling_for_enormous_n(func):
    p = func(k=2, n=10**9, delta_star=5000.0, sigma=3000.0, alpha=0.05)
    assert p == pytest.approx(0.95, abs=1e-4)


def test_probability_uses_chi2_tail_beyond_asymptotic_switch():
    n = 600_000  # k = 2 gives nu = 1_199_998, above ASYMPTOTIC_NU_SWITCH
    d_inf = core.constrained_critical_value((1 / math.sqrt(2),), math.inf, 0.05)
    ratio = d_inf * math.sqrt(2 / n)  # makes u == 1
    p = constrained_mcb_correct_useful_probability(k=2, n=n, delta_star=ratio, sigma=1.0, alpha=0.05)
    assert 0.40 < p < 0.55  # a step function would give 0 or 0.95


@pytest.mark.parametrize("func", PROB_FUNCS)
def test_probability_is_zero_when_margin_is_tiny(func):
    assert func(k=8, n=2, delta_star=1e-6, sigma=1.0, alpha=0.05) == 0.0


def test_probability_accepts_numpy_integers():
    p = constrained_mcb_correct_useful_probability(
        k=np.int64(8), n=np.int64(6), delta_star=5000.0, sigma=3000.0, alpha=0.05
    )
    assert p == pytest.approx(0.9088, abs=5e-4)


@pytest.mark.parametrize("func", PROB_FUNCS)
@pytest.mark.parametrize(
    ("override", "message"),
    [
        ({"k": 1}, "k must be at least 2."),
        ({"n": 1}, "n must be at least 2."),
        ({"delta_star": 0.0}, "delta_star must be positive and finite."),
        ({"delta_star": math.inf}, "delta_star must be positive and finite."),
        ({"sigma": -1.0}, "sigma must be positive and finite."),
        ({"sigma": math.nan}, "sigma must be positive and finite."),
        ({"alpha": 0.0}, "alpha must be in (0, 1)."),
        ({"alpha": 1.0}, "alpha must be in (0, 1)."),
    ],
)
def test_probability_rejects_bad_inputs(func, override, message):
    with pytest.raises(ValueError, match=re.escape(message)):
        func(**({"n": 6, **BOOK} | override))
