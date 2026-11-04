from __future__ import annotations

import math
import re

import numpy as np
import pytest

import mcb.core as core
from mcb import (
    constrained_mcb_correct_useful_probability,
    constrained_mcb_sample_size,
    unconstrained_mcb_edwards_hsu_correct_useful_probability,
    unconstrained_mcb_edwards_hsu_sample_size,
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


@pytest.mark.parametrize("func", PROB_FUNCS)
def test_probability_handles_absurdly_large_margin(func):
    # nu = 2e7 takes the chi-square path; u ~ 1e163 would overflow u**2 as a Python float.
    assert func(k=2, n=10**7, delta_star=1e160, sigma=1.0, alpha=0.05) == pytest.approx(0.95, abs=1e-4)


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


PAIRS = [
    (constrained_mcb_sample_size, constrained_mcb_correct_useful_probability),
    (unconstrained_mcb_edwards_hsu_sample_size, unconstrained_mcb_edwards_hsu_correct_useful_probability),
]


@pytest.mark.parametrize(("size_func", "expected"), [(PAIRS[0][0], 6), (PAIRS[1][0], 8)])
def test_sample_size_reproduces_hsu_appendix_c4(size_func, expected):
    assert size_func(target=0.90, **BOOK) == expected


def _linear_scan(prob, target):
    n = 2
    while prob(n) < target:
        n += 1
    return n


# (2, 0.5, 0.01): constrained P(2) = 0.018 already meets target (returns 2);
# Edwards-Hsu P(2) = 0.0035 dips then rises, exercising the non-monotone start.
@pytest.mark.parametrize(("size_func", "prob_func"), PAIRS)
@pytest.mark.parametrize(
    ("k", "ratio", "target"),
    [(2, 0.5, 0.01), (2, 0.5, 0.9), (3, 1.0, 0.8), (8, 5 / 3, 0.5), (8, 5 / 3, 0.9)],
)
def test_search_matches_linear_scan(size_func, prob_func, k, ratio, target):
    expected = _linear_scan(lambda n: prob_func(k=k, n=n, delta_star=ratio, sigma=1.0), target)
    assert size_func(k=k, delta_star=ratio, sigma=1.0, target=target) == expected


@pytest.mark.parametrize(
    ("size_func", "prob_func", "crit_func"),
    [
        (*PAIRS[0], core.constrained_critical_value),
        (*PAIRS[1], core.unconstrained_edwards_hsu_critical_value),
    ],
)
@pytest.mark.parametrize("k", [2, 3])
def test_sample_size_near_asymptotic_switch(size_func, prob_func, crit_func, k):
    n_sw = math.ceil(core.ASYMPTOTIC_NU_SWITCH / k) + 1  # first n with k(n-1) >= switch
    d_inf = crit_func((1 / math.sqrt(2),) * (k - 1), math.inf, 0.05)
    ratio = d_inf * math.sqrt(2 / n_sw)  # u ~ 1 at the switch, so P ~ 0.475

    def p(n):
        return prob_func(k=k, n=n, delta_star=ratio, sigma=1.0)

    p0, p1, p2, p3 = (p(n_sw + i) for i in (-2, -1, 0, 1))
    across = p2 - p1
    assert across >= 0.0
    # Spec says "about 5x"; measured ratio is about 4-5 for k = 3, so 10x leaves headroom.
    assert across <= 10 * max(p1 - p0, p3 - p2)
    n = size_func(k=k, delta_star=ratio, sigma=1.0, target=0.45)
    assert p(n) >= 0.45 > p(n - 1)


def test_sample_size_returns_two_when_margin_is_large():
    assert constrained_mcb_sample_size(k=3, delta_star=10.0, sigma=1.0) == 2


def test_target_ceiling_follows_alpha():
    assert constrained_mcb_sample_size(k=8, delta_star=5000.0, sigma=3000.0, alpha=0.05, target=0.92) == 7
    with pytest.raises(ValueError, match=re.escape("target must be in (0, 1 - alpha)")):
        constrained_mcb_sample_size(k=8, delta_star=5000.0, sigma=3000.0, alpha=0.10, target=0.92)


@pytest.mark.parametrize("size_func", [pair[0] for pair in PAIRS])
@pytest.mark.parametrize("target", [0.0, 0.95, 0.99])
def test_sample_size_rejects_unreachable_target(size_func, target):
    with pytest.raises(ValueError, match=re.escape("target must be in (0, 1 - alpha)")):
        size_func(target=target, **BOOK)


def test_sample_size_validates_inputs_before_searching():
    with pytest.raises(ValueError, match=re.escape("k must be at least 2.")):
        constrained_mcb_sample_size(k=1, delta_star=1.0, sigma=1.0)


def test_smallest_n_raises_at_cap():
    # A tiny margin hits the cap too, so the message must not blame only the target.
    with pytest.raises(RuntimeError, match=re.escape("delta_star is too small relative to sigma")):
        core._smallest_n(lambda n: 0.0, target=0.5, alpha=0.05)


# For k = 2, constrained d = t quantile at alpha, which is 0 at alpha = 0.5 and negative above.
@pytest.mark.parametrize("alpha", [0.5, 0.6])
def test_probability_rejects_alpha_with_nonpositive_critical_value(alpha):
    with pytest.raises(ValueError, match=re.escape("alpha is too large")):
        constrained_mcb_correct_useful_probability(k=2, n=600_000, delta_star=1.0, sigma=1.0, alpha=alpha)


def test_sample_size_rejects_alpha_with_nonpositive_critical_value():
    with pytest.raises(ValueError, match=re.escape("alpha is too large")):
        constrained_mcb_sample_size(k=2, delta_star=1.0, sigma=1.0, alpha=0.6, target=0.3)
