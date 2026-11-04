"""Hsu's Multiple Comparisons with the Best (MCB) for the unbalanced one-way model.

This module implements:

1. Constrained MCB for parameters
       theta_i = mu_i - max_{j != i} mu_j
   using Hsu's unbalanced-design construction.

2. Unconstrained MCB (Edwards-Hsu) for parameters
       eta_i = mu_i - max_j mu_j
   using Hsu's unbalanced-design extension.

3. Sample-size computation for balanced designs (Appendix C): the probability
   of correct and useful inference, and the smallest per-arm n reaching a
   target, for both methods above.

The critical-value computation follows the formulas and quadrature strategy
in Appendix D.2 of Hsu (1996):
- 48-point Gauss-Hermite quadrature for the inside integral.
- 48-point Gauss-Legendre quadrature for the outside integral.
- The outer integral on (0, infinity) is truncated to the central
  1 - 1e-6 probability interval of S = hat(sigma) / sigma.
- For nu >= 1e6, the code uses the infinite-df approximation because the
  finite-df correction is below the numerical solver tolerance and the finite
  quadrature path becomes unstable for extremely concentrated S.
- Root-finding uses the Illinois variant of regula falsi, matching the
  bracketing philosophy described by Hsu.

The code is written for revenue/session A/B/n settings, but the statistical
model is the standard homoscedastic one-way normal model. For zero-inflated
revenue/session data, exact finite-sample coverage is not guaranteed; the
method is most defensible when arm sizes are large enough that sample means
are approximately normal.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import inf, isinf, lgamma, log, sqrt
from typing import Callable, Mapping, Sequence

import math
import numpy as np
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss
from scipy.stats import chi2, norm, t

__all__ = [
    "MCBResult",
    "OneWaySummary",
    "constrained_critical_value",
    "constrained_mcb",
    "constrained_mcb_correct_useful_probability",
    "constrained_mcb_sample_size",
    "pooled_mse_from_group_stats",
    "pooled_mse_from_samples",
    "print_result",
    "summarize_from_means_ns_and_mse",
    "summarize_from_samples",
    "unconstrained_edwards_hsu_critical_value",
    "unconstrained_mcb_edwards_hsu",
    "unconstrained_mcb_edwards_hsu_correct_useful_probability",
    "unconstrained_mcb_edwards_hsu_sample_size",
]


# -----------------------------------------------------------------------------
# Numerical constants (chosen to mirror Appendix D.2)
# -----------------------------------------------------------------------------
GAUSS_HERMITE_POINTS = 48
GAUSS_LEGENDRE_POINTS = 48
OUTER_TRUNCATION_MASS = 1.0 - 1e-6
PROB_TOL = 2e-5
HALF_WIDTH_TOL = 5e-4
ASYMPTOTIC_NU_SWITCH = 1_000_000.0
MAX_ROOT_ITERS = 200

_GH_X, _GH_W = hermgauss(GAUSS_HERMITE_POINTS)
_GL_X, _GL_W = leggauss(GAUSS_LEGENDRE_POINTS)
_STD_NORMAL_MEASURE_NODES = np.sqrt(2.0) * _GH_X
_STD_NORMAL_MEASURE_WEIGHTS = _GH_W / np.sqrt(np.pi)


# -----------------------------------------------------------------------------
# Basic data containers
# -----------------------------------------------------------------------------
@dataclass(frozen=True)
class OneWaySummary:
    """Sufficient statistics for the pooled-variance one-way model.

    Parameters
    ----------
    labels:
        Variant names in display order.
    means:
        Sample mean revenue/session per variant.
    ns:
        Sample sizes per variant.
    mse:
        Pooled within-variant variance estimate.
    df:
        Error degrees of freedom, sum_i (n_i - 1).
    """

    labels: tuple[str, ...]
    means: tuple[float, ...]
    ns: tuple[int, ...]
    mse: float
    df: int

    @property
    def k(self) -> int:
        """Number of treatment arms in the summary."""
        return len(self.means)

    @property
    def sigma_hat(self) -> float:
        """Estimated common standard deviation, ``sqrt(mse)``."""
        return float(np.sqrt(self.mse))


@dataclass(frozen=True)
class MCBResult:
    """Output of a constrained or unconstrained MCB analysis.

    Attributes
    ----------
    method:
        Human-readable description of the procedure that produced the result.
    parameter:
        Target parameter family, for example ``mu_i - max_{j != i} mu_j``.
    alpha:
        Familywise error rate used to construct the simultaneous intervals.
    labels:
        Variant labels in display order.
    estimates:
        Point estimates for the target parameter, one per variant.
    lower:
        Simultaneous lower bounds for the target parameter.
    upper:
        Simultaneous upper bounds for the target parameter.
    critical_values:
        Mapping from label to the arm-specific critical value used internally.
    candidate_best_set:
        Labels that remain compatible with being best under the selected method.
    notes:
        Additional interpretation details about the result.
    """

    method: str
    parameter: str
    alpha: float
    labels: tuple[str, ...]
    estimates: tuple[float, ...]
    lower: tuple[float, ...]
    upper: tuple[float, ...]
    critical_values: Mapping[str, float]
    candidate_best_set: tuple[str, ...]
    notes: str = ""

    def as_dict(self) -> dict[str, dict[str, float | str]]:
        """Return the interval table keyed by variant label."""
        rows: dict[str, dict[str, float | str]] = {}
        for label, est, lo, hi in zip(self.labels, self.estimates, self.lower, self.upper):
            rows[label] = {
                "estimate": est,
                "lower": lo,
                "upper": hi,
            }
        return rows


# -----------------------------------------------------------------------------
# Model summaries
# -----------------------------------------------------------------------------
def pooled_mse_from_samples(samples: Mapping[str, Sequence[float]]) -> float:
    """Compute pooled within-group MSE from raw observations.

    Parameters
    ----------
    samples:
        Dict mapping variant label -> iterable of observations.

    Returns
    -------
    float
        Pooled within-group mean squared error.
    """
    if len(samples) < 2:
        raise ValueError("Need at least two variants to compute pooled MSE.")

    arrays = [np.asarray(group, dtype=float) for group in samples.values()]
    ns = tuple(int(arr.size) for arr in arrays)

    if any(n < 2 for n in ns):
        raise ValueError("Each variant must have at least two observations.")

    sse = sum(float(((arr - arr.mean()) ** 2).sum()) for arr in arrays)
    df = int(sum(n - 1 for n in ns))
    return float(sse / df)


def pooled_mse_from_group_stats(ns: Sequence[int], variances: Sequence[float]) -> float:
    """Compute pooled within-group MSE from per-group sample variances.

    Parameters
    ----------
    ns:
        Sample sizes per variant.
    variances:
        Per-variant sample variances computed with denominator ``n_i - 1``.

    Returns
    -------
    float
        Pooled within-group mean squared error.
    """
    ns = tuple(int(n) for n in ns)
    variances = tuple(float(v) for v in variances)

    if len(ns) != len(variances):
        raise ValueError("ns and variances must have the same length.")
    if len(ns) < 2:
        raise ValueError("Need at least two variants to compute pooled MSE.")
    if any(n < 2 for n in ns):
        raise ValueError("Each variant must have at least two observations.")
    if any(v < 0 for v in variances):
        raise ValueError("variances must be nonnegative.")

    numerator = sum((n - 1) * variance for n, variance in zip(ns, variances))
    denominator = sum(n - 1 for n in ns)
    return float(numerator / denominator)


def summarize_from_samples(samples: Mapping[str, Sequence[float]]) -> OneWaySummary:
    """Create a OneWaySummary from raw per-variant observations.

    Parameters
    ----------
    samples:
        Dict mapping variant label -> iterable of revenue/session observations.

    Returns
    -------
    OneWaySummary
        Sample means, group sizes, pooled MSE, and pooled error degrees of freedom.
    """
    if len(samples) < 2:
        raise ValueError("Need at least two variants for MCB.")

    labels = tuple(samples.keys())
    arrays = [np.asarray(samples[label], dtype=float) for label in labels]
    ns = tuple(int(arr.size) for arr in arrays)

    if any(n < 2 for n in ns):
        raise ValueError("Each variant must have at least two observations.")

    means = tuple(float(arr.mean()) for arr in arrays)
    df = int(sum(n - 1 for n in ns))
    mse = pooled_mse_from_samples(samples)
    return OneWaySummary(labels=labels, means=means, ns=ns, mse=mse, df=df)


def summarize_from_means_ns_and_mse(
    means: Sequence[float],
    ns: Sequence[int],
    mse: float,
    labels: Sequence[str] | None = None,
) -> OneWaySummary:
    """Create a ``OneWaySummary`` from means, sample sizes, and pooled MSE.

    Parameters
    ----------
    means:
        Per-arm sample means.
    ns:
        Per-arm sample sizes.
    mse:
        Pooled within-arm mean squared error.
    labels:
        Optional display labels. Default labels are generated when omitted.

    Returns
    -------
    OneWaySummary
        Sample means, group sizes, pooled MSE, and pooled error degrees of freedom.
    """
    means = tuple(float(x) for x in means)
    ns = tuple(int(n) for n in ns)
    if len(means) != len(ns):
        raise ValueError("means and ns must have the same length.")
    if len(means) < 2:
        raise ValueError("Need at least two variants for MCB.")
    if any(n < 2 for n in ns):
        raise ValueError("Each variant must have at least two observations.")
    if mse <= 0:
        raise ValueError("mse must be positive.")

    if labels is None:
        labels = tuple(f"Variant {i+1}" for i in range(len(means)))
    labels = tuple(labels)
    if len(labels) != len(means):
        raise ValueError("labels must have the same length as means.")

    df = int(sum(n - 1 for n in ns))
    return OneWaySummary(labels=labels, means=means, ns=ns, mse=float(mse), df=df)


# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------
def _positive_part(x: float | np.ndarray) -> float | np.ndarray:
    return np.maximum(x, 0.0)


def _negative_part_magnitude(x: float | np.ndarray) -> float | np.ndarray:
    # Hsu's x^- denotes the magnitude of the negative part, max(-x, 0).
    return np.maximum(-np.asarray(x), 0.0)


def _effective_nu(nu: float) -> float:
    return inf if isinf(nu) or nu >= ASYMPTOTIC_NU_SWITCH else nu


def _t_upper_quantile(beta: float, nu: float) -> float:
    if isinf(nu):
        return float(norm.ppf(1.0 - beta))
    return float(t.ppf(1.0 - beta, df=nu))


def _gamma_sigmahat_over_sigma_density(s: np.ndarray, nu: float) -> np.ndarray:
    # Density of S = sqrt(ChiSquare_nu / nu), s > 0.
    log_coef = log(2.0) + (nu / 2.0) * log(nu / 2.0) - lgamma(nu / 2.0)
    return np.exp(log_coef + (nu - 1.0) * np.log(s) - nu * np.square(s) / 2.0)


def _central_s_interval(nu: float) -> tuple[float, float]:
    tail = (1.0 - OUTER_TRUNCATION_MASS) / 2.0
    lo = sqrt(chi2.ppf(tail, df=nu) / nu)
    hi = sqrt(chi2.ppf(1.0 - tail, df=nu) / nu)
    return lo, hi


def _integrate_wrt_standard_normal_measure(values_at_nodes: np.ndarray) -> float:
    return float(np.dot(_STD_NORMAL_MEASURE_WEIGHTS, values_at_nodes))


# -----------------------------------------------------------------------------
# Appendix D.2 quadrature-based probability evaluators
# -----------------------------------------------------------------------------
def _outer_integral(inside_for_s: Callable[[float], float], nu: float, s_upper: float = inf) -> float:
    """Integrate ``inside_for_s`` against the density of S = sigma_hat / sigma on (0, s_upper)."""
    nu_eff = _effective_nu(nu)
    if isinf(nu_eff):
        h1 = inside_for_s(1.0)
        if isinf(s_upper):
            return h1
        if isinf(nu):  # sigma known: S == 1
            return h1 if s_upper >= 1.0 else 0.0
        # Not a step function: at nu ~ 1e6, S still has sd ~ 7e-4.
        return h1 * float(chi2.cdf(nu * s_upper * s_upper, df=nu))

    a, b = _central_s_interval(nu_eff)
    b = min(b, s_upper)
    if b <= a:
        return 0.0
    s_nodes = 0.5 * (b - a) * _GL_X + 0.5 * (b + a)
    s_weights = 0.5 * (b - a) * _GL_W
    density = _gamma_sigmahat_over_sigma_density(s_nodes, nu_eff)
    integrand = np.array([inside_for_s(s) for s in s_nodes]) * density
    return float(np.dot(s_weights, integrand))


def _probability_one_sided_dunnett(
    lambdas: np.ndarray, critical_value: float, nu: float, s_upper: float = inf
) -> float:
    """Probability in Hsu's Eq. (4.11), with the outer integral cut at ``s_upper``."""
    lambdas = np.asarray(lambdas, dtype=float)
    denom = np.sqrt(1.0 - np.square(lambdas))

    def inside_for_s(s: float) -> float:
        z = _STD_NORMAL_MEASURE_NODES
        args = (lambdas[:, None] * z[None, :] + critical_value * s) / denom[:, None]
        prod = norm.cdf(args).prod(axis=0)
        return _integrate_wrt_standard_normal_measure(prod)

    return _outer_integral(inside_for_s, nu, s_upper)


def _probability_two_sided_dunnett(
    lambdas: np.ndarray, critical_value: float, nu: float, s_upper: float = inf
) -> float:
    """Edwards-Hsu / two-sided MCC probability, with the outer integral cut at ``s_upper``."""
    lambdas = np.asarray(lambdas, dtype=float)
    denom = np.sqrt(1.0 - np.square(lambdas))

    def inside_for_s(s: float) -> float:
        z = _STD_NORMAL_MEASURE_NODES
        upper = (lambdas[:, None] * z[None, :] + critical_value * s) / denom[:, None]
        lower = (lambdas[:, None] * z[None, :] - critical_value * s) / denom[:, None]
        prod = (norm.cdf(upper) - norm.cdf(lower)).prod(axis=0)
        return _integrate_wrt_standard_normal_measure(prod)

    return _outer_integral(inside_for_s, nu, s_upper)


# -----------------------------------------------------------------------------
# Illinois / regula falsi root finder (mirrors Appendix D.2 bracketing logic)
# -----------------------------------------------------------------------------
def _illinois_root(
    f,
    a: float,
    b: float,
    *,
    prob_tol: float = PROB_TOL,
    half_width_tol: float = HALF_WIDTH_TOL,
    max_iter: int = MAX_ROOT_ITERS,
) -> float:
    fa = float(f(a))
    fb = float(f(b))
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    if fa * fb > 0.0:
        raise RuntimeError("Root is not bracketed.")

    last_replaced = 0  # -1 if left replaced last, +1 if right replaced last
    c = a

    for _ in range(max_iter):
        c = (a * fb - b * fa) / (fb - fa)
        fc = float(f(c))

        if abs(fc) <= prob_tol or 0.5 * abs(b - a) <= half_width_tol:
            return c

        if fa * fc < 0.0:
            # Root is in [a, c]; replace b.
            b, fb = c, fc
            if last_replaced == 1:
                fa *= 0.5
            last_replaced = 1
        elif fb * fc < 0.0:
            # Root is in [c, b]; replace a.
            a, fa = c, fc
            if last_replaced == -1:
                fb *= 0.5
            last_replaced = -1
        else:
            return c

    raise RuntimeError("Illinois method did not converge within max_iter.")


# -----------------------------------------------------------------------------
# Critical values
# -----------------------------------------------------------------------------
@lru_cache(maxsize=512)
def constrained_critical_value(lambdas: tuple[float, ...], nu: float, alpha: float) -> float:
    """Compute the constrained MCB critical value ``d^i``.

    Parameters
    ----------
    lambdas:
        Tuple of lambda_ij = (1 + n_i / n_j)^(-1/2) over j != i.
    nu:
        Error degrees of freedom.
    alpha:
        Familywise error rate.

    Returns
    -------
    float
        One-sided Dunnett-with-control critical value for the focal arm.

    Notes
    -----
    When there are only two arms, the result reduces to the upper ``alpha``
    quantile of the corresponding ``t`` distribution.
    """
    lambdas_arr = np.asarray(lambdas, dtype=float)
    k = lambdas_arr.size + 1
    nu_eff = _effective_nu(float(nu))

    lower = _t_upper_quantile(alpha, nu_eff)
    upper = _t_upper_quantile(alpha / (k - 1), nu_eff)

    if k == 2:
        return lower

    target = 1.0 - alpha
    f = lambda d: _probability_one_sided_dunnett(lambdas_arr, d, nu_eff) - target
    return _illinois_root(f, lower, upper)


@lru_cache(maxsize=512)
def unconstrained_edwards_hsu_critical_value(lambdas: tuple[float, ...], nu: float, alpha: float) -> float:
    """Compute the Edwards-Hsu critical value ``|d|^i``.

    Parameters
    ----------
    lambdas:
        Tuple of lambda_ij = (1 + n_i / n_j)^(-1/2) over j != i.
    nu:
        Error degrees of freedom.
    alpha:
        Familywise error rate.

    Returns
    -------
    float
        Two-sided Dunnett-with-control critical value for the focal arm.

    Notes
    -----
    When there are only two arms, the result reduces to the upper ``alpha / 2``
    quantile of the corresponding ``t`` distribution.
    """
    lambdas_arr = np.asarray(lambdas, dtype=float)
    k = lambdas_arr.size + 1
    nu_eff = _effective_nu(float(nu))

    lower = _t_upper_quantile(alpha / 2.0, nu_eff)
    upper = _t_upper_quantile(alpha / (2.0 * (k - 1)), nu_eff)

    if k == 2:
        return lower

    target = 1.0 - alpha
    f = lambda d: _probability_two_sided_dunnett(lambdas_arr, d, nu_eff) - target
    return _illinois_root(f, lower, upper)


# -----------------------------------------------------------------------------
# Main MCB procedures
# -----------------------------------------------------------------------------
def _pairwise_se(summary: OneWaySummary, i: int, j: int) -> float:
    return summary.sigma_hat * math.sqrt(1.0 / summary.ns[i] + 1.0 / summary.ns[j])


def _lambdas_for_control(summary: OneWaySummary, control_index: int) -> tuple[float, ...]:
    ni = summary.ns[control_index]
    return tuple((1.0 + ni / summary.ns[j]) ** (-0.5) for j in range(summary.k) if j != control_index)


def constrained_mcb(summary: OneWaySummary, alpha: float = 0.05) -> MCBResult:
    """Hsu's constrained MCB for the unbalanced one-way model.

    The target parameter family is ``theta_i = mu_i - max_{j != i} mu_j``.

    A positive lower bound is impossible; each interval contains 0 when the
    corresponding arm might still be best.

    Parameters
    ----------
    summary:
        Pooled one-way summary statistics.
    alpha:
        Familywise error rate for the simultaneous intervals.

    Returns
    -------
    MCBResult
        Constrained MCB estimates, simultaneous bounds, critical values, and
        the candidate-best set.
    """
    k = summary.k
    means = np.asarray(summary.means, dtype=float)

    d_values = [
        constrained_critical_value(_lambdas_for_control(summary, i), float(summary.df), alpha)
        for i in range(k)
    ]

    # Upper bounds D_i^+
    d_plus = []
    for i in range(k):
        bounds = [
            means[i] - means[j] + d_values[i] * _pairwise_se(summary, i, j)
            for j in range(k)
            if j != i
        ]
        d_plus.append(float(_positive_part(min(bounds))))

    G = [i for i, up in enumerate(d_plus) if up > 0.0]

    # Lower bounds D_i^-
    d_minus: list[float] = []
    for i in range(k):
        if G == [i]:
            d_minus.append(0.0)
            continue

        candidate_js = [j for j in G if j != i]
        if not candidate_js:
            # This only happens when G == [i], already handled above.
            d_minus.append(0.0)
            continue

        bounds = [
            means[i] - means[j] - d_values[j] * _pairwise_se(summary, i, j)
            for j in candidate_js
        ]
        d_minus.append(float(min(bounds)))

    estimates = [
        means[i] - max(means[j] for j in range(k) if j != i)
        for i in range(k)
    ]

    return MCBResult(
        method="Hsu constrained MCB (unbalanced one-way)",
        parameter="mu_i - max_{j != i} mu_j",
        alpha=alpha,
        labels=summary.labels,
        estimates=tuple(float(x) for x in estimates),
        lower=tuple(d_minus),
        upper=tuple(d_plus),
        critical_values={summary.labels[i]: float(d_values[i]) for i in range(k)},
        candidate_best_set=tuple(summary.labels[i] for i in G),
        notes=(
            "Candidate-best set G consists of variants with strictly positive "
            "upper constrained MCB bounds."
        ),
    )


def unconstrained_mcb_edwards_hsu(summary: OneWaySummary, alpha: float = 0.05) -> MCBResult:
    """Unbalanced Edwards-Hsu unconstrained MCB.

    The direct target parameter family is ``eta_i = mu_i - max_j mu_j``, which
    is always ``<= 0`` and measures how far variant ``i`` is below the global
    best.

    This is the cleanest direct unconstrained MCB construction in Hsu's chapter.

    Parameters
    ----------
    summary:
        Pooled one-way summary statistics.
    alpha:
        Familywise error rate for the simultaneous intervals.

    Returns
    -------
    MCBResult
        Edwards-Hsu estimates, simultaneous bounds, critical values, and the
        candidate-best screening set.
    """
    k = summary.k
    means = np.asarray(summary.means, dtype=float)

    dabs_values = [
        unconstrained_edwards_hsu_critical_value(_lambdas_for_control(summary, i), float(summary.df), alpha)
        for i in range(k)
    ]

    # Candidate-best set S
    S = []
    for i in range(k):
        vals = [
            means[i] - means[j] + dabs_values[i] * _pairwise_se(summary, i, j)
            for j in range(k)
            if j != i
        ]
        if min(vals) > 0.0:
            S.append(i)

    lower = []
    upper = []
    for i in range(k):
        L_ij = []
        U_ij = []
        for j in S:
            if i == j:
                L_ij.append(0.0)
                U_ij.append(0.0)
            else:
                x = means[i] - means[j]
                se = _pairwise_se(summary, i, j)
                L_ij.append(float(x - dabs_values[j] * se))
                U_ij.append(float(-_negative_part_magnitude(x + dabs_values[j] * se)))

        lower.append(float(min(L_ij)))
        upper.append(float(max(U_ij)))

    estimates = [float(means[i] - np.max(means)) for i in range(k)]

    return MCBResult(
        method="Edwards-Hsu unconstrained MCB (unbalanced one-way)",
        parameter="mu_i - max_j mu_j",
        alpha=alpha,
        labels=summary.labels,
        estimates=tuple(estimates),
        lower=tuple(lower),
        upper=tuple(upper),
        critical_values={summary.labels[i]: float(dabs_values[i]) for i in range(k)},
        candidate_best_set=tuple(summary.labels[i] for i in S),
        notes=(
            "Candidate-best set S consists of variants whose one-variant-as-control "
            "two-sided MCC upper screens remain positive against every rival."
        ),
    )


# -----------------------------------------------------------------------------
# Sample size (Hsu 1996, Appendix C)
# -----------------------------------------------------------------------------
def _validate_design_inputs(k: int, delta_star: float, sigma: float, alpha: float) -> None:
    if k < 2:
        raise ValueError("k must be at least 2.")
    if not (math.isfinite(delta_star) and delta_star > 0):
        raise ValueError("delta_star must be positive and finite.")
    if not (math.isfinite(sigma) and sigma > 0):
        raise ValueError("sigma must be positive and finite.")
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must be in (0, 1).")


def _correct_useful_probability(
    k: int, n: int, delta_star: float, sigma: float, alpha: float, *, two_sided: bool
) -> float:
    k, n = int(k), int(n)
    _validate_design_inputs(k, delta_star, sigma, alpha)
    if n < 2:
        raise ValueError("n must be at least 2.")

    nu = float(k * (n - 1))
    lambdas = (1.0 / sqrt(2.0),) * (k - 1)
    if two_sided:
        crit = unconstrained_edwards_hsu_critical_value(lambdas, nu, alpha)
        prob = _probability_two_sided_dunnett
    else:
        crit = constrained_critical_value(lambdas, nu, alpha)
        prob = _probability_one_sided_dunnett
    if crit <= 0.0:
        raise ValueError("alpha is too large: the critical value is not positive.")
    # Narrowness crit * sigma_hat * sqrt(2/n) < delta_star  <=>  S < u.
    u = sqrt(n / 2.0) * (delta_star / sigma) / crit
    return prob(np.asarray(lambdas), crit, nu, s_upper=u)


def constrained_mcb_correct_useful_probability(
    k: int, n: int, delta_star: float, sigma: float, alpha: float = 0.05
) -> float:
    """Hsu's probability of correct and useful inference for balanced constrained MCB.

    Computes Hsu (1996), Theorem C.2.1, Eq. (C.3): the probability that, in a
    balanced experiment with ``n`` observations per arm, the pivot event from
    the proof of Theorem 4.1.1 holds and the half-width
    ``d * sigma_hat * sqrt(2 / n)`` is below ``delta_star``. That event implies
    the constrained MCB intervals cover their targets and are that narrow, so
    the value is a conservative lower bound on the probability of correct and
    useful inference. When the event happens:

    - every arm more than ``2 * delta_star`` below the true best gets upper
      bound 0, so it leaves the candidate-best set;
    - the sample-best arm is correctly asserted to be within ``delta_star``
      of the true best.

    The value is at most ``1 - alpha``, up to the critical-value solver
    tolerance, and approaches it as ``n`` grows.

    Parameters
    ----------
    k:
        Number of arms.
    n:
        Common per-arm sample size (balanced design).
    delta_star:
        Indifference margin, in the units of the response.
    sigma:
        Planning value (a guess) for the common standard deviation.
    alpha:
        Familywise error rate of the intervals.

    Returns
    -------
    float
        Probability of correct and useful inference.
    """
    return _correct_useful_probability(k, n, delta_star, sigma, alpha, two_sided=False)


def unconstrained_mcb_edwards_hsu_correct_useful_probability(
    k: int, n: int, delta_star: float, sigma: float, alpha: float = 0.05
) -> float:
    """Hsu's probability of correct and useful inference for balanced Edwards-Hsu MCB.

    Hsu (1996), Appendix C, sizes Edwards-Hsu MCB with the two-sided MCC
    computation, Eq. (C.5): the probability that the pivot event of Theorem
    4.2.3 holds and the half-width ``|d| * sigma_hat * sqrt(2 / n)`` is below
    ``delta_star``. That event implies the Edwards-Hsu intervals cover their
    targets and are that narrow, so the value is a conservative lower bound on
    the probability of correct and useful inference. When the event happens:

    - every arm more than ``2 * delta_star`` below the true best leaves the
      screening set S;
    - the sample-best arm is correctly asserted to be within ``delta_star``
      of the true best.

    The value is at most ``1 - alpha``, up to the critical-value solver
    tolerance, and approaches it as ``n`` grows. It needs more ``n`` than constrained
    MCB because ``|d| > d``, in exchange for two-sided bounds on each arm's
    gap to the best.

    Parameters
    ----------
    k:
        Number of arms.
    n:
        Common per-arm sample size (balanced design).
    delta_star:
        Indifference margin, in the units of the response.
    sigma:
        Planning value (a guess) for the common standard deviation.
    alpha:
        Familywise error rate of the intervals.

    Returns
    -------
    float
        Probability of correct and useful inference.
    """
    return _correct_useful_probability(k, n, delta_star, sigma, alpha, two_sided=True)


_MAX_DESIGN_N = 2**31


def _smallest_n(prob: Callable[[int], float], target: float, alpha: float) -> int:
    # Any dip in P(n) starts at n = 2: P may fall at first, then only rises. So once
    # P(2) < target, the set {n : P(n) >= target} is upward-closed and bisection finds its minimum.
    if prob(2) >= target:
        return 2
    lo, hi = 2, 4
    while prob(hi) < target:
        if hi >= _MAX_DESIGN_N:
            raise RuntimeError(
                f"target={target} not reached by n={hi} at alpha={alpha}: either "
                "delta_star is too small relative to sigma, or target is within the "
                "critical-value solver tolerance of 1 - alpha."
            )
        lo, hi = hi, 2 * hi
    while hi - lo > 1:  # invariant: prob(lo) < target <= prob(hi)
        mid = (lo + hi) // 2
        if prob(mid) >= target:
            hi = mid
        else:
            lo = mid
    return hi


def _sample_size(
    prob_func: Callable[..., float], k: int, delta_star: float, sigma: float, alpha: float, target: float
) -> int:
    _validate_design_inputs(int(k), delta_star, sigma, alpha)
    if not 0.0 < target < 1.0 - alpha:
        raise ValueError(
            "target must be in (0, 1 - alpha); the probability of correct and useful "
            "inference never reaches 1 - alpha."
        )
    return _smallest_n(lambda n: prob_func(k, n, delta_star, sigma, alpha), target, alpha)


def constrained_mcb_sample_size(
    k: int, delta_star: float, sigma: float, alpha: float = 0.05, target: float = 0.90
) -> int:
    """Smallest balanced per-arm sample size for constrained MCB (Hsu 1996, Appendix C).

    Returns the smallest ``n >= 2`` per arm for which
    :func:`constrained_mcb_correct_useful_probability` is at least ``target``.
    With at least that probability, every arm more than ``2 * delta_star``
    below the true best is ruled out, and the sample-best arm is correctly
    asserted to be within ``delta_star`` of the true best.

    Parameters
    ----------
    k:
        Number of arms.
    delta_star:
        Indifference margin, in the units of the response.
    sigma:
        Planning value (a guess) for the common standard deviation.
    alpha:
        Familywise error rate of the intervals.
    target:
        Required probability of correct and useful inference; must be below
        ``1 - alpha``.

    Returns
    -------
    int
        Common per-arm sample size. Call the probability function at this
        ``n`` to see the probability actually achieved.
    """
    return _sample_size(constrained_mcb_correct_useful_probability, k, delta_star, sigma, alpha, target)


def unconstrained_mcb_edwards_hsu_sample_size(
    k: int, delta_star: float, sigma: float, alpha: float = 0.05, target: float = 0.90
) -> int:
    """Smallest balanced per-arm sample size for Edwards-Hsu MCB (Hsu 1996, Appendix C).

    Returns the smallest ``n >= 2`` per arm for which
    :func:`unconstrained_mcb_edwards_hsu_correct_useful_probability` is at
    least ``target``. With at least that probability, every arm more than
    ``2 * delta_star`` below the true best leaves the screening set, and the
    sample-best arm is correctly asserted to be within ``delta_star`` of the
    true best.

    Parameters
    ----------
    k:
        Number of arms.
    delta_star:
        Indifference margin, in the units of the response.
    sigma:
        Planning value (a guess) for the common standard deviation.
    alpha:
        Familywise error rate of the intervals.
    target:
        Required probability of correct and useful inference; must be below
        ``1 - alpha``.

    Returns
    -------
    int
        Common per-arm sample size. Call the probability function at this
        ``n`` to see the probability actually achieved.
    """
    return _sample_size(
        unconstrained_mcb_edwards_hsu_correct_useful_probability, k, delta_star, sigma, alpha, target
    )


# -----------------------------------------------------------------------------
# Pretty-print helper for demos
# -----------------------------------------------------------------------------
def print_result(result: MCBResult, digits: int = 4) -> None:
    """Pretty-print an ``MCBResult`` for demos or exploratory work.

    Parameters
    ----------
    result:
        The result object to display.
    digits:
        Number of decimal places used for interval endpoints and critical values.
    """
    print(result.method)
    print(f"Parameter: {result.parameter}")
    print(f"Familywise alpha: {result.alpha}")
    print(f"Candidate-best set: {result.candidate_best_set}")
    print("Critical values:")
    for label, d in result.critical_values.items():
        print(f"  {label}: {d:.{digits}f}")
    print("Intervals:")
    for label, est, lo, hi in zip(result.labels, result.estimates, result.lower, result.upper):
        print(f"  {label:<18s} estimate={est:.{digits}f}  CI=[{lo:.{digits}f}, {hi:.{digits}f}]")
