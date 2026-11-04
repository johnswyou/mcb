# MCB Sample-Size Determination Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add Hsu's (1996, Appendix C) "probability of correct and useful inference" and smallest-`n` search for constrained MCB and Edwards–Hsu MCB, for balanced designs.

**Architecture:** Everything lives in `src/mcb/core.py`, next to the Appendix D.2 quadrature it reuses. The two existing probability evaluators gain an optional upper cut `s_upper` on the outer integral. The duplicated outer-integration code moves into one shared helper. Two thin public probability functions sit on top. Two sample-size functions run a doubling-plus-bisection search over `n`.

**Tech Stack:** Python ≥ 3.10, numpy, scipy, pytest, Sphinx with myst-parser.

**Spec:** `specs/completed/2026-11-04-sample-size-design.md`. Read it first; it holds the math, the guarantees and the reasoning behind each decision below.

## Global Constraints

- No new dependencies; numpy and scipy only.
- Existing critical values must not change. The current suite (2,329 tests, about 2.5 min) must pass untouched.
- Balanced designs only. The API takes one common `n` per arm.
- Never name anything "power".
- `ValueError` messages are full sentences ending in a period, matching `summarize_from_means_ns_and_mse`.
- Read the Docs uses `fail_on_warning: true`. Every new `.md` under `docs/` must be in a toctree.
- Environment: system `python3` has no numpy. Set up once with `uv venv` then `uv pip install -e '.[docs,plot,test]'`, and run everything with `.venv/bin/python`. Do not create or commit `uv.lock`.
- Work on a feature branch, not `main`. End every commit message with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Review Focus

1. A margin tiny relative to σ at small `n` (outer cut `u` below the S interval) returns `0.0`, never NaN or a negative value. Pinned in Task 1.
2. numpy integers for `k` and `n` (as from arrays or dataframes) are accepted and give the same result. Pinned in Task 1.
3. Enormous `n` (ν around 2·10⁹) stays at the 1 − α ceiling, with no overflow or NaN from `chi2.cdf`. Pinned in Task 1.
4. A margin large relative to σ returns `n = 2` without searching. Pinned in Task 2.
5. α other than 0.05 moves the target ceiling: target 0.92 is rejected at α = 0.10 and gives 7 at α = 0.05 for the book example. Pinned in Task 2.

---

### Task 1: Cut outer integral and the two probability functions

**Files:**
- Modify: `src/mcb/core.py`
  - line 36: add `Callable` to the `typing` import
  - lines 343–388: the two evaluators
  - `__all__`
  - new section "Sample size (Hsu 1996, Appendix C)" just before `# Pretty-print helper for demos` (line ~690)
- Modify: `src/mcb/__init__.py` (imports and `__all__`)
- Create: `tests/test_sample_size.py`

**Interfaces:**
- Produces, in `mcb.core` (private):
  - `_outer_integral(inside_for_s: Callable[[float], float], nu: float, s_upper: float = inf) -> float`
  - `_probability_one_sided_dunnett(lambdas, critical_value, nu, s_upper: float = inf) -> float`
  - `_probability_two_sided_dunnett(lambdas, critical_value, nu, s_upper: float = inf) -> float`
  - `_validate_design_inputs(k: int, delta_star: float, sigma: float, alpha: float) -> None`
- Produces, public and exported from `mcb`:
  - `constrained_mcb_correct_useful_probability(k: int, n: int, delta_star: float, sigma: float, alpha: float = 0.05) -> float`
  - `unconstrained_mcb_edwards_hsu_correct_useful_probability(k: int, n: int, delta_star: float, sigma: float, alpha: float = 0.05) -> float`

- [ ] **Step 1: Write the failing tests** in `tests/test_sample_size.py`

```python
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
```

- [ ] **Step 2: Run them and confirm they fail**

Run: `.venv/bin/python -m pytest tests/test_sample_size.py -q`
Expected: collection error, `ImportError: cannot import name 'constrained_mcb_correct_useful_probability'`.

- [ ] **Step 3: Add `_outer_integral` and route both evaluators through it**

Both evaluators keep their own `inside_for_s` and take the new keyword `s_upper: float = inf`. Replace their identical outer-integration tails with `return _outer_integral(inside_for_s, nu, s_upper)`. With the default `s_upper`, the result must be bit-identical to today's code, because `min(b, inf) == b`.

```python
def _outer_integral(inside_for_s, nu, s_upper=inf):
    nu_eff = _effective_nu(nu)
    if isinf(nu_eff):
        h1 = inside_for_s(1.0)
        if isinf(s_upper):
            return h1
        if isinf(nu):                      # sigma known: S == 1
            return h1 if s_upper >= 1.0 else 0.0
        return h1 * float(chi2.cdf(nu * s_upper**2, df=nu))   # not a step: S still has spread
    a, b = _central_s_interval(nu_eff)
    b = min(b, s_upper)
    if b <= a:
        return 0.0
    # ... existing Gauss–Legendre nodes, weights and density over [a, b], unchanged
```

- [ ] **Step 4: Add `_validate_design_inputs` and the two public probability functions** in the new "Sample size (Hsu 1996, Appendix C)" section

- `_validate_design_inputs` raises, in this order:
  - `"k must be at least 2."`
  - `"delta_star must be positive and finite."`
  - `"sigma must be positive and finite."`
  - `"alpha must be in (0, 1)."`
  - Use `math.isfinite(x) and x > 0` so NaN fails.
- Each public function:
  1. Coerces `k, n = int(k), int(n)`, as `summarize_from_means_ns_and_mse` does.
  2. Validates, then raises `"n must be at least 2."` if needed.
  3. Computes `nu = float(k * (n - 1))`, `lambdas = (1 / sqrt(2),) * (k - 1)`, the critical value (`constrained_critical_value` or `unconstrained_edwards_hsu_critical_value`), and `u = sqrt(n / 2) * (delta_star / sigma) / crit`.
  4. Returns the matching evaluator called with `s_upper=u`.
- Share the body through one private helper with a `two_sided: bool` flag.
- Numpy-style docstrings must state:
  - the two assertions from the spec ("What the user is buying"), including the δ* versus 2δ* distinction;
  - that σ is a planning guess;
  - that the value is a conservative lower bound and is always below 1 − α.
- Add both names to `core.__all__`, and to the imports and `__all__` in `src/mcb/__init__.py`.

- [ ] **Step 5: Run the new tests**

Run: `.venv/bin/python -m pytest tests/test_sample_size.py -q`
Expected: all pass.

- [ ] **Step 6: Run the full suite to prove critical values are untouched**

Run: `.venv/bin/python -m pytest -q`
Expected: 2,329 existing tests plus the new ones pass.

- [ ] **Step 7: Commit**

```bash
git add src/mcb/core.py src/mcb/__init__.py tests/test_sample_size.py
git commit -m "Add MCB probability of correct and useful inference (Hsu App. C)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Smallest-n search

**Files:**
- Modify: `src/mcb/core.py`
  - Appendix C section from Task 1
  - `__all__`
  - module docstring: add "3. Sample-size computation for balanced designs (Appendix C)" to the numbered list
- Modify: `src/mcb/__init__.py`
- Test: `tests/test_sample_size.py`

**Interfaces:**
- Consumes: `_validate_design_inputs`, plus both `*_correct_useful_probability` functions from Task 1.
- Produces:
  - `_MAX_DESIGN_N = 2**31` (private constant)
  - `_smallest_n(prob: Callable[[int], float], target: float, alpha: float) -> int` (private)
  - `constrained_mcb_sample_size(k: int, delta_star: float, sigma: float, alpha: float = 0.05, target: float = 0.90) -> int` (public)
  - `unconstrained_mcb_edwards_hsu_sample_size(k: int, delta_star: float, sigma: float, alpha: float = 0.05, target: float = 0.90) -> int` (public)

- [ ] **Step 1: Write the failing tests**, appended to `tests/test_sample_size.py`, with the two new names added to the `from mcb import (...)` block

```python
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


@pytest.mark.parametrize("k", [2, 3])
def test_sample_size_near_asymptotic_switch(k):
    n_sw = math.ceil(core.ASYMPTOTIC_NU_SWITCH / k) + 1  # first n with k(n-1) >= switch
    d_inf = core.constrained_critical_value((1 / math.sqrt(2),) * (k - 1), math.inf, 0.05)
    ratio = d_inf * math.sqrt(2 / n_sw)  # u ~ 1 at the switch, so P ~ 0.475

    def p(n):
        return constrained_mcb_correct_useful_probability(k=k, n=n, delta_star=ratio, sigma=1.0)

    p0, p1, p2, p3 = (p(n_sw + i) for i in (-2, -1, 0, 1))
    across = p2 - p1
    assert across >= 0.0
    # Spec says "about 5x"; measured ratio is about 4-5 for k = 3, so 10x leaves headroom.
    assert across <= 10 * max(p1 - p0, p3 - p2)
    n = constrained_mcb_sample_size(k=k, delta_star=ratio, sigma=1.0, target=0.45)
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
    with pytest.raises(RuntimeError, match="not reached"):
        core._smallest_n(lambda n: 0.0, target=0.5, alpha=0.05)
```

- [ ] **Step 2: Run them and confirm they fail**

Run: `.venv/bin/python -m pytest tests/test_sample_size.py -q`
Expected: collection error, `ImportError: cannot import name 'constrained_mcb_sample_size'`.

- [ ] **Step 3: Implement `_MAX_DESIGN_N` and `_smallest_n`**

Correctness relies on P(n) falling only at the very start and only while P < 0.02 (spec, "Numerics: Search for n"). So when P(2) < target, the set of passing `n` is upward-closed.

```python
def _smallest_n(prob, target, alpha):
    if prob(2) >= target:
        return 2
    lo, hi = 2, 4
    while prob(hi) < target:
        if hi >= _MAX_DESIGN_N:
            raise RuntimeError(
                f"target={target} not reached by n={hi} at alpha={alpha}; "
                "target is too close to 1 - alpha for the critical-value tolerance."
            )
        lo, hi = hi, 2 * hi
    while hi - lo > 1:          # invariant: prob(lo) < target <= prob(hi)
        mid = (lo + hi) // 2
        if prob(mid) >= target:
            hi = mid
        else:
            lo = mid
    return hi
```

- [ ] **Step 4: Implement the two public sample-size functions**

Each function:
1. Calls `_validate_design_inputs(int(k), delta_star, sigma, alpha)`.
2. Raises `ValueError("target must be in (0, 1 - alpha); the probability of correct and useful inference never reaches 1 - alpha.")` unless `0 < target < 1 - alpha`.
3. Returns `_smallest_n(lambda n: <matching probability function>(k, n, delta_star, sigma, alpha), target, alpha)`.

Docstrings:
- point to the matching probability function for the achieved value;
- repeat the δ* meaning;
- say the design is balanced (per-arm `n`).

Add both names to `core.__all__` and to `src/mcb/__init__.py`, and update the module docstring.

- [ ] **Step 5: Run the new tests**

Run: `.venv/bin/python -m pytest tests/test_sample_size.py -q --durations=5`
Expected: all pass, and no single test over about 10 s.

- [ ] **Step 6: Run the full suite**

Run: `.venv/bin/python -m pytest -q`
Expected: all pass.

- [ ] **Step 7: Commit**

```bash
git add src/mcb/core.py src/mcb/__init__.py tests/test_sample_size.py
git commit -m "Add balanced MCB sample-size search (Hsu App. C)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Documentation

**Files:**
- Create: `docs/sample_size.md`
- Modify:
  - `docs/index.md`: toctree, task list, highlights
  - `docs/api.md`: autosummary
  - `README.md`: new section after the plotting paragraph in "Usage"
  - `tests/README.md`: contents list
- Maybe modify: `docs/conf.py` (`exclude_patterns`, only if Step 4 requires it)

**Interfaces:**
- Consumes: the four public functions from Tasks 1–2, with their exact names and signatures.

- [ ] **Step 1: Write `docs/sample_size.md`**

Use these sections, with math in the `$…$` / `$$…$$` style the other pages use:

1. **Why not F-test power.** Hsu's argument: F-test power counts wrong assertions as successes.
2. **What the guarantee says.** On the correct-and-useful event:
   - every arm more than 2δ* below the best is ruled out (leaves `G` or `S`);
   - the sample-best arm is within δ* of the best.

   The Edwards–Hsu version follows from the pivot event in Theorem 4.2.3; say so.
3. **Formulas.** C.3 and C.5, with ν = k(n − 1) and u = √(n/2)(δ*/σ)/d.
4. **Choosing δ* and a planning σ.**
5. **The 1 − α ceiling,** and why `target` must be below it.
6. **Scope.** Balanced designs only, and the analysis functions still accept unequal `n`.
7. **Numerical note.** The small jump at ν = 10⁶ (a few units of `n` out of 10⁵–10⁶).
8. **Worked example.** The alloy-steel case from Appendix C.4, using the README snippet below. Show outputs 6 and 8, plus `constrained_mcb_correct_useful_probability(...)` at n = 6 giving about 0.909.

- [ ] **Step 2: Wire up the index, API page, README and tests README**

- `docs/index.md`:
  - add `sample_size` to the toctree after `mcb_critical_values`;
  - add a bullet for planning to the "organized around" list;
  - add one highlight: "Plan balanced experiments with Hsu's sample-size computation."
- `docs/api.md`: add the four functions after `unconstrained_mcb_edwards_hsu` in the autosummary.
- `tests/README.md`: add `` `test_sample_size.py` — Appendix C probabilities, book example, search and switch checks``.
- `README.md`, new section "Planning sample size":

```python
from mcb import constrained_mcb_sample_size, unconstrained_mcb_edwards_hsu_sample_size

# Hsu (1996) Appendix C.4: 8 alloys, planning sigma 3000 psi, margin delta* = 5000 psi.
constrained_mcb_sample_size(k=8, delta_star=5000, sigma=3000, alpha=0.05, target=0.90)              # 6
unconstrained_mcb_edwards_hsu_sample_size(k=8, delta_star=5000, sigma=3000, alpha=0.05, target=0.90)  # 8
```

- [ ] **Step 3: Check that the snippet runs**

Run: `.venv/bin/python -c "from mcb import constrained_mcb_sample_size as c, unconstrained_mcb_edwards_hsu_sample_size as u; print(c(k=8, delta_star=5000, sigma=3000), u(k=8, delta_star=5000, sigma=3000))"`
Expected: `6 8`

- [ ] **Step 4: Build the docs with warnings as errors**

Run: `.venv/bin/python -m sphinx -W -b html docs docs/_build/html`
Expected: `build succeeded`.

If the build fails because Sphinx picked up the gitignored `docs/references/` copy of the book, add `"references"` to `exclude_patterns` in `docs/conf.py`, then rerun. Fix any warning coming from the new page itself.

- [ ] **Step 5: Commit**

```bash
git add docs/sample_size.md docs/index.md docs/api.md README.md tests/README.md docs/conf.py
git commit -m "Document MCB sample-size computation

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

`git add` of an unchanged `docs/conf.py` is a harmless no-op.
