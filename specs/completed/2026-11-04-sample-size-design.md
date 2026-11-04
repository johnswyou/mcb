# Sample-size determination for MCB: design

Status: draft for review, 2026-11-04.
Source: Hsu (1996), Appendix C (C.2, C.3, C.4) and Appendix D.3; Hsu (1988).

## Intent

The package can analyze an experiment, but it can't plan one. Users who want to know "how many units per arm?" must guess or use F-test power. Hsu argues F-test power is the wrong criterion: it counts outcomes where the intervals miss the truth, including outcomes with wrong directional assertions, as successes.

This feature adds Hsu's criterion. Pick the per-arm sample size `n` so that, with probability at least `target` (Hsu's 1 − β), the MCB intervals both **cover the truth** and are **narrow enough to be useful**.

## What the user is buying

The user supplies a margin δ* > 0 (Hsu's indifference margin) and a planning value σ for the common standard deviation. Let `w` be the interval half-width: `d·σ̂·√(2/n)` for constrained MCB, `|d|·σ̂·√(2/n)` for Edwards–Hsu.

On the "correct and useful" event (simultaneous coverage and `w < δ*`), both methods deliver the same two assertions:

1. Every arm more than **2δ*** below the true best is ruled out. For constrained MCB its upper bound is 0, so it leaves `G`. For Edwards–Hsu it leaves `S`.
2. The sample-best arm is correctly asserted to be within **δ*** of the true best (its lower bound exceeds −δ*).

For constrained MCB these two statements come from the book (C.2). For Edwards–Hsu the book only says the computation is "essentially the same as" two-sided MCC (C.3). I derived both statements from the pivot event `E` in the proof of Theorem 4.2.3; the builder should re-check them before they go into the docs.

- (1): under `E`, μ̂_i − μ̂_(k) + w < η_i + 2w, which is < 0 when η_i < −2δ*, so i ∉ S.
- (2): for every j ∈ S, L_[k]j = μ̂_[k] − μ̂_j − w ≥ −w, and −w > −δ*.

Edwards–Hsu costs more `n` than constrained MCB because |d| > d. In return it gives two-sided bounds on how far each arm is below the best.

## Goals

- For constrained MCB and Edwards–Hsu, provide:
  - (a) the probability of correct and useful inference at a given `n`, the equivalent of Hsu's `pmcb` program;
  - (b) the smallest `n` that reaches `target`.
- Reproduce the book's example (C.4): k = 8, σ = 3000, δ* = 5000, α = 0.05, target = 0.90 gives **n = 6** for constrained MCB and **n = 8** for Edwards–Hsu. The Edwards–Hsu value is the book's P{C} curve; the integrand is identical, see Math.
- Handle A/B-scale problems in a few seconds. Small δ*/σ means `n` in the 10⁴–10⁸ range per arm. The module docstring targets revenue/session A/B/n tests, so this matters (my inference, not stated by you).

## Non-goals

- **Unequal allocations.** Appendix C covers only the balanced model (C.1), and the book gives no theory for unequal `n`. The API takes a common `n`.
- Tukey MCA (C.1) and Dunnett MCC (C.3) sample sizes. The package is MCB-only. Two-sided MCC uses the same integral as Edwards–Hsu, so it's a one-line addition later if wanted.
- A separate known-σ (indifference-zone, n ≥ 2(dσ/δ*)²) API. Large `n` already reaches it through the ν → ∞ path.
- Two-stage designs (Matejcik and Nelson 1995), plotting the probability curve, a CLI, and a new example script.

## Anti-goals

- Calling this "power", or offering F-test-power sizing.
- Replacing the exact integral with Tukey's Bonferroni shortcut (narrowness at 1 − (β − α)) or any other approximation that silently changes `n`.
- Changing any existing critical value. The Appendix E regression tests must pass untouched.

## Constraints

- Dependencies stay numpy and scipy only.
- Reuse the existing Appendix D.2 quadrature in `src/mcb/core.py`; don't add a second integrator.
- The docs build on Read the Docs with `fail_on_warning: true`. Every new `.md` under `docs/` must be in a toctree. This design file therefore lives in `specs/`, outside the Sphinx tree.
- Locally, Sphinx may also pick up the gitignored `docs/references/` copy of the book, because `docs/conf.py` doesn't exclude it. The builder should check whether `sphinx -W` fails locally for that reason before blaming the new pages. Adding `"references"` to `exclude_patterns` is acceptable if needed.

## Math (what gets computed)

Balanced design with k arms and n per arm:

```
ν  = k(n − 1)
d  = constrained_critical_value((1/√2,)*(k−1), ν, α)              # constrained
|d|= unconstrained_edwards_hsu_critical_value((1/√2,)*(k−1), ν, α) # Edwards–Hsu
u  = √(n/2) · (δ*/σ) / d          (use |d| for Edwards–Hsu)

P_constrained(n) = ∫₀ᵘ ∫ Φ(z + √2·d·s)^(k−1)                      dΦ(z) γ_ν(s) ds    (C.3)
P_edwards_hsu(n) = ∫₀ᵘ ∫ [Φ(z + √2|d|s) − Φ(z − √2|d|s)]^(k−1)    dΦ(z) γ_ν(s) ds    (C.5)
```

γ_ν is the density of S = σ̂/σ = √(χ²_ν/ν).

With λ = 1/√2, the existing evaluators' integrand Φ((λz + ds)/√(1 − λ²)) equals Φ(z + √2·d·s) exactly. So C.3 is `_probability_one_sided_dunnett` and C.5 is `_probability_two_sided_dunnett`, each with the outer integral cut at `u` instead of ∞.

Over (0, ∞) each integral equals 1 − α by the definition of d. So P(n) < 1 − α for every n, and P(n) → 1 − α as n → ∞. `target` must be < 1 − α.

These probabilities are lower bounds on P(coverage and narrowness), since that event contains the pivot event. Sizing on them is conservative, as in the book.

## Public API

Four functions in `src/mcb/core.py`, added to `__all__` there and in `src/mcb/__init__.py`:

```python
def constrained_mcb_correct_useful_probability(
    k: int, n: int, delta_star: float, sigma: float, alpha: float = 0.05
) -> float: ...

def unconstrained_mcb_edwards_hsu_correct_useful_probability(
    k: int, n: int, delta_star: float, sigma: float, alpha: float = 0.05
) -> float: ...

def constrained_mcb_sample_size(
    k: int, delta_star: float, sigma: float, alpha: float = 0.05, target: float = 0.90
) -> int: ...

def unconstrained_mcb_edwards_hsu_sample_size(
    k: int, delta_star: float, sigma: float, alpha: float = 0.05, target: float = 0.90
) -> int: ...
```

- Names follow the existing `constrained_*` / `unconstrained_*_edwards_hsu*` pattern.
- `*_sample_size` returns the smallest per-arm `n ≥ 2` with probability ≥ `target`. To see the achieved probability, call the matching `*_probability` function at that `n`.
- Docstrings use the numpy style the module already uses. Each states the two assertions from "What the user is buying", including the δ* versus 2δ* distinction, and notes that σ is a planning guess.

## Numerics

**Outer-limit cut.** Add a keyword `s_upper: float = inf` to `_probability_one_sided_dunnett` and `_probability_two_sided_dunnett`. With the default, behavior is byte-for-byte unchanged.

- Finite ν path: Gauss–Legendre over `[a, min(b, s_upper)]`, where `[a, b]` is `_central_s_interval`. Return 0 if `s_upper ≤ a`.
- Asymptotic path (`_effective_nu` returns ∞): return `inside_for_s(1.0) * P(S < s_upper)`.
  - For finite raw ν, P(S < u) = `chi2.cdf(ν·u², ν)`.
  - For ν = ∞, it is 1 if u ≥ 1, else 0.

Don't use the plain step function 1[u ≥ 1] when ν is large but finite. At ν ≈ 10⁶, S still has sd ≈ 7·10⁻⁴, which is wider than the change in `u` from one step in `n`. A step function would bias the answer.

**Search for n.**

1. If P(2) ≥ target, return 2.
2. Otherwise double `n` (2, 4, 8, …) until P(hi) ≥ target.
3. Bisect on integers between `lo` (P < target) and `hi` (P ≥ target).

Correctness argument: P(n) isn't monotone. My spike scanned k ∈ {2, 3, 5, 8}, δ*/σ ∈ {0.25, 0.5, 1, 5/3, 3} and n ≤ 60, both methods. P only decreases at the very start (n from 2 to about 6), after which it increases. Every dip starts at n = 2; the final review measured dips from P(2) up to 0.048 (k = 2), so the dips are small but not below 0.02 as first thought. This matches theory: P(S < u) falls with ν while u < ½. So when P(2) < target, the set {n : P(n) ≥ target} is upward-closed, and bisection returns the smallest n.

**Cap.** If doubling passes an internal cap (2³¹ is fine), raise `RuntimeError` naming the target and α. This can happen when `target` is within the critical-value solver tolerance (`PROB_TOL = 2e-5`) of 1 − α. The cap is not a parameter.

**Known small jump at the ν switch.** The reviewer measured this. Just below `ASYMPTOTIC_NU_SWITCH`, the existing solver's `d` sits about 1e-5 above d(∞) at α = 0.05, and about 4e-5 above at α = 0.01. That is within the solver's tolerance, but the cut probability magnifies it by about √(2ν) ≈ 1400.

- For k ≥ 3, P jumps by about 1e-3 across the switch. Ordinary steps between consecutive n are about 2e-4 to 1.5e-3.
- P stays monotone and bisection stays correct. Near the switch, the returned `n` is a few units off out of 10⁵–10⁶.
- Accept this and mention it in the docs. Don't tighten the critical-value solver or give sizing its own `d`, since that would break the anti-goal or make sizing and analysis use different `d`.

**Cost.** The spike measured about 40 ms per evaluation, including the critical value. An A/B-scale search needs about 40–60 evaluations, so about 2–3 s. The `lru_cache` on critical values (512 entries) is enough.

## Validation

Raise `ValueError`, in the message style of `summarize_from_means_ns_and_mse`, for:

- `k < 2`
- `n < 2`
- `delta_star` or `sigma` not finite and positive
- `alpha` not in (0, 1)
- `target` not in (0, 1 − alpha). The message says the probability can never reach 1 − α, so a larger target is unreachable.

## Docs

- New page `docs/sample_size.md`. It covers: the criterion and why not F-test power; the two assertions and δ* vs 2δ*; formulas C.3 and C.5; choosing δ* and a planning σ; the 1 − α ceiling; the balanced-only scope; and the alloy-steel example (C.4) with code reproducing 6 and 8.
- Add the page to the toctree in `docs/index.md`. Add the four functions to the autosummary in `docs/api.md`. Add a "planning" line to the index's task list and highlights.
- `README.md`: a short "Planning sample size" snippet.
- `tests/README.md`: list the new test file.

## Tests: `tests/test_sample_size.py`

1. **Book example:** constrained gives 6 and Edwards–Hsu gives 8 (k = 8, δ* = 5000, σ = 3000, α = 0.05, target = 0.90).
2. **Probability curve:** check against reference values from an independent integrator (`scipy.integrate.quad`), never from the code under test. Constrained, k = 8, δ*/σ = 5/3, α = 0.05, n = 2…7, already computed with `quad`: 0.039, 0.122, 0.364, 0.710, 0.909, 0.948, with tolerance about 1e-3. For Edwards–Hsu, compute the equivalent with `quad`. The spike, using the repo quadrature, gave n = 4…8: 0.122, 0.352, 0.673, 0.884, 0.943, which matches Figure C.1 by eye.
3. **Ceiling:** P at very large n is within 1e-4 of 1 − α, and `target ≥ 1 − alpha` raises `ValueError`.
4. **Search equals linear scan** on a small grid that covers k = 2 and the non-monotone start. Keep the test under about 10 s.
5. **A/B scale:** pick inputs where the answer lands near `ASYMPTOTIC_NU_SWITCH`, for k = 2 and at least one k ≥ 3. Check:
   - P(n) ≥ target > P(n − 1);
   - P is non-decreasing across the switch;
   - the step across the switch is no more than about 5× the neighbouring steps (see "Known small jump" in Numerics). An absolute limit won't work, because ordinary steps grow with k.
6. **Validation:** parametrized `ValueError` cases, following the pattern in `test_core.py`.
7. **Existing tests** pass unchanged. That proves the default `s_upper` path is untouched.

## Calls I made (please check)

1. Balanced designs only.
2. MCB only. No Tukey MCA or Dunnett MCC.
3. Inputs are `delta_star` and `sigma` in raw units, rather than a single ratio δ*/σ.
4. The margin is Hsu's δ* (the within-δ* selection margin), not δ = 2δ*.
5. `*_sample_size` returns a plain `int`. Probabilities come from separate public functions, not a result dataclass.
6. Code goes in `core.py`, next to the quadrature it reuses, rather than in a new `sample_size.py`.
7. The function names, especially `*_correct_useful_probability`.
8. Defaults are α = 0.05 and target = 0.90, the book's example.
9. Large-ν handling uses `h(1)·chi2.cdf`, not a step function.
10. This file lives in `specs/` so it stays out of the Read the Docs build.
11. Accept the small probability jump at ν = 10⁶ (a few units of `n` out of 10⁵–10⁶) instead of changing the critical-value solver.

## Left to the builder

- Exact test grids and tolerances within the limits above.
- Where in the toctree the new page goes, and the page's wording.
- Whether `k` and `n` are coerced with `int()` (as `summarize_from_means_ns_and_mse` does) or rejected when not integers.
- The exact doubling cap constant.
