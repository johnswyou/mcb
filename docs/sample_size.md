# Sample Size

The package can plan a balanced experiment before you run constrained MCB or Edwards-Hsu MCB. The method is Hsu's (Hsu 1988; Hsu 1996, Appendix C): choose the per-arm sample size $n$ so that, with a probability you pick, the MCB intervals are both **correct** (they cover the truth) and **useful** (they are narrow enough to conclude something).

```python
from mcb import (
    constrained_mcb_correct_useful_probability,
    constrained_mcb_sample_size,
    unconstrained_mcb_edwards_hsu_correct_useful_probability,
    unconstrained_mcb_edwards_hsu_sample_size,
)
```

## Why not F-test power

A common way to size an experiment is the power of the $F$-test. Hsu argues this is the wrong criterion for multiple comparisons. The power of the $F$-test includes outcomes where the confidence set misses the true parameter, and in those outcomes the method can make wrong assertions, such as $\mu_i > \mu_j$ when $\mu_i < \mu_j$. Power counts those as successes. Hsu's criterion counts an outcome as a success only if the intervals cover the truth *and* are narrow.

## What the guarantee says

You choose a margin $\delta^* > 0$ in the units of the response, and a planning value $\sigma$ for the common standard deviation. Let the interval half-width be $w = d\hat{\sigma}\sqrt{2/n}$ for constrained MCB, or $w = |d|\hat{\sigma}\sqrt{2/n}$ for Edwards-Hsu.

On the correct-and-useful event (the pivot event holds and $w < \delta^*$), both methods make two assertions, and both are correct:

1. Every arm more than $2\delta^*$ below the true best is ruled out. For constrained MCB its upper bound is 0, so it leaves the candidate-best set $G$. For Edwards-Hsu it leaves the screening set $S$.
2. The sample-best arm is within $\delta^*$ of the true best: its lower bound is above $-\delta^*$.

For constrained MCB these are Hsu's statements in Appendix C.2. For Edwards-Hsu, Hsu says only that the computation is "essentially the same" as for two-sided MCC (Appendix C.3); both statements follow from the pivot event $E$ in the proof of Theorem 4.2.3.

Edwards-Hsu needs more $n$ than constrained MCB because $|d| > d$. In return, it gives two-sided bounds on how far each arm is below the best.

## Formulas

With $k$ arms and $n$ per arm, let $\nu = k(n-1)$, let $\gamma_\nu$ be the density of $S = \hat{\sigma}/\sigma$, and let

$$
u = \sqrt{n/2}\,\frac{\delta^*/\sigma}{d}
$$

(with $|d|$ in place of $d$ for Edwards-Hsu). The critical values are the balanced ones, with every $\lambda = 1/\sqrt{2}$.

`constrained_mcb_correct_useful_probability` computes Hsu's Eq. (C.3):

$$
\int_0^u \int_{-\infty}^{\infty} \left[\Phi(z + \sqrt{2}\,d\,s)\right]^{k-1} d\Phi(z)\,\gamma_\nu(s)\,ds.
$$

`unconstrained_mcb_edwards_hsu_correct_useful_probability` computes Eq. (C.5):

$$
\int_0^u \int_{-\infty}^{\infty} \left[\Phi(z + \sqrt{2}\,|d|\,s) - \Phi(z - \sqrt{2}\,|d|\,s)\right]^{k-1} d\Phi(z)\,\gamma_\nu(s)\,ds.
$$

The upper limit $u$ is the narrowness condition $w < \delta^*$ rewritten as $S < u$. Both values are conservative: they are lower bounds on the probability of coverage and narrowness.

The `*_sample_size` functions return the smallest $n \ge 2$ at which the matching probability reaches `target`.

## Choosing $\delta^*$ and $\sigma$

- $\delta^*$ is the gap to the best that you are willing to treat as "as good as the best". It is half the gap $\delta = 2\delta^*$ at which an arm is guaranteed to be ruled out.
- $\sigma$ is a guess made before the data exist, for example from a pilot study or past experiments. Only the ratio $\delta^*/\sigma$ matters. If the real $\sigma$ turns out larger, the achieved probability is lower than planned.

## The $1 - \alpha$ ceiling

If the outer integral ran over $(0, \infty)$, it would equal $1 - \alpha$, because that is how the critical value is defined. The cut at $u$ only removes probability, so the probability of correct and useful inference is at most $1 - \alpha$ (up to the critical-value solver's tolerance of about $2 \times 10^{-5}$), and approaches it as $n$ grows. `target` must therefore be below $1 - \alpha$; the functions raise `ValueError` otherwise. Hsu's example uses $\alpha = 0.05$ with a target of $0.90$.

## Scope

The sample-size computation covers balanced designs only, as Appendix C does: one common $n$ for every arm. The analysis functions (`constrained_mcb`, `unconstrained_mcb_edwards_hsu`) still accept unequal sample sizes.

## Numerical note

The critical values switch to their infinite-degrees-of-freedom form at $\nu = 10^6$ (see {doc}`mcb_critical_values`). Above the switch, the narrowness condition uses the exact chi-square distribution of $S$ rather than treating $\hat{\sigma}$ as exactly $\sigma$; the coverage part is evaluated at $\hat{\sigma} = \sigma$, which changes the probability by about $10^{-4}$. The critical value just below the switch differs from the infinite-df one by slightly more than the true finite-df correction, within the solver's tolerance. As a result the probability jumps by about $10^{-3}$ across the switch for $k \ge 3$. The probability still increases with $n$, and a sample size near the switch is off by a few units out of $10^5$ to $10^6$.

## Worked example

Hsu (1996, Appendix C.4) takes Scheffé's design problem: compare the tensile strength of eight alloys ($k = 8$), with $\sigma$ expected to be about 3000 psi. Take $\delta^* = 5000$ psi, $\alpha = 0.05$ and a target of $0.90$:

```python
constrained_mcb_sample_size(k=8, delta_star=5000, sigma=3000, alpha=0.05, target=0.90)              # 6
unconstrained_mcb_edwards_hsu_sample_size(k=8, delta_star=5000, sigma=3000, alpha=0.05, target=0.90)  # 8

constrained_mcb_correct_useful_probability(k=8, n=6, delta_star=5000, sigma=3000, alpha=0.05)        # about 0.909
```

Six alloy specimens each are enough to select an alloy within 5000 psi of the true best with probability at least 0.90, matching Hsu's Figure C.1. Edwards-Hsu needs eight each.
