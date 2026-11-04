# Sample Size

Before you launch an A/B/n test, you have to decide how many users each arm gets. With too few, the test ends with intervals so wide that you can't pick a winner or drop anything. With too many, you spend traffic and calendar time you didn't need. This page builds up, one step at a time, how `mcb` makes that choice for a test you will analyze with MCB. The method is Hsu's (Hsu 1988; Hsu 1996, Appendix C).

If you only want the recipe, skip to [Plan the test](#plan-the-test).

## A running example

A streaming service wants to test four new recommendation algorithms against the current one, so there are $k = 5$ arms. Users are split evenly and at random across the arms. The metric is minutes watched per user in the first week, and more is better. Last month's data put the standard deviation of that metric at about 40 minutes.

After the test, the team will run `constrained_mcb` on the results. For each arm $i$, it returns a confidence interval for

$$
\theta_i = \mu_i - \max_{j \neq i} \mu_j,
$$

how far arm $i$'s true mean is above or below the best of the other arms. All $k$ intervals are correct at the same time with probability at least $1 - \alpha$. An arm whose upper bound is 0 is ruled out as the best, and the arms left form the candidate-best set $G$. [Constrained MCB](constrained_mcb.md) has the details.

The question for this page is how many users each arm needs.

## What a successful test looks like

When the test ends, the team wants to be able to say two things:

1. Ship this one: the observed winner is at most a little worse than the true best arm.
2. Drop these: the arms that are clearly worse than the best are ruled out.

"A little worse" needs a number. You can't ask for the exact best arm. If two algorithms differ by 0.01 minutes, no test of practical size can tell them apart, and nobody would care which one shipped. So you pick a margin $\delta^* > 0$, in the units of the metric: the shortfall from the best that you are willing to treat as "as good as the best". Here, say $\delta^* = 2$ minutes.

The margin splits the arms into three zones, by how far each arm's true mean is below the best one:

```text
       clearly worse          gray zone       good enough
  <--------------------|-----------------|---------------|
                  best - 2δ*         best - δ*         best
```

A well-planned test promises, with high probability, that

- the observed winner is in the good-enough zone, and the intervals say so;
- every arm in the clearly-worse zone is ruled out.

It makes no promise about arms in the gray zone: they may or may not be ruled out. The factor of 2 is explained below.

## Two ways a test can let you down

The intervals can fail you in two different ways.

**They can be wrong.** By bad luck, the data mislead and some interval misses its true value. MCB caps the probability of this at $\alpha$, whatever the sample size. More data does not lower it. It is the price of a $1 - \alpha$ confidence level.

**They can be too wide to be useful.** The intervals are right, but so wide that nothing is ruled out and the winner's lower bound sits far below 0. More data fixes this.

How wide are they? In a balanced test with $n$ users per arm, every constrained MCB bound is an estimated difference between two arm means, plus or minus the same amount

$$
w = d\,\hat{\sigma}\sqrt{2/n},
$$

then clipped at 0. Each factor has a job:

- $\sqrt{2/n}$ is the standard error of a difference between two arm means, in units of $\sigma$. It is the factor you control, and it shrinks like $1/\sqrt{n}$.
- $\hat{\sigma}$ is the standard deviation estimated from the test data. Before the test it is unknown: it will land near the true $\sigma$, a little above or a little below.
- $d$ is the MCB critical value. It depends on $k$, $\alpha$ and the error degrees of freedom $\nu = k(n-1)$, and it shrinks slightly as $n$ grows. For $k = 5$ and $\alpha = 0.05$ it is about 2.16 once $n$ is large.

## Why $w < \delta^*$ is the right bar

If the intervals are correct and $w < \delta^*$, both promises hold. Here is why.

The computation works with one specific event from Hsu's proof that the MCB intervals are correct (Theorem 4.1.1). It has probability exactly $1 - \alpha$, and whenever it happens, every interval covers its target. Call it the correctness event. Let $b$ be the true best arm. The event says that no arm looks better against arm $b$ than it really is by $w$ or more:

$$
(\hat{\mu}_i - \hat{\mu}_b) - (\mu_i - \mu_b) < w \quad \text{for every arm } i \neq b.
$$

**Clearly worse arms are ruled out.** Arm $i$'s upper bound is $(\hat{\mu}_i - \max_{j \neq i}\hat{\mu}_j + w)^+$. The best of the other sample means is at least $\hat{\mu}_b$, so

$$
\text{upper bound}_i
\le (\hat{\mu}_i - \hat{\mu}_b + w)^+
< (\mu_i - \mu_b + 2w)^+.
$$

If arm $i$ is more than $2\delta^*$ below the best, then $\mu_i - \mu_b < -2\delta^*$. With $2w < 2\delta^*$, the inside is negative, the upper bound is 0, and arm $i$ is ruled out. The factor of 2 appears because the error is paid twice: once in the estimated difference, and once in the margin $w$ added to it.

**The winner is good enough, and the intervals show it.** The observed winner's sample mean is at least every other arm's, so each estimated difference in its lower bound is at least 0, and its lower bound is at least $-w > -\delta^*$. Its interval is correct, so its true $\theta$ is above $-\delta^*$ as well: the winner is within $\delta^*$ of the true best.

So the plan is to choose $n$ so that, with high probability, the intervals are correct **and** $w < \delta^*$. Hsu calls this the probability of *correct and useful* inference. `mcb` computes it, and finds the smallest $n$ at which it reaches a probability you choose, the `target`.

## The ceiling: `target` must be below $1 - \alpha$

The correctness event has probability $1 - \alpha$. Also requiring $w < \delta^*$ can only lower that, so the probability of correct and useful inference is always below $1 - \alpha$. It approaches $1 - \alpha$ as $n$ grows. `target` must therefore be below $1 - \alpha$, and the functions raise `ValueError` otherwise. Hsu's example uses $\alpha = 0.05$ with a target of $0.90$.

## Why not just plug in $\sigma$?

If $\hat{\sigma}$ were sure to equal $\sigma$, you could solve $w < \delta^*$ directly:

$$
n \ge 2\left(\frac{d\,\sigma}{\delta^*}\right)^2.
$$

For the running example ($d \approx 2.16$, $\sigma = 40$, $\delta^* = 2$), this gives $n = 3735$. But $\hat{\sigma}$ is random. At this $n$, $w$ is just under $\delta^*$ when $\hat{\sigma} = \sigma$, so the test is useful only when $\hat{\sigma}$ comes in at or below $\sigma$, which happens about half the time. The probability of correct and useful inference at $n = 3735$ is 0.48, a coin flip.

The exact computation accounts for the randomness in $\hat{\sigma}$. It asks for $n = 3797$, 62 more users per arm, which lifts the probability to 0.90.

## How the probability grows with $n$

`constrained_mcb_correct_useful_probability` gives the probability at any $n$:

```python
from mcb import constrained_mcb_correct_useful_probability

for n in [3650, 3700, 3750, 3797, 3850, 3900]:
    print(n, constrained_mcb_correct_useful_probability(k=5, n=n, delta_star=2, sigma=40, alpha=0.05))
```

| $n$ per arm | 3650 | 3700 | 3750 | 3797 | 3850 | 3900 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Probability | 0.014 | 0.179 | 0.626 | 0.901 | 0.949 | 0.950 |

The probability climbs from near 0 to its ceiling of 0.95 over a 7% change in $n$. Compare a small experiment, Hsu's own example (Hsu 1996, Appendix C.4, taken from Scheffé): comparing the tensile strength of eight alloys ($k = 8$), with $\sigma$ about 3000 psi, $\delta^* = 5000$ psi and $\alpha = 0.05$:

| $n$ per arm | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Probability | 0.039 | 0.122 | 0.364 | 0.710 | 0.909 | 0.948 | 0.950 |

Here the climb spans a factor of 4 in $n$. The plug-in formula gives $n = 5$, where the probability is only 0.71. The exact answer is $n = 6$, matching Hsu's Figure C.1.

The two curves differ because of how well $\hat{\sigma}$ estimates $\sigma$. With $\nu = k(n-1)$ degrees of freedom, $\hat{\sigma}$ is typically off by about $1/\sqrt{2\nu}$ of $\sigma$.

- In the alloy test at $n = 6$, $\nu = 40$, so $\hat{\sigma}$ is typically off by about 11%. The plan needs room for $\hat{\sigma}$ coming in high, so $n$ goes well above the plug-in value, and the choice of `target` matters.
- In the A/B test at $n \approx 3800$, $\nu \approx 19{,}000$, so $\hat{\sigma}$ is typically off by about 0.5%. For practical purposes $\hat{\sigma} = \sigma$, and the plug-in formula is within 2% of the exact answer. `target` barely moves $n$: a target of 0.80 needs 3773 per arm, and 0.90 needs 3797.

At A/B scale, then, $n$ is set by $\delta^*$, $\sigma$, $k$ and $\alpha$, roughly through $n \approx 2(d\sigma/\delta^*)^2$.

## What if your guess of $\sigma$ is wrong?

At A/B scale, the noise in $\hat{\sigma}$ is tiny, so the risk moves to your planning value. If the true $\sigma$ is 41 minutes rather than 40, $\hat{\sigma}$ comes in near 41, $w$ comes in near 2.03, and the probability that $w < 2$ is essentially 0 (it is 0.001).

That sounds worse than it is. The intervals are still correct with probability at least $1 - \alpha$; they are only a little wider than planned. The margin you can claim scales with the ratio of the two values: about $\delta^* \times (\text{true } \sigma / \text{planning } \sigma)$. Plan with $\sigma = 40$, get 42, and the winner is shown to be within about 2.1 minutes of the best instead of 2.

So:

- Use a planning $\sigma$ you are fairly sure is not too low, such as an upper estimate from past data rather than the average.
- After the test, read the intervals you got. They show the margin the data actually support.

## Plan the test

For the running example:

```python
from mcb import constrained_mcb_correct_useful_probability, constrained_mcb_sample_size

n = constrained_mcb_sample_size(k=5, delta_star=2, sigma=40, alpha=0.05, target=0.90)
print(n)  # 3797
print(constrained_mcb_correct_useful_probability(k=5, n=n, delta_star=2, sigma=40, alpha=0.05))  # about 0.901
```

That is 3797 users per arm, or 18,985 in total. At 2,500 new users a day, the test runs for about 8 days. What it buys, with probability at least 0.90:

- every algorithm more than 4 minutes below the best is ruled out;
- the observed winner is shown to be within 2 minutes of the best.

### How the answer responds to each input

Each row changes one input from the running example ($k = 5$, $\delta^* = 2$, $\sigma = 40$, $\alpha = 0.05$, target 0.90):

| Change | $n$ per arm | Total users |
| --- | ---: | ---: |
| None (running example) | 3797 | 18,985 |
| $\delta^* = 1$ | 15,061 | 75,305 |
| $\delta^* = 4$ | 966 | 4,830 |
| $\sigma = 50$ | 5913 | 29,565 |
| $k = 2$ | 2241 | 4,482 |
| $k = 10$ | 4724 | 47,240 |
| $k = 20$ | 5576 | 111,520 |
| $\alpha = 0.10$, target 0.85 | 2757 | 13,785 |
| Target 0.80 | 3773 | 18,865 |

Three things to take from the table:

- $n$ scales with $(\sigma/\delta^*)^2$. Halve $\delta^*$ and you need four times the users. $\delta^*$ is the input worth arguing about.
- Per-arm $n$ grows slowly with $k$, because $k$ enters only through $d$. Total traffic, $k \cdot n$, grows faster.
- $\alpha$ moves $n$ much more than `target` does, because $\alpha$ sets $d$.

### Constrained MCB or Edwards-Hsu?

If you will analyze with `unconstrained_mcb_edwards_hsu`, plan with `unconstrained_mcb_edwards_hsu_sample_size`. It makes the same two promises, with the screening set $S$ in place of $G$. Its intervals use the larger critical value $|d|$, so it needs more users: 4842 per arm in the running example, against 3797. In return it gives two-sided bounds on how far each arm is below the best.

Plan with the procedure and the $\alpha$ you will use in the analysis. The guarantee is about that analysis.

## Check the plan by simulation

The 0.90 is a probability over repeated tests. You can watch it happen: make up a truth, simulate the test many times, analyze each one with `constrained_mcb`, and count how often both promises come true.

```python
import numpy as np
from mcb import constrained_mcb, summarize_from_samples

rng = np.random.default_rng(1)
true_means = {"Control": 50.0, "A": 50.0, "B": 51.0, "C": 53.5, "D": 54.5}
best = max(true_means.values())
n, sigma, delta_star = 3797, 40.0, 2.0

reps, successes = 1000, 0
for _ in range(reps):
    samples = {arm: rng.normal(mu, sigma, n) for arm, mu in true_means.items()}
    result = constrained_mcb(summarize_from_samples(samples), alpha=0.05)
    winner = max(samples, key=lambda arm: samples[arm].mean())
    lower = result.lower[result.labels.index(winner)]
    ruled_out = not any(true_means[arm] < best - 2 * delta_star for arm in result.candidate_best_set)
    near_best = lower > -delta_star and true_means[winner] > best - delta_star
    successes += ruled_out and near_best
print(successes / reps)  # 0.997
```

Here `Control` and `A` are clearly worse (more than 4 minutes below the best), `C` and `D` are good enough, and `B` is in the gray zone, so the check asks nothing of it. Both promises came true in 997 of 1000 simulated tests. The 0.90 is a floor that holds whatever the true means are. This truth has well-separated arms, so it does better than the floor. Try your own guesses at the truth to see what your plan delivers in the cases you care about.

## Choosing the inputs

- `k`: the number of arms, including the control.
- `delta_star`: the shortfall from the best you would accept in the arm you ship. This is a product decision, not a statistical one. Arms more than $2\delta^*$ below the best are guaranteed to be ruled out.
- `sigma`: a planning value for the standard deviation of the per-user metric, in the same units as `delta_star`. Get it from historical data, an A/A test or a pre-period. For a conversion rate near $p$, $\sigma \approx \sqrt{p(1-p)}$. If unsure, round up. Only the ratio $\delta^*/\sigma$ matters.
- `alpha`: the familywise error rate you will use in the analysis.
- `target`: the probability of correct and useful inference you want. It must be below $1 - \alpha$.

The guarantee rests on the one-way normal model: independent users, a common variance across arms, and normal errors, which at A/B scale means approximately normal arm means. See [Assumptions and scope](assumptions.md).

## Why not a power calculation?

The usual way to size an experiment is the power of a test, such as the $F$-test that all arms have equal means. Power is the probability that the test rejects. That includes outcomes where the data mislead and the follow-up comparisons make wrong claims, such as $\mu_i > \mu_j$ when in fact $\mu_i < \mu_j$. Power counts those outcomes as successes. Sizing for narrow intervals alone has the same flaw: it counts narrow intervals that miss the truth. Hsu's criterion counts a test as a success only if its intervals are correct and narrow.

Power also answers a different question: "will we detect that the arms differ?" A test with a winner to ship needs to answer "will we find an arm close to the best, and drop the clearly worse ones?" In Hsu's alloy example, $F$-test power asks for 5 specimens per alloy. Correct and useful MCB needs 6.

## How the probability is computed

This section is for readers who want the integrals. Let $S = \hat{\sigma}/\sigma$, with density $\gamma_\nu$ for $\nu = k(n-1)$, and let

$$
u = \sqrt{n/2}\,\frac{\delta^*/\sigma}{d}.
$$

The narrowness condition $w < \delta^*$ is the same as $S < u$. The critical values are the balanced ones, with every $\lambda = 1/\sqrt{2}$.

`constrained_mcb_correct_useful_probability` computes Hsu's Eq. (C.3):

$$
\int_0^u \int_{-\infty}^{\infty} \left[\Phi(z + \sqrt{2}\,d\,s)\right]^{k-1} d\Phi(z)\,\gamma_\nu(s)\,ds.
$$

Read it from the inside out:

- $z$ is the standardized error in the true best arm's sample mean, $\sqrt{n}(\hat{\mu}_b - \mu_b)/\sigma$. Given $z$ and $S = s$, the correctness event needs each of the other $k - 1$ standardized errors to stay below $z + \sqrt{2}\,d\,s$. They are independent, so that has probability $\Phi(z + \sqrt{2}\,d\,s)^{k-1}$.
- Averaging over $z$ gives the probability of the correctness event when $S = s$.
- Averaging over all $s > 0$ would give exactly $1 - \alpha$, because that is how $d$ is defined. Stopping the outer integral at $u$ keeps only the values of $\hat{\sigma}$ that make $w < \delta^*$.

`unconstrained_mcb_edwards_hsu_correct_useful_probability` computes Eq. (C.5), the same integral with $|d|$ in place of $d$ (including in $u$) and the two-sided condition that each standardized error stays within $\sqrt{2}\,|d|\,s$ of $z$:

$$
\int_0^u \int_{-\infty}^{\infty} \left[\Phi(z + \sqrt{2}\,|d|\,s) - \Phi(z - \sqrt{2}\,|d|\,s)\right]^{k-1} d\Phi(z)\,\gamma_\nu(s)\,ds.
$$

For constrained MCB, the two promises are Hsu's statements in Appendix C.2. For Edwards-Hsu, Hsu says only that the computation is "essentially the same" as for two-sided MCC (Appendix C.3); both promises follow from the correctness event in the proof of Theorem 4.2.3, by the same argument as above.

The correctness event implies that every interval covers its target, but the intervals can also cover without it. So both values are conservative: they are lower bounds on the probability that the intervals cover and are narrow.

The `*_sample_size` functions return the smallest $n \ge 2$ at which the matching probability reaches `target`.

## Scope and numerical details

- The computation covers balanced designs only, as Appendix C does: one common $n$ for every arm. The analysis functions (`constrained_mcb`, `unconstrained_mcb_edwards_hsu`) still accept unequal sample sizes.
- The ceiling of $1 - \alpha$ holds up to the critical-value solver's tolerance of about $2 \times 10^{-5}$.
- The critical values switch to their infinite-degrees-of-freedom form at $\nu = 10^6$ (see {doc}`mcb_critical_values`). Above the switch, the narrowness condition still uses the exact chi-square distribution of $S$, while the correctness part is evaluated at $\hat{\sigma} = \sigma$, which changes the probability by about $10^{-4}$. The critical value just below the switch differs from the infinite-df one by slightly more than the true finite-df correction, within the solver's tolerance. As a result, for $k \ge 3$ the probability jumps by about $10^{-3}$ across the switch. It still increases with $n$, and a sample size near the switch is off by a few units out of $10^5$ to $10^6$.
