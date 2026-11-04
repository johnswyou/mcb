# Sample Size

`mcb` can choose the number of users per arm for an A/B/n test that will be analyzed with MCB. It finds the smallest per-arm sample size for which, with a chosen probability, the MCB intervals at the end of the test are both correct and narrow enough to be useful: they identify an arm close to the best and rule out the arms that are clearly worse. The method is Hsu's (Hsu 1988; Hsu 1996, Appendix C).

## Running example

A streaming service plans to test four new recommendation algorithms against the current one, so there are $k = 5$ arms. Users are assigned to the arms evenly and at random, and about 2,500 new users arrive each day. The metric is minutes watched per user in the first week, with larger values better. Last month's data put its standard deviation at about 40 minutes.

The results will be analyzed with `constrained_mcb`. For each arm it gives a confidence interval for

$$
\theta_i = \mu_i - \max_{j \neq i} \mu_j,
$$

the difference between arm $i$'s true mean and the largest of the other arms' true means. The $k$ intervals hold simultaneously with probability at least $1 - \alpha$. An arm whose upper bound is 0 is ruled out as the best, and the arms that remain form the candidate-best set $G$ (see [Constrained MCB](constrained_mcb.md)).

The planning question is how many users each arm needs.

## The answer in brief

```python
from mcb import constrained_mcb_sample_size

n = constrained_mcb_sample_size(k=5, delta_star=2, sigma=40, alpha=0.05, target=0.90)
print(n)  # 3797
```

The test needs 3797 users per arm, 18,985 in total, or about 8 days of traffic. With probability at least 0.90, the algorithm with the highest sample mean is shown to be within 2 minutes (`delta_star`) of the true best, and every algorithm more than 4 minutes below the best is ruled out.

For a test this large, a simple formula comes within 2% of that answer:

$$
n \approx 2\left(\frac{d\,\sigma}{\delta^*}\right)^2,
$$

where $d$ is the MCB critical value, about 2.16 for $k = 5$ and $\alpha = 0.05$. The input that most needs care is the planning value of $\sigma$. If the true standard deviation turns out larger, the intervals are still valid, but the margin they can show grows in proportion.

The rest of the page explains what the guarantee means, why the exact answer is larger than the formula, how to choose the inputs and how to check the plan by simulation. The derivations are collected at the end.

## What the guarantee means

At the end of the test the team wants to draw two conclusions:

1. The arm with the highest sample mean (the *observed winner*) is not much worse than the true best arm, so it can be shipped.
2. The arms that are clearly worse than the best are ruled out.

Both conclusions need a tolerance. Identifying the exact best arm is not a realistic goal: two algorithms that differ by 0.01 minutes cannot be separated by a test of practical size, and the difference would not matter. The tolerance is a margin $\delta^* > 0$ (`delta_star`), in the units of the metric. A shortfall from the best smaller than $\delta^*$ is treated as negligible. The example uses $\delta^* = 2$ minutes.

The margin divides the arms into three zones according to how far each arm's true mean is below the best:

```text
       clearly worse          gray zone       good enough
  <--------------------|-----------------|---------------|
                  best - 2δ*         best - δ*         best
```

With the sample size computed here, the following holds with at least the chosen probability: the observed winner is in the good-enough zone and its interval shows this, and every arm in the clearly-worse zone is ruled out. Arms in the gray zone may or may not be ruled out. The next section explains why the outer boundary is at $2\delta^*$.

Readers used to power calculations may look for a minimum detectable effect. The closest counterpart is $2\delta^*$: an arm that far below the best is guaranteed to be ruled out. $\delta^*$ itself is a tolerance on the arm that ships, not an effect to detect.

## What the sample size controls

Two things can prevent these conclusions. First, an interval can miss its true value. MCB keeps the probability of this at or below $\alpha$ at every sample size, so adding users does not change it. Second, the intervals can all be correct but too wide to rule any arm out or to show that the winner is close to the best. Width is what the sample size controls.

In a balanced test with $n$ users per arm, every constrained MCB bound is an estimated difference between two arm means, plus or minus the same amount

$$
w = d\,\hat{\sigma}\sqrt{2/n},
$$

clipped at 0. It has three factors:

- $\sqrt{2/n}$ is the standard error of a difference between two arm means, in units of $\sigma$. It is the only factor set directly by the design, and it decreases like $1/\sqrt{n}$.
- $\hat{\sigma}$ is the pooled standard deviation estimated from the test data. It is random. Before the test, all that is known is that it will be close to $\sigma$.
- $d$ is the MCB critical value. It depends on $k$, $\alpha$ and the error degrees of freedom $\nu = k(n-1)$, and it decreases slightly as $n$ grows. For $k = 5$ and $\alpha = 0.05$ it is about 2.16 when $n$ is large.

The margin and $w$ do different jobs. $\delta^*$ and $2\delta^*$ are lines drawn before the test, measured from the true best arm. $w$ is a ruler that the test supplies, measured back from the highest sample mean. The position it reaches is called the *cut line* here: an arm whose sample mean is at or behind the cut line has an upper bound of 0 and is ruled out.

The figures in this section use the true means from [Checking the plan by simulation](#checking-the-plan-by-simulation). With $\delta^* = 2$, `Control` and `A` are clearly worse, `B` is in the gray zone, and `C` and `D` are good enough.

```{figure} figures/sample_size_one_test.svg
:alt: Five arms on an axis of minutes watched, over three zones measured from the true best of 54.5. Rings show true means and dots show sample means. C leads at 53.54 and the cut line sits w = 1.98 behind it at 51.57. Control and A are behind the line and ruled out. B, C and D are kept.
:width: 100%

One simulated test with 3797 users per arm. This test was chosen because it shows the harder case, in which the true best arm does not lead.
```

`D` is the true best, but `C` has the highest sample mean. That is acceptable, because `C`'s true mean is in the good-enough zone. Both clearly worse arms are behind the cut line and are ruled out. `B` is in the gray zone, and here noise carries its sample mean past the line, so it stays in.

More users shorten the ruler. The next figure repeats this test at three sample sizes, with its random errors scaled to each $n$ and $\hat{\sigma}$ held fixed, so $n$ is the only thing that changes. The bars under each panel compare $w$ with $\delta^*$.

```{figure} figures/sample_size_more_users.svg
:alt: Three panels of the same test at 1000, 3797 and 15,000 users per arm. As n grows, sample means move toward true means and w falls from 3.86 to 1.98 to 1.00, compared with a fixed delta-star of 2. At 1000, A is clearly worse but stays in. At 3797, Control and A are out. At 15,000, gray-zone arm B is out too.
:width: 100%

The same test at three sample sizes.
```

With 1000 users per arm, $w$ is wider than $\delta^*$, and `A` stays in even though it is clearly worse. With 3797, $w$ just fits inside $\delta^*$, which is what the plan buys. With 15,000, $w$ is about half of $\delta^*$, and `B` is ruled out too. Because $w$ decreases like $1/\sqrt{n}$, halving it takes about four times the users.

If the intervals are correct and $w < \delta^*$, both conclusions follow:

- Every arm more than $2\delta^*$ below the best is ruled out. To be sure of that, an arm's true gap has to cover two distances. The arm must end up at least $w$ behind the leader to fall behind the cut line. It also needs a further $w$ to absorb luck, because when the intervals are correct, noise can make its gap to the true best look smaller than it really is by less than $w$. A true gap of more than $2\delta^*$ is more than $2w$, so it covers both. The leader is never behind the true best arm, so this holds whichever arm leads.
- The observed winner trails no other arm, so its lower bound is at least $-w$, which is above $-\delta^*$. Because its interval is correct, its true shortfall from the best is less than $\delta^*$.

In the running example, $w$ is about 2 minutes. `A` is truly 4.5 minutes behind the best, so even if luck hides nearly 2 of them, it still looks about 2.5 behind and is ruled out. `B` is truly 3.5 behind, so if luck hides nearly 2, it looks only about 1.5 behind and can stay in. This is why the gray zone exists: an arm between $\delta^*$ and $2\delta^*$ below the best may or may not be ruled out, so the plan promises nothing about it. $2\delta^*$ is where the guarantee starts, not where arms start to be ruled out.

The next figure shows 16 simulated tests at the planned sample size. They are consecutive draws, not selected, sorted by the leader's sample mean, so the cut lines form one edge and the leaders form a second edge about $w$ to its right.

```{figure} figures/sample_size_many_tests.svg
:alt: Sixteen rows, one per simulated test, sorted by the leader's sample mean. In each row a gray bar runs from the left edge to that test's cut line, so the bar ends form one stepped edge, and the leaders form a second edge about 2 minutes to its right. Control and A sit on the gray bar in every row. B sits on it in all but two rows. D, a dark blue square, usually leads, and C, a light blue circle, leads in two rows.
:width: 100%

Sixteen simulated tests with 3797 users per arm. Each gray bar ends at that test's cut line.
```

The cut line moves with whichever arm leads, and `C` leads in 2 of these tests. `Control` and `A` are behind the line in all 16. `B` usually is, but it stays in twice, which is what the gray zone means. [Checking the plan by simulation](#checking-the-plan-by-simulation) repeats this over 1000 tests.

Hsu calls the probability that the intervals are correct and $w < \delta^*$ the probability of *correct and useful* inference. `mcb` computes it, and the sample-size functions find the smallest $n$ at which it reaches `target`. Strictly, the computation uses an event, the *correctness event* (defined in [Derivations](#derivations)), that has probability exactly $1 - \alpha$ and guarantees that every interval is correct. Intervals can also be correct without it, so the computed probability is a lower bound and the plan is conservative.

Because the correctness event has probability $1 - \alpha$, the probability of correct and useful inference is always below $1 - \alpha$ and approaches it as $n$ grows. `target` must therefore be below $1 - \alpha$; otherwise the functions raise `ValueError`. With $\alpha = 0.05$, Hsu's example uses a target of 0.90.

## Why not a power calculation

Experiments are often sized by power: the probability that a test, such as the $F$-test that all arms have equal means, rejects. Power answers a different question, whether the test will detect that the arms differ. The question here is whether the test will identify an arm close to the best and rule out the clearly worse ones.

Power also counts some failures as successes. A test can reject on misleading data, after which follow-up comparisons claim $\mu_i > \mu_j$ when in fact $\mu_i < \mu_j$. Sizing for narrow intervals alone has the same problem, because it counts narrow intervals that miss the truth. Hsu's criterion counts a test as a success only when its intervals are correct and narrow.

The two approaches give different sizes. For the running example, an $F$-test with power 0.90 whenever the best and worst arms differ by at least 4 minutes ($2\delta^*$) needs 3082 users per arm. The MCB plan needs 3797, and the extra users buy an answer to which arm to ship, which the $F$-test does not give.

## The known-$\sigma$ formula and its limits

### Why the exact answer is larger

If $\hat{\sigma}$ were equal to $\sigma$, the condition $w < \delta^*$ could be solved for $n$ directly, which gives the formula from [the answer in brief](#the-answer-in-brief):

$$
n \ge 2\left(\frac{d\,\sigma}{\delta^*}\right)^2.
$$

For the running example ($d \approx 2.16$, $\sigma = 40$, $\delta^* = 2$) this gives $n = 3735$. But $\hat{\sigma}$ is random. At this $n$, $w < \delta^*$ holds only when $\hat{\sigma}$ comes out at or below roughly $\sigma$, which happens about half the time. Combined with the 0.95 probability of the correctness event, the probability of correct and useful inference at $n = 3735$ is 0.48.

The exact computation accounts for the distribution of $\hat{\sigma}$ and gives $n = 3797$. The extra 62 users per arm raise the probability to 0.90.

### Large tests and small tests

`constrained_mcb_correct_useful_probability` gives the probability at any $n$:

```python
from mcb import constrained_mcb_correct_useful_probability

for n in [3650, 3700, 3750, 3797, 3850, 3900]:
    print(n, constrained_mcb_correct_useful_probability(k=5, n=n, delta_star=2, sigma=40, alpha=0.05))
```

| $n$ per arm | 3650 | 3700 | 3750 | 3797 | 3850 | 3900 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Probability | 0.014 | 0.179 | 0.626 | 0.901 | 0.949 | 0.950 |

Over this 7% range of $n$, the probability goes from near 0 to its limit of 0.95. Small experiments behave differently. Hsu's example (1996, Appendix C.4, taken from Scheffé) compares the tensile strength of eight alloys ($k = 8$), with $\sigma$ about 3000 psi, $\delta^* = 5000$ psi and $\alpha = 0.05$:

| $n$ per arm | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Probability | 0.039 | 0.122 | 0.364 | 0.710 | 0.909 | 0.948 | 0.950 |

Here the probability rises over a fourfold range of $n$. The known-$\sigma$ formula gives $n = 5$, where the probability is only 0.71. The exact answer is $n = 6$, which agrees with Hsu's Figure C.1. (An $F$-test power calculation, for differences of 10,000 psi with probability 0.90, gives 5.)

The difference between the two cases is how precisely $\hat{\sigma}$ is estimated. With $\nu = k(n-1)$ degrees of freedom, the standard deviation of $\hat{\sigma}/\sigma$ is about $1/\sqrt{2\nu}$.

- In the alloy example at $n = 6$, $\nu = 40$ and this is about 11%. The sample size has to allow for $\hat{\sigma}$ coming out well above $\sigma$, so the exact $n$ is larger than the known-$\sigma$ value and depends on `target`.
- In the A/B example at $n \approx 3800$, $\nu \approx 19{,}000$ and it is about 0.5%. $\hat{\sigma}$ is then close to $\sigma$, the known-$\sigma$ formula is within 2% of the exact answer, and `target` has little effect: 0.80 gives 3773 per arm and 0.90 gives 3797.

For tests of A/B size, then, $n$ is set mainly by $\delta^*$, $\sigma$, $k$ and $\alpha$ through the known-$\sigma$ formula.

### When the planning $\sigma$ is wrong

For large tests, the planning value of $\sigma$ is a bigger source of error than the sampling variation of $\hat{\sigma}$. If the true $\sigma$ is 41 minutes instead of 40, $\hat{\sigma}$ will be close to 41 and $w$ close to 2.03, and the probability that $w < 2$ is 0.001.

The intervals are still correct with probability at least $1 - \alpha$, and they are only slightly wider than planned. In this regime the margin achieved is approximately $\delta^* \times (\text{true } \sigma / \text{planning } \sigma)$. With a planning value of 40 and a true $\sigma$ of 42, the winner is shown to be within about 2.1 minutes of the best instead of 2.

It is therefore better to use a planning $\sigma$ that is unlikely to be too low, such as an upper estimate from past data rather than the average. After the test, the intervals show the margin that the data support.

## Choosing the inputs

- `k`: the number of arms, including the control.
- `delta_star`: the largest shortfall from the best that is acceptable in the arm that ships. It is set by what difference matters for the product. Arms more than $2\delta^*$ below the best are guaranteed to be ruled out.
- `sigma`: a planning value for the standard deviation of the per-user metric, in the same units as `delta_star`. Historical data, an A/A test or a pre-period can supply it. For a conversion rate near $p$, $\sigma \approx \sqrt{p(1-p)}$. When in doubt, use a larger value. Only the ratio $\delta^*/\sigma$ enters the calculation.
- `alpha`: the familywise error rate that will be used in the analysis.
- `target`: the required probability of correct and useful inference. It must be below $1 - \alpha$.

The guarantee rests on the one-way normal model: independent users, a common variance across arms, and normal errors. At A/B scale, approximately normal arm means are usually an adequate substitute for normal errors. See [Assumptions and scope](assumptions.md).

### Sensitivity to the inputs

Each row changes one input from the running example ($k = 5$, $\delta^* = 2$, $\sigma = 40$, $\alpha = 0.05$, target 0.90). The $\alpha = 0.10$ row also lowers the target to 0.85, because the target must be below $1 - \alpha = 0.90$.

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

$n$ is roughly proportional to $(\sigma/\delta^*)^2$, so halving $\delta^*$ multiplies it by about four. Per-arm $n$ grows slowly with $k$, because $k$ enters only through $d$, while total traffic $k \cdot n$ grows faster. $\alpha$ changes $n$ much more than `target` does, because $\alpha$ determines $d$.

### Planning for Edwards-Hsu MCB

For an analysis with `unconstrained_mcb_edwards_hsu`, use `unconstrained_mcb_edwards_hsu_sample_size`, and `unconstrained_mcb_edwards_hsu_correct_useful_probability` for the probability. The same two conclusions hold, with the Edwards-Hsu screening set in place of $G$. Edwards-Hsu intervals use a two-sided critical value, which Hsu writes $|d|$. It is a separate constant, not the absolute value of $d$, and it is larger, so more users are needed: 4842 per arm in the running example, compared with 3797. In exchange, the analysis gives two-sided bounds on how far each arm is below the best (see [Edwards-Hsu Unconstrained MCB](edwards_hsu.md)).

Plan with the procedure and $\alpha$ that the analysis will use, since the guarantee applies to that analysis.

## Checking the plan by simulation

The 0.90 is a probability over repeated tests, so it can be checked by simulation. Two checks are useful.

The first reproduces the 0.90 itself. The correctness event depends only on the estimation errors $e_i = \hat{\mu}_i - \mu_i$ and on $\hat{\sigma}$, not on the true means. With arm 0 as the true best, it requires $e_i - e_0 < w$ for every other arm $i$. The simulation can therefore draw the errors and $\hat{\sigma}$ directly:

```python
import numpy as np
from mcb import constrained_critical_value

rng = np.random.default_rng(1)
k, n, sigma, delta_star = 5, 3797, 40.0, 2.0
nu = k * (n - 1)
d = constrained_critical_value((1 / np.sqrt(2),) * (k - 1), nu=nu, alpha=0.05)  # balanced design

reps = 100_000
errors = rng.normal(0, sigma / np.sqrt(n), size=(reps, k))  # arm 0 is the true best
sigma_hat = sigma * np.sqrt(rng.chisquare(nu, reps) / nu)
w = d * sigma_hat * np.sqrt(2 / n)
correct = np.all(errors[:, 1:] - errors[:, [0]] < w[:, None], axis=1)
print(np.mean(correct & (w < delta_star)))  # 0.90155
```

The second fixes a set of true means, simulates the whole test, analyzes each simulated test with `constrained_mcb`, and counts how often both conclusions are drawn and correct. Expect a number above 0.90: the 0.90 is a lower bound that applies to every configuration of true means, and the conclusions can also hold when the correctness event fails.

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

In this configuration, `Control` and `A` are more than 4 minutes below the best, `C` and `D` are within 2 minutes of it, and `B` is in the gray zone, so the check places no requirement on `B`. Both conclusions held in 997 of the 1000 simulated tests. Other configurations can be checked by changing `true_means`.

## Derivations

### Why the two conclusions follow

The computation uses an event from the proof of Hsu's Theorem 4.1.1, called the correctness event here. It has probability exactly $1 - \alpha$, and when it occurs every interval covers its target. With $b$ the true best arm, the event is

$$
(\hat{\mu}_i - \hat{\mu}_b) - (\mu_i - \mu_b) < w \quad \text{for every arm } i \neq b,
$$

that is, the estimated difference between any arm and arm $b$ overstates the true difference by less than $w$.

For an arm $i \neq b$, the upper bound is $(\hat{\mu}_i - \max_{j \neq i}\hat{\mu}_j + w)^+$, where $x^+ = \max(x, 0)$. The largest of the other sample means is at least $\hat{\mu}_b$, so when the correctness event occurs,

$$
\text{upper bound}_i
\le (\hat{\mu}_i - \hat{\mu}_b + w)^+
\le (\mu_i - \mu_b + 2w)^+.
$$

If arm $i$ is more than $2\delta^*$ below the best and $w < \delta^*$, then $\mu_i - \mu_b < -2\delta^* < -2w$, the expression inside is negative, and the upper bound is 0. Arm $i$ is ruled out.

The observed winner's sample mean is at least as large as every other arm's, so each estimated difference in its lower bound is at least 0, and the lower bound is at least $-w > -\delta^*$. Because its interval is correct, the winner's true $\theta$ is also above $-\delta^*$, which means it is within $\delta^*$ of the true best.

The correctness event implies that every interval covers its target, but the intervals can also cover their targets when it does not occur. The computed probabilities are therefore lower bounds on the probability that the intervals cover their targets and are narrow.

### How the probability is computed

Let $S = \hat{\sigma}/\sigma$. Then $\nu S^2$ has a chi-square distribution with $\nu = k(n-1)$ degrees of freedom, and $\gamma_\nu$ denotes the density of $S$. Let

$$
u = \sqrt{n/2}\,\frac{\delta^*/\sigma}{d}.
$$

The condition $w < \delta^*$ is equivalent to $S < u$. In a balanced design every $\lambda_{ij}$ equals $1/\sqrt{2}$, so every arm has the same critical value $d$ (see [MCB Critical Values](mcb_critical_values.md)).

`constrained_mcb_correct_useful_probability` computes Hsu's Eq. (C.3):

$$
\int_0^u \int_{-\infty}^{\infty} \left[\Phi(z + \sqrt{2}\,d\,s)\right]^{k-1} d\Phi(z)\,\gamma_\nu(s)\,ds.
$$

Working from the inside out:

- $z$ is the standardized error in the true best arm's sample mean, $\sqrt{n}(\hat{\mu}_b - \mu_b)/\sigma$. Given $z$ and $S = s$, the correctness event requires each of the other $k - 1$ standardized errors to be below $z + \sqrt{2}\,d\,s$. The errors are independent, so this has probability $\Phi(z + \sqrt{2}\,d\,s)^{k-1}$.
- Averaging over $z$ gives the probability of the correctness event when $S = s$.
- Averaging over all $s > 0$ would give exactly $1 - \alpha$, because that is how $d$ is defined. Stopping the outer integral at $u$ keeps only the values of $\hat{\sigma}$ for which $w < \delta^*$.

`unconstrained_mcb_edwards_hsu_correct_useful_probability` computes Eq. (C.5). It is the same integral with $|d|$ in place of $d$ (including in $u$) and the two-sided condition that each standardized error is within $\sqrt{2}\,|d|\,s$ of $z$:

$$
\int_0^u \int_{-\infty}^{\infty} \left[\Phi(z + \sqrt{2}\,|d|\,s) - \Phi(z - \sqrt{2}\,|d|\,s)\right]^{k-1} d\Phi(z)\,\gamma_\nu(s)\,ds.
$$

For constrained MCB, the two conclusions are Hsu's statements in Appendix C.2. For Edwards-Hsu, Hsu says only that the computation is "essentially the same" as for two-sided multiple comparisons with a control (MCC, Appendix C.3). Both conclusions follow from the correctness event in the proof of Theorem 4.2.3, by the same argument as above.

The `*_sample_size` functions return the smallest $n \ge 2$ at which the corresponding probability reaches `target`.

## Scope and numerical details

- The computation covers balanced designs only, as Appendix C does: one common $n$ for every arm. The analysis functions (`constrained_mcb`, `unconstrained_mcb_edwards_hsu`) accept unequal sample sizes.
- The upper limit of $1 - \alpha$ holds up to the critical-value solver's tolerance of about $2 \times 10^{-5}$.
- The critical values switch to their infinite-degrees-of-freedom form at $\nu = 10^6$ (see [MCB Critical Values](mcb_critical_values.md)). For $k \ge 3$ the probability jumps by about $10^{-3}$ across the switch, mainly because the critical value just below the switch differs from its infinite-df value by slightly more than the true finite-df correction, within the solver's tolerance. The probability still increases with $n$, and a sample size near the switch can be off by a few units out of $10^5$ to $10^6$.
