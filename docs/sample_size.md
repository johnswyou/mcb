# Sample Size

`mcb` can choose the number of users per arm for an A/B/n test that will be analyzed with MCB. The calculation follows Hsu (1988) and Hsu (1996, Appendix C). It finds the smallest per-arm sample size for which, with a chosen probability, the MCB intervals at the end of the test are correct and narrow enough to identify a near-best arm and rule out the arms that are clearly worse.

The sections below derive the criterion. The function calls are in [Computing the sample size](#computing-the-sample-size).

## Running example

A streaming service plans to test four new recommendation algorithms against the current one, so there are $k = 5$ arms. Users are assigned to the arms evenly and at random. The metric is minutes watched per user in the first week, with larger values better. Last month's data put its standard deviation at about 40 minutes.

The results will be analyzed with `constrained_mcb`, which gives a confidence interval for each arm's

$$
\theta_i = \mu_i - \max_{j \neq i} \mu_j,
$$

the difference between arm $i$'s true mean and the largest of the other arms' true means. The $k$ intervals hold simultaneously with probability at least $1 - \alpha$. An arm whose upper bound is 0 is ruled out as the best, and the remaining arms form the candidate-best set $G$ (see [Constrained MCB](constrained_mcb.md)).

The planning question is how many users each arm needs.

## What the plan guarantees

At the end of the test the team wants to draw two conclusions: that the arm with the highest sample mean (the observed winner) is not much worse than the true best arm, so it can be shipped, and that the arms clearly worse than the best are ruled out.

Both conclusions need a tolerance. Identifying the exact best arm is not a realistic goal: two algorithms that differ by 0.01 minutes cannot be separated by a test of practical size, and the difference would not matter. The tolerance is a margin $\delta^* > 0$, in the units of the metric. A shortfall from the best smaller than $\delta^*$ is treated as negligible. The example uses $\delta^* = 2$ minutes.

The margin divides the arms into three zones according to how far each arm's true mean is below the best:

```text
       clearly worse          gray zone       good enough
  <--------------------|-----------------|---------------|
                  best - 2δ*         best - δ*         best
```

With the sample size computed here, the following holds with at least the chosen probability: the observed winner is in the good-enough zone and its interval shows this, and every arm in the clearly-worse zone is ruled out. Arms in the gray zone may or may not be ruled out. The reason for the factor of 2 is given below.

## Coverage and interval width

Two things can prevent these conclusions. First, an interval can miss its true value. MCB limits the probability of this to at most $\alpha$ at every sample size, so adding users does not change it. Second, the intervals can all be correct but too wide to rule any arm out or to put the winner's lower bound near 0. The width decreases as $n$ increases, and this is what the sample size controls.

In a balanced test with $n$ users per arm, every constrained MCB bound is an estimated difference between two arm means plus or minus the same amount

$$
w = d\,\hat{\sigma}\sqrt{2/n},
$$

clipped at 0. The three factors are:

- $\sqrt{2/n}$, the standard error of a difference between two arm means in units of $\sigma$. It is the only factor set directly by the design, and it decreases like $1/\sqrt{n}$.
- $\hat{\sigma}$, the pooled standard deviation estimated from the test data. It is random. Before the test, all that is known is that it will be close to $\sigma$.
- $d$, the MCB critical value. It depends on $k$, $\alpha$ and the error degrees of freedom $\nu = k(n-1)$, and it decreases slightly as $n$ grows. For $k = 5$ and $\alpha = 0.05$ it is about 2.16 when $n$ is large.

## Consequences of $w < \delta^*$

If the intervals are correct and $w < \delta^*$, both conclusions follow.

The computation uses an event from the proof of Hsu's Theorem 4.1.1, called the correctness event here. It has probability exactly $1 - \alpha$, and when it occurs every interval covers its target. With $b$ the true best arm, the event is

$$
(\hat{\mu}_i - \hat{\mu}_b) - (\mu_i - \mu_b) < w \quad \text{for every arm } i \neq b,
$$

that is, the estimated difference between any arm and arm $b$ overstates the true difference by less than $w$.

For an arm $i \neq b$, the upper bound is $(\hat{\mu}_i - \max_{j \neq i}\hat{\mu}_j + w)^+$. The largest of the other sample means is at least $\hat{\mu}_b$, so

$$
\text{upper bound}_i
\le (\hat{\mu}_i - \hat{\mu}_b + w)^+
< (\mu_i - \mu_b + 2w)^+.
$$

If arm $i$ is more than $2\delta^*$ below the best, then $\mu_i - \mu_b < -2\delta^*$. Since $2w < 2\delta^*$, the expression inside is negative, the upper bound is 0, and arm $i$ is ruled out. The bound picks up $w$ twice, once from the error in the estimated difference and once from the margin added to it, which is where the factor of 2 comes from.

The observed winner's sample mean is at least as large as every other arm's, so each estimated difference in its lower bound is at least 0, and the lower bound is at least $-w > -\delta^*$. Because the interval is correct, the winner's true $\theta$ is also above $-\delta^*$, which means it is within $\delta^*$ of the true best.

The sample size is therefore chosen so that the correctness event and $w < \delta^*$ occur together with high probability. Hsu calls this the probability of *correct and useful* inference. `mcb` computes it and finds the smallest $n$ at which it reaches a specified value, `target`.

## Upper limit on `target`

Because the correctness event has probability $1 - \alpha$, the probability of correct and useful inference is below $1 - \alpha$ for every $n$, and it approaches $1 - \alpha$ as $n$ increases. `target` must therefore be less than $1 - \alpha$; otherwise the functions raise `ValueError`. Hsu's example uses $\alpha = 0.05$ and a target of $0.90$.

## Comparison with the known-$\sigma$ formula

If $\hat{\sigma}$ were equal to $\sigma$, the condition $w < \delta^*$ could be solved for $n$ directly:

$$
n \ge 2\left(\frac{d\,\sigma}{\delta^*}\right)^2.
$$

For the running example ($d \approx 2.16$, $\sigma = 40$, $\delta^* = 2$) this gives $n = 3735$. Because $\hat{\sigma}$ is random, at this $n$ the condition $w < \delta^*$ holds only when $\hat{\sigma}$ comes out at or below roughly $\sigma$, which happens about half the time. The probability of correct and useful inference at $n = 3735$ is 0.48.

The exact computation, which accounts for the distribution of $\hat{\sigma}$, gives $n = 3797$. The extra 62 users per arm raise the probability to 0.90.

## Probability as a function of $n$

`constrained_mcb_correct_useful_probability` gives the probability at any $n$:

```python
from mcb import constrained_mcb_correct_useful_probability

for n in [3650, 3700, 3750, 3797, 3850, 3900]:
    print(n, constrained_mcb_correct_useful_probability(k=5, n=n, delta_star=2, sigma=40, alpha=0.05))
```

| $n$ per arm | 3650 | 3700 | 3750 | 3797 | 3850 | 3900 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Probability | 0.014 | 0.179 | 0.626 | 0.901 | 0.949 | 0.950 |

Over this 7% range of $n$, the probability goes from near 0 to its limit of 0.95. For comparison, Hsu's example (1996, Appendix C.4, taken from Scheffé) is a small experiment comparing the tensile strength of eight alloys ($k = 8$), with $\sigma$ about 3000 psi, $\delta^* = 5000$ psi and $\alpha = 0.05$:

| $n$ per arm | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Probability | 0.039 | 0.122 | 0.364 | 0.710 | 0.909 | 0.948 | 0.950 |

Here the probability rises over a fourfold range of $n$. The known-$\sigma$ formula gives $n = 5$, where the probability is 0.71. The exact answer is $n = 6$, which agrees with Hsu's Figure C.1.

The difference between the two cases is the precision of $\hat{\sigma}$. With $\nu = k(n-1)$ degrees of freedom, the standard deviation of $\hat{\sigma}/\sigma$ is about $1/\sqrt{2\nu}$.

- In the alloy example at $n = 6$, $\nu = 40$ and this is about 11%. The sample size has to allow for $\hat{\sigma}$ coming out well above $\sigma$, so the exact $n$ is larger than the known-$\sigma$ value and depends on `target`.
- In the A/B example at $n \approx 3800$, $\nu \approx 19{,}000$ and it is about 0.5%. $\hat{\sigma}$ is then close to $\sigma$, the known-$\sigma$ formula is within 2% of the exact answer, and `target` has little effect: 0.80 gives 3773 per arm and 0.90 gives 3797.

For tests of this size, $n$ is determined mainly by $\delta^*$, $\sigma$, $k$ and $\alpha$, approximately through $n \approx 2(d\sigma/\delta^*)^2$.

## Effect of an incorrect planning $\sigma$

For large tests, the planning value of $\sigma$ is a bigger source of error than the sampling variation of $\hat{\sigma}$. If the true $\sigma$ is 41 minutes instead of 40, $\hat{\sigma}$ will be close to 41 and $w$ close to 2.03, and the probability that $w < 2$ is 0.001.

The intervals are still correct with probability at least $1 - \alpha$, and they are only slightly wider than planned. In this regime the margin achieved is approximately $\delta^* \times (\text{true } \sigma / \text{planning } \sigma)$. With a planning value of 40 and a true $\sigma$ of 42, the winner is shown to be within about 2.1 minutes of the best instead of 2.

It is therefore better to use a planning $\sigma$ that is unlikely to be too low, such as an upper estimate from past data rather than the average. After the test, the intervals show the margin that the data support.

## Computing the sample size

For the running example:

```python
from mcb import constrained_mcb_correct_useful_probability, constrained_mcb_sample_size

n = constrained_mcb_sample_size(k=5, delta_star=2, sigma=40, alpha=0.05, target=0.90)
print(n)  # 3797
print(constrained_mcb_correct_useful_probability(k=5, n=n, delta_star=2, sigma=40, alpha=0.05))  # about 0.901
```

The test needs 3797 users per arm, 18,985 in total, or about 8 days at 2,500 new users a day. With probability at least 0.90, every algorithm more than 4 minutes below the best is ruled out, and the observed winner is shown to be within 2 minutes of the best.

### Sensitivity to the inputs

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

$n$ is roughly proportional to $(\sigma/\delta^*)^2$, so halving $\delta^*$ multiplies it by about four. Per-arm $n$ grows slowly with $k$, because $k$ enters only through $d$, while total traffic $k \cdot n$ grows faster. $\alpha$ changes $n$ much more than `target` does, because $\alpha$ determines $d$.

### Planning for Edwards-Hsu MCB

For an analysis with `unconstrained_mcb_edwards_hsu`, use `unconstrained_mcb_edwards_hsu_sample_size`. The same two conclusions hold, with the screening set $S$ in place of $G$. Edwards-Hsu intervals use the larger critical value $|d|$, so more users are needed: 4842 per arm in the running example, compared with 3797. In exchange, the analysis gives two-sided bounds on how far each arm is below the best.

The procedure and $\alpha$ used for planning should be the ones used in the analysis, since the guarantee applies to that analysis.

## Checking the plan by simulation

The 0.90 is a probability over repeated tests. It can be checked by simulation: fix a set of true means, simulate the test many times, analyze each simulated test with `constrained_mcb`, and count how often both conclusions are drawn and correct.

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

In this configuration, `Control` and `A` are more than 4 minutes below the best, `C` and `D` are within 2 minutes of it, and `B` is in the gray zone, so the check places no requirement on `B`. Both conclusions held in 997 of the 1000 simulated tests. The 0.90 is a lower bound that applies to every configuration of true means. This configuration has well-separated arms and exceeds it. Other configurations can be checked by changing `true_means`.

## Choosing the inputs

- `k`: the number of arms, including the control.
- `delta_star`: the largest shortfall from the best that is acceptable in the arm that ships. It is set by what difference matters for the product. Arms more than $2\delta^*$ below the best are guaranteed to be ruled out.
- `sigma`: a planning value for the standard deviation of the per-user metric, in the same units as `delta_star`. Historical data, an A/A test or a pre-period can supply it. For a conversion rate near $p$, $\sigma \approx \sqrt{p(1-p)}$. When in doubt, use a larger value. Only the ratio $\delta^*/\sigma$ enters the calculation.
- `alpha`: the familywise error rate that will be used in the analysis.
- `target`: the required probability of correct and useful inference. It must be below $1 - \alpha$.

The guarantee rests on the one-way normal model: independent users, a common variance across arms, and normal errors. At A/B scale, approximately normal arm means are usually an adequate substitute for normal errors. See [Assumptions and scope](assumptions.md).

## Comparison with power calculations

Experiments are often sized by the power of a test, such as the $F$-test that all arms have equal means. Power is the probability that the test rejects. That probability includes outcomes in which the data mislead and the follow-up comparisons make wrong claims, such as $\mu_i > \mu_j$ when in fact $\mu_i < \mu_j$, and power counts those outcomes as successes. Sizing for narrow intervals alone has the same problem, because it counts narrow intervals that miss the truth. Hsu's criterion counts a test as a success only if its intervals are correct and narrow.

Power also addresses a different question: whether the test will detect that the arms differ. The planning question here is whether the test will identify an arm close to the best and rule out the clearly worse ones. In Hsu's alloy example, the $F$-test power calculation (for detecting differences of 10,000 psi with probability 0.90) gives 5 specimens per alloy, and the correct and useful MCB calculation gives 6.

## How the probability is computed

Let $S = \hat{\sigma}/\sigma$, with density $\gamma_\nu$ for $\nu = k(n-1)$, and let

$$
u = \sqrt{n/2}\,\frac{\delta^*/\sigma}{d}.
$$

The condition $w < \delta^*$ is equivalent to $S < u$. The critical values are the balanced ones, with every $\lambda = 1/\sqrt{2}$.

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

For constrained MCB, the two conclusions are Hsu's statements in Appendix C.2. For Edwards-Hsu, Hsu says only that the computation is "essentially the same" as for two-sided MCC (Appendix C.3). Both conclusions follow from the correctness event in the proof of Theorem 4.2.3, by the same argument as above.

The correctness event implies that every interval covers its target, but the intervals can also cover their targets when it does not occur. Both values are therefore conservative: they are lower bounds on the probability that the intervals cover their targets and are narrow.

The `*_sample_size` functions return the smallest $n \ge 2$ at which the corresponding probability reaches `target`.

## Scope and numerical details

- The computation covers balanced designs only, as Appendix C does: one common $n$ for every arm. The analysis functions (`constrained_mcb`, `unconstrained_mcb_edwards_hsu`) accept unequal sample sizes.
- The upper limit of $1 - \alpha$ holds up to the critical-value solver's tolerance of about $2 \times 10^{-5}$.
- The critical values switch to their infinite-degrees-of-freedom form at $\nu = 10^6$ (see {doc}`mcb_critical_values`). Above the switch, the narrowness condition still uses the exact chi-square distribution of $S$, while the correctness part is evaluated at $\hat{\sigma} = \sigma$, which changes the probability by about $10^{-4}$. The critical value just below the switch differs from the infinite-df one by slightly more than the true finite-df correction, within the solver's tolerance. As a result, for $k \ge 3$ the probability jumps by about $10^{-3}$ across the switch. It still increases with $n$, and a sample size near the switch is off by a few units out of $10^5$ to $10^6$.
