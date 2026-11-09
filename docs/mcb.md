# Overview

An experiment compares $k$ arms (treatments, variants) on a response for which larger is better. Often the questions that matter are which arms could be the best and how far the others are from it. In that case, how two arms that are both far from the best compare with each other is of no interest. Multiple Comparisons with the Best (MCB) answers the two questions with one confidence interval per arm.

The intervals are *simultaneous*: with probability at least $1 - \alpha$, all $k$ of them contain their true values at once. $\alpha$ is the *familywise error rate*.

MCB compares each arm only with the best, so its intervals are narrower than those of a method that covers all $k(k-1)/2$ pairwise differences. Each MCB bound is a difference of two sample means plus or minus a critical value times a standard error. In the three-arm [SAT example](sat_example.md) at $\alpha = 0.10$, the critical values are 1.50 to 1.61 for constrained MCB and 1.86 to 1.95 for Edwards-Hsu MCB, against 2.05 for the Tukey-Kramer all-pairwise method (Hsu 1996, Section 4.2.5).

## Two target parameters

`mcb` implements two procedures. They differ in what each arm's interval is an interval for.

`constrained_mcb` targets arm $i$'s mean minus the best of the *other* arms' means:

$$
\theta_i = \mu_i - \max_{j \neq i} \mu_j.
$$

If $\theta_i > 0$, arm $i$ is the best and leads the runner-up by $\theta_i$. If $\theta_i < 0$, arm $i$ is not the best and is $-\theta_i$ below it.

`unconstrained_mcb_edwards_hsu` targets arm $i$'s mean minus the best mean, where the best may be arm $i$ itself:

$$
\eta_i = \mu_i - \max_{1 \leq j \leq k} \mu_j = \min(\theta_i, 0).
$$

$\eta_i$ is 0 for the best arm and negative for every other arm. For an arm that is not the best, the two targets are equal.

## Which procedure to use

| | `constrained_mcb` | `unconstrained_mcb_edwards_hsu` |
| --- | --- | --- |
| Target | $\theta_i$ | $\eta_i$ |
| Interval contains 0 | always | only if the arm is not shown to be below the best |
| Critical value | one-sided, smaller | two-sided, larger |
| For an arm ruled out as the best | at most how far below the best it is | at most how far below and, when the data allow, at least how far below |
| Per-arm sample size for the same guarantee | smaller | larger |

The SAT example shows the trade-off. For the arm `Other`, constrained MCB gives the interval $[-78.97,\ 0]$: the arm is not the best and is at most 78.97 points below it. Edwards-Hsu MCB gives $[-84.87,\ -22.56]$: the arm is between 22.56 and 84.87 points below the best. In exchange, Edwards-Hsu MCB places the other arms less precisely: Computer science is at most 41.44 points below the best, against 35.43 under constrained MCB.

Use `constrained_mcb` unless you need to state that an arm is at least some amount below the best.

Both functions return an `MCBResult` with the same fields, including `candidate_best_set`, the arms that the data do not rule out as the best.
