# Model and notation

## The one-way normal model

Arm $i$ has $n_i$ observations

$$
Y_{ia} = \mu_i + \varepsilon_{ia},
\qquad
i = 1,\ldots,k,
\qquad
a = 1,\ldots,n_i,
$$

where the errors $\varepsilon_{ia}$ are independent $N(0,\sigma^2)$ and $\sigma^2$ is unknown. The design is *balanced* if every $n_i$ is the same and *unbalanced* otherwise. The analysis functions handle both.

Larger responses are better throughout the package and these pages. If smaller is better, multiply the response by $-1$ before the analysis.

## Sample summaries

The procedures use the data only through the arm means, the sample sizes and the pooled variance estimate:

$$
\hat{\mu}_i = \bar{Y}_i = \frac{1}{n_i}\sum_{a=1}^{n_i}Y_{ia},
\qquad
\hat{\sigma}^2 = MSE
=
\frac{1}{\nu}\sum_{i=1}^k\sum_{a=1}^{n_i}(Y_{ia} - \bar{Y}_i)^2,
\qquad
\nu = \sum_{i=1}^k(n_i - 1).
$$

$\nu$ is the error degrees of freedom. A `OneWaySummary` holds these quantities. There are two ways to build one:

- `summarize_from_samples`, from the raw observations;
- `summarize_from_means_ns_and_mse`, from the arm means, the sample sizes and the MSE. If only per-arm sample variances are known, `pooled_mse_from_group_stats` computes the MSE from them.

## Target parameters

$$
\theta_i = \mu_i - \max_{j \neq i}\mu_j,
\qquad
\eta_i = \mu_i - \max_{1 \leq j \leq k}\mu_j = \min(\theta_i, 0).
$$

Constrained MCB gives intervals for $\theta_i$, arm $i$'s mean minus the best of the other arms' means. Edwards-Hsu MCB gives intervals for $\eta_i$, arm $i$'s mean minus the best mean. [Overview](mcb.md) explains how to read them. The `estimates` field of a result is the target with each $\mu$ replaced by $\hat{\mu}$.

## Standard errors and $\lambda_{ij}$

The estimated standard error of $\hat{\mu}_i - \hat{\mu}_j$ is

$$
\hat{\sigma}\sqrt{n_i^{-1} + n_j^{-1}}.
$$

Every bound in both procedures is a difference of two sample means plus or minus a critical value times this standard error. Each arm $i$ has its own critical value, written $d^i$ for constrained MCB and $|d|^i$ for Edwards-Hsu MCB. Apart from $\alpha$, $k$ and $\nu$, it depends on the sample sizes through

$$
\lambda_{ij}
= \left(1 + \frac{n_i}{n_j}\right)^{-1/2}
= \sqrt{\frac{n_j}{n_i+n_j}},
\qquad j \neq i.
$$

In a balanced design every $\lambda_{ij}$ equals $1/\sqrt{2}$ and all arms share one critical value. In an unbalanced design the critical values differ from arm to arm. [MCB critical values](mcb_critical_values.md) gives the computation.

## Candidate-best set

Both procedures return a `candidate_best_set`, the arms that the data do not rule out as the best. It is a confidence set for the best arm: it contains the true best arm with probability at least $1-\alpha$. Constrained MCB calls the set $G$ and Edwards-Hsu MCB calls it $S$. Each is defined on its procedure's page.

Membership means the data cannot exclude the arm. It says nothing about how probable it is that the arm is the best.
