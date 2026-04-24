# Model and Notation

The package uses Hsu's unbalanced one-way normal model. Arm $i$ has observations

$$
Y_{ia} = \mu_i + \varepsilon_{ia},
\qquad
i = 1,\ldots,k,
\qquad
a = 1,\ldots,n_i,
$$

with independent errors

$$
\varepsilon_{ia} \stackrel{iid}{\sim} N(0,\sigma^2).
$$

The larger-is-better convention is used throughout these docs and in the package. If smaller values are better in your application, multiply the response by `-1` before forming the summary.

## Sample Summaries

The observed arm mean is

$$
\hat{\mu}_i = \bar{Y}_i = \frac{1}{n_i}\sum_{a=1}^{n_i}Y_{ia}.
$$

The pooled within-arm variance estimate is

$$
MSE
= \hat{\sigma}^2
=
\frac{
\sum_{i=1}^k\sum_{a=1}^{n_i}(Y_{ia} - \bar{Y}_i)^2
}{
\nu
},
\qquad
\nu = \sum_{i=1}^k(n_i - 1).
$$

`mcb` can create this summary three ways:

- `summarize_from_samples`: start from raw observations.
- `summarize_from_means_ns_and_mse`: start from per-arm means, sample sizes, and pooled MSE.
- `pooled_mse_from_group_stats`: compute pooled MSE from per-arm sample variances.

## Target Parameters

There are two related but different MCB targets.

The constrained procedure targets

$$
\theta_i = \mu_i - \max_{j \neq i}\mu_j.
$$

This compares arm $i$ with the best of the other arms. If $\theta_i > 0$, arm $i$ is better than every competitor. If $\theta_i < 0$, arm $i$ is not best. If $\theta_i$ is close to 0, arm $i$ is close to the best competitor.

The Edwards-Hsu unconstrained procedure targets

$$
\eta_i = \mu_i - \max_{1 \leq j \leq k}\mu_j.
$$

This compares arm $i$ with the global best arm, including itself. Therefore $\eta_i \leq 0$ always, and $\eta_i = 0$ exactly for arms tied for the true best.

## Pairwise Standard Errors

The standard error for comparing arms $i$ and $j$ is

$$
\hat{\sigma}\sqrt{n_i^{-1} + n_j^{-1}}.
$$

For unbalanced critical values, Hsu rewrites the standardized comparison using

$$
\lambda_{ij}
= \left(1 + \frac{n_i}{n_j}\right)^{-1/2}
= \sqrt{\frac{n_j}{n_i+n_j}},
\qquad j \neq i.
$$

Because the vector $(\lambda_{ij}: j \neq i)$ depends on the focal arm $i$, unbalanced designs generally have different critical values for different arms.

## Candidate-Best Sets

Both implemented procedures return `candidate_best_set`, but the exact construction differs:

- `constrained_mcb` returns $G$, the arms whose constrained upper MCB bound is positive.
- `unconstrained_mcb_edwards_hsu` returns $S$, the arms that pass the Edwards-Hsu screening inequality.

These sets are confidence-set outputs. Membership means the data do not rule the arm out as best at the requested familywise level. It is not a posterior probability and not a guarantee that the arm is best.
