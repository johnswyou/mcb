# MCB procedures in `mcb`

This page summarizes the two procedures implemented by the package and how they map to Hsu's unbalanced one-way model.

## Setup

Let variant $i$ denote treatment arm $i$, with observations

$$
Y_{ia} = \mu_i + \varepsilon_{ia},
\qquad
\varepsilon_{ia} \stackrel{iid}{\sim} N(0, \sigma^2),
$$

and pooled within-arm variance estimate $MSE = \hat{\sigma}^2$. In the package, you can supply this pooled variance directly through `summarize_from_means_ns_and_mse`, compute it from raw samples with `pooled_mse_from_samples`, or compute it from per-arm sample variances with `pooled_mse_from_group_stats`.

The two parameterizations to keep distinct are

$$
\theta_i = \mu_i - \max_{j \neq i} \mu_j
$$

and

$$
\eta_i = \mu_i - \max_{1 \leq j \leq k} \mu_j.
$$

The constrained procedure targets $\theta_i$. The direct Edwards-Hsu unconstrained procedure targets $\eta_i$.

## Constrained MCB

Use `constrained_mcb` when your main question is, "which variants can still plausibly be the best, and how close might they be to the best alternative?"

For each arm $i$, Hsu defines a one-sided Dunnett-with-control critical value $d^i$ that depends on the familywise error rate $\alpha$, the pooled error degrees of freedom

$$
\nu = \sum_i (n_i - 1),
$$

and the sample-size pattern through

$$
\lambda_{ij} = \left(1 + \frac{n_i}{n_j}\right)^{-1/2},
\qquad j \neq i.
$$

Once $d^i$ is known, the constrained upper bound is

$$
D_i^+ =
\left[
\min_{j \neq i}
\left\{
\hat{\mu}_i - \hat{\mu}_j + d^i \hat{\sigma} \sqrt{n_i^{-1} + n_j^{-1}}
\right\}
\right]^+,
$$

where $x^+ = \max(x, 0)$. The candidate-best set is

$$
G = \{ i : D_i^+ > 0 \}.
$$

If $D_i^+ = 0$, variant $i$ is ruled out as the true winner at familywise level $\alpha$.

The lower constrained bound is

$$
D_i^- =
\begin{cases}
0, & G = \{ i \}, \\[4pt]
\min_{j \in G,\ j \neq i}
\left\{
\hat{\mu}_i - \hat{\mu}_j - d^j \hat{\sigma} \sqrt{n_i^{-1} + n_j^{-1}}
\right\}, & \text{otherwise.}
\end{cases}
$$

The simultaneous intervals $[D_i^-, D_i^+]$ cover all $\theta_i = \mu_i - \max_{j \neq i} \mu_j$ jointly with probability at least $1 - \alpha$.

### Interpretation

- If the interval is entirely negative, the variant is confidently worse than the true best.
- If the interval touches 0, the data do not rule out that the variant could still be best.
- If sample sizes differ, different arms can have different critical values.

## Unconstrained Edwards-Hsu MCB

Use `unconstrained_mcb_edwards_hsu` when you want direct statements about how far each arm is below the global best.

This procedure targets

$$
\eta_i = \mu_i - \max_{1 \leq j \leq k} \mu_j \leq 0.
$$

For each arm $i$, let $|d|^i$ be the two-sided Dunnett-with-control critical value computed from the same $\lambda_{ij}$ values. Define the screening set

$$
S =
\left\{
i :
\min_{j \neq i}
\left\{
\hat{\mu}_i - \hat{\mu}_j + |d|^i \hat{\sigma} \sqrt{n_i^{-1} + n_j^{-1}}
\right\}
> 0
\right\}.
$$

Only arms in $S$ remain plausible winners.

For $j \in S$, define

$$
L_{ij} =
\begin{cases}
0, & i = j, \\[4pt]
\hat{\mu}_i - \hat{\mu}_j - |d|^j \hat{\sigma} \sqrt{n_i^{-1} + n_j^{-1}}, & i \neq j,
\end{cases}
\qquad
L_i = \min_{j \in S} L_{ij},
$$

and

$$
U_{ij} =
\begin{cases}
0, & i = j, \\[4pt]
-\left(\hat{\mu}_i - \hat{\mu}_j + |d|^j \hat{\sigma} \sqrt{n_i^{-1} + n_j^{-1}}\right)^-, & i \neq j,
\end{cases}
\qquad
U_i = \max_{j \in S} U_{ij},
$$

where $x^- = \max(-x, 0)$ is the magnitude of the negative part.

The intervals $[L_i, U_i]$ cover all $\eta_i = \mu_i - \max_j \mu_j$ jointly with probability at least $1 - \alpha$.

### Interpretation

- If $U_i < 0$, arm $i$ is decisively inferior to the true best.
- If $U_i = 0$, the arm remains compatible with being the global winner.
- These intervals are naturally phrased in terms of the gap to the overall best arm, not the best competing arm.

## Critical values in the package

The package computes both critical values numerically from Hsu's Appendix D.2 construction:

- 48-point Gauss-Hermite quadrature for the inner normal integral.
- 48-point Gauss-Legendre quadrature for the outer integral over $S = \hat{\sigma} / \sigma$.
- Root bracketing from $t$-distribution quantiles with an Illinois/regula falsi solver.
- A switch to the infinite-degrees-of-freedom approximation when $\nu > 240`.

See [MCB critical values](mcb_critical_values.md) for the formulas used by `constrained_critical_value` and `unconstrained_edwards_hsu_critical_value`.
