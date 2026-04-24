# Edwards-Hsu Unconstrained MCB

`unconstrained_mcb_edwards_hsu` implements the direct Edwards-Hsu unconstrained MCB method for the unbalanced one-way model. The source construction is Hsu (1996), Section 4.2.4.2, equations (4.44) and (4.45).

Use it when you want direct simultaneous intervals for how far each arm is below the global best.

## Target

The target parameter is

$$
\eta_i = \mu_i - \max_{1 \leq j \leq k}\mu_j.
$$

This parameter is always nonpositive:

$$
\eta_i \leq 0.
$$

An arm has $\eta_i = 0$ if it is tied for the true best. If $\eta_i < 0$, the arm is below the best by $|\eta_i|$ units.

## Arm-Specific Critical Values

For each focal arm $i$, let $|d|_i$ be the two-sided Dunnett-with-control critical value computed as if arm $i$ were the control:

$$
P\left\{
\left|
\hat{\mu}_i - \mu_i - (\hat{\mu}_j - \mu_j)
\right|
<
|d|_i\hat{\sigma}\sqrt{n_i^{-1} + n_j^{-1}}
\text{ for all } j \neq i
\right\}
= 1-\alpha.
$$

Like the constrained critical value, $|d|_i$ is arm-specific in an unbalanced design.

## Screening Set

The Edwards-Hsu screening set is

$$
S =
\left\{
i:
\min_{j \neq i}
\left[
\hat{\mu}_i - \hat{\mu}_j
+ |d|_i\hat{\sigma}\sqrt{n_i^{-1} + n_j^{-1}}
\right]
> 0
\right\}.
$$

Only arms in $S$ remain compatible with being the global best.

## Interval Bounds

For each possible best arm $j \in S$, define

$$
L_{ij} =
\begin{cases}
0, & i=j, \\[6pt]
\hat{\mu}_i - \hat{\mu}_j
- |d|_j\hat{\sigma}\sqrt{n_i^{-1} + n_j^{-1}},
& i \neq j,
\end{cases}
$$

and

$$
U_{ij} =
\begin{cases}
0, & i=j, \\[6pt]
-\left(
\hat{\mu}_i - \hat{\mu}_j
+ |d|_j\hat{\sigma}\sqrt{n_i^{-1} + n_j^{-1}}
\right)^-,
& i \neq j,
\end{cases}
$$

where $x^- = \max(-x,0)$.

The final bounds are

$$
L_i = \min_{j \in S} L_{ij},
\qquad
U_i = \max_{j \in S} U_{ij}.
$$

The intervals

$$
[L_i,U_i],
\qquad i=1,\ldots,k,
$$

jointly cover all $\eta_i = \mu_i - \max_j\mu_j$ with probability at least $1-\alpha$ under the one-way normal model.

## Interpretation

Read the Edwards-Hsu interval as a gap-to-best interval:

- If $U_i < 0$, arm $i$ is inferred to be below the true best.
- If $U_i = 0$, arm $i$ remains compatible with being globally best.
- The estimate reported by `mcb` is $\hat{\mu}_i - \max_j \hat{\mu}_j$, so the sample-best arm has estimate 0.

The Edwards-Hsu result targets $\eta_i$, not $\theta_i$. Hsu notes that when an interval for $\eta_i$ is entirely negative, the same numerical interval can be read as an interval for $\mu_i - \max_{j \neq i}\mu_j$, because the best arm cannot be $i$ in that case. The package reports the direct $\eta_i$ intervals.
