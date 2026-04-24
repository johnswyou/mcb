# Constrained MCB

`constrained_mcb` implements Hsu's constrained Multiple Comparisons with the Best procedure for the unbalanced one-way model. The source construction is Hsu (1996), Section 4.1.4, equations (4.10) through (4.16).

Use it when you want to identify arms that are not compatible with being best, while also showing how close the remaining arms may be to their best competitor.

## Target

The target parameter is

$$
\theta_i = \mu_i - \max_{j \neq i}\mu_j.
$$

The simultaneous intervals are constrained to contain 0. This reflects the role of 0 as the decision boundary:

- $\theta_i > 0$ means arm $i$ is better than every other arm.
- $\theta_i < 0$ means some other arm is better than arm $i$.
- $\theta_i$ near 0 means arm $i$ is close to the best competing arm.

## Arm-Specific Critical Values

For each focal arm $i$, let $d_i$ be the one-sided Dunnett-with-control critical value computed as if arm $i$ were the control. It is chosen so that

$$
P\left\{
\hat{\mu}_i - \mu_i
>
\hat{\mu}_j - \mu_j
- d_i\hat{\sigma}\sqrt{n_i^{-1} + n_j^{-1}}
\text{ for all } j \neq i
\right\}
= 1-\alpha.
$$

In an unbalanced design, $d_i$ depends on the full sample-size pattern through the values $\lambda_{ij}$, so the critical value can differ by arm.

## Upper Bounds and Candidate-Best Set

The constrained upper bound is

$$
D_i^+
=
\left[
\min_{j \neq i}
\left\{
\hat{\mu}_i - \hat{\mu}_j
+ d_i\hat{\sigma}\sqrt{n_i^{-1} + n_j^{-1}}
\right\}
\right]^+,
$$

where $x^+ = \max(x,0)$.

The candidate-best set is

$$
G = \{i : D_i^+ > 0\}.
$$

If $D_i^+ = 0$, arm $i$ is inferred not to be the true best at familywise level $\alpha$. If $D_i^+ > 0$, arm $i$ remains compatible with being best.

## Lower Bounds

Once $G$ is known, the lower constrained bound is

$$
D_i^- =
\begin{cases}
0, & G = \{i\}, \\[6pt]
\displaystyle
\min_{j \in G,\ j \neq i}
\left\{
\hat{\mu}_i - \hat{\mu}_j
- d_j\hat{\sigma}\sqrt{n_i^{-1} + n_j^{-1}}
\right\},
& \text{otherwise.}
\end{cases}
$$

The critical value in the lower bound is $d_j$, not $d_i$, because the possible best arm $j$ is treated as the control in that comparison.

The result is a set of simultaneous intervals

$$
[D_i^-,D_i^+],
\qquad i=1,\ldots,k,
$$

with joint coverage at least $1-\alpha$ for all $\theta_i$ under the one-way normal model.

## Interpretation

Read the constrained interval against 0:

- If the upper bound is 0, the arm is ruled out as best.
- If the interval includes values near 0, the arm could be best or close to best.
- If $G$ contains exactly one arm, that arm alone remains in the constrained candidate-best set.

Constrained MCB does not provide a positive lower bound on how much the selected best arm exceeds the others. Hsu's point is that giving up that lower-bound claim produces sharper inference for identifying arms that are not best.
