# Constrained MCB

`constrained_mcb` computes simultaneous confidence intervals $[D_i^-, D_i^+]$ for

$$
\theta_i = \mu_i - \max_{j \neq i}\mu_j,
\qquad i = 1,\ldots,k,
$$

arm $i$'s mean minus the best of the other arms' means. The intervals are called *constrained* because every one of them contains 0. The construction is Hsu (1996), Section 4.1.4, equations (4.10) to (4.16). The notation is defined in [Model and notation](model_notation.md).

## Critical values

For each arm $i$, $d^i$ is the critical value of Dunnett's one-sided method for comparing every other arm with a control, with arm $i$ as the control. It solves

$$
P\left\{
\hat{\mu}_i - \mu_i
>
\hat{\mu}_j - \mu_j
- d^i\hat{\sigma}\sqrt{n_i^{-1} + n_j^{-1}}
\text{ for all } j \neq i
\right\}
= 1-\alpha.
$$

In words: with probability $1-\alpha$, no other arm's estimation error exceeds arm $i$'s by $d^i$ standard errors or more. [MCB critical values](mcb_critical_values.md) describes the computation.

## Upper bounds and the candidate-best set

The upper bound is

$$
D_i^+
=
\left[
\min_{j \neq i}
\left\{
\hat{\mu}_i - \hat{\mu}_j
+ d^i\hat{\sigma}\sqrt{n_i^{-1} + n_j^{-1}}
\right\}
\right]^+,
\qquad x^+ = \max(x,0).
$$

$D_i^+$ is 0 when some other arm's sample mean exceeds arm $i$'s by at least $d^i$ standard errors. Arm $i$ is then ruled out as the best. The arms that are not ruled out form the candidate-best set

$$
G = \{i : D_i^+ > 0\}.
$$

## Lower bounds

The lower bound is

$$
D_i^- =
\begin{cases}
0, & G = \{i\}, \\[6pt]
\displaystyle
\min_{j \in G,\ j \neq i}
\left\{
\hat{\mu}_i - \hat{\mu}_j
- d^j\hat{\sigma}\sqrt{n_i^{-1} + n_j^{-1}}
\right\},
& \text{otherwise.}
\end{cases}
$$

Arm $i$ is compared only with the arms in $G$, because only they could be the best. The critical value is $d^j$, not $d^i$: in each comparison the candidate $j$ is the control.

$D_i^-$ is never positive. Arm $j$ is in $G$ only if $\hat{\mu}_j - \hat{\mu}_i + d^j\hat{\sigma}\sqrt{n_i^{-1} + n_j^{-1}} > 0$, which makes every term in the minimum negative. Since $D_i^+$ is never negative, every interval contains 0.

The intervals $[D_i^-, D_i^+]$ contain all the $\theta_i$ simultaneously with probability at least $1-\alpha$ (Hsu's Theorem 4.1.2).

## Reading the intervals

Each interval has one of three forms.

- $[D_i^-,\ 0]$ with $D_i^- < 0$. Arm $i$ is not the best. It is at most $-D_i^-$ below the best.
- $[0,\ D_i^+]$. Arm $i$ is the best, and it leads the runner-up by at most $D_i^+$. This form occurs only when $G = \{i\}$.
- $[D_i^-,\ D_i^+]$ with $D_i^- < 0 < D_i^+$. Arm $i$ could be the best. It is at most $-D_i^-$ below the best of the other arms and at most $D_i^+$ above it.

A lower bound near 0 shows that an arm is close to the best even when it cannot be shown to be the best. [Assumptions and scope](assumptions.md#ties) covers how ties change the first two readings.

Because every interval contains 0, constrained MCB never states that an arm is at least some amount below the best, or that the best arm leads by at least some amount. [Edwards-Hsu MCB](edwards_hsu.md) can make such statements, with larger critical values.

## Balanced designs

When every arm has $n$ observations, all arms share one critical value $d$ and one margin $w = d\,\hat{\sigma}\sqrt{2/n}$, and the bounds reduce to

$$
D_i^- = \min\left(0,\ \hat{\mu}_i - \max_{j \neq i}\hat{\mu}_j - w\right),
\qquad
D_i^+ = \max\left(0,\ \hat{\mu}_i - \max_{j \neq i}\hat{\mu}_j + w\right)
$$

(Hsu's equation (4.5)). The [sample-size computation](sample_size.md) uses this form.
