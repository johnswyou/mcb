# Edwards-Hsu unconstrained MCB

`unconstrained_mcb_edwards_hsu` computes simultaneous confidence intervals $[L_i, U_i]$ for

$$
\eta_i = \mu_i - \max_{1 \leq j \leq k}\mu_j,
\qquad i = 1,\ldots,k,
$$

arm $i$'s mean minus the best mean. $\eta_i$ is 0 for the best arm and negative for every other arm. These intervals are called *unconstrained* because they need not contain 0, so they can show that an arm is at least some amount below the best. The construction is Hsu (1996), Section 4.2.4.2, equations (4.44) and (4.45). The notation is defined in [Model and notation](model_notation.md).

## Critical values

For each arm $i$, $|d|^i$ is the critical value of Dunnett's two-sided method for comparing every other arm with a control, with arm $i$ as the control. It solves

$$
P\left\{
\left|
\hat{\mu}_i - \mu_i - (\hat{\mu}_j - \mu_j)
\right|
<
|d|^i\hat{\sigma}\sqrt{n_i^{-1} + n_j^{-1}}
\text{ for all } j \neq i
\right\}
= 1-\alpha.
$$

$|d|^i$ is Hsu's symbol for this constant. It is not the absolute value of the constrained critical value $d^i$, and it is larger than $d^i$. [MCB critical values](mcb_critical_values.md) describes the computation.

## Screening set

The candidate-best set of this procedure is the screening set

$$
S =
\left\{
i:
\min_{j \neq i}
\left[
\hat{\mu}_i - \hat{\mu}_j
+ |d|^i\hat{\sigma}\sqrt{n_i^{-1} + n_j^{-1}}
\right]
> 0
\right\}.
$$

Arm $i$ is left out of $S$ when some other arm's sample mean exceeds its own by at least $|d|^i$ standard errors.

## Interval bounds

Each arm $j \in S$ could be the best. If arm $j$ is the best, then $\eta_i = \mu_i - \mu_j$, which is bounded by

$$
L_{ij} =
\begin{cases}
0, & i=j, \\[6pt]
\hat{\mu}_i - \hat{\mu}_j
- |d|^j\hat{\sigma}\sqrt{n_i^{-1} + n_j^{-1}},
& i \neq j,
\end{cases}
\qquad
U_{ij} =
\begin{cases}
0, & i=j, \\[6pt]
\min\left(0,\
\hat{\mu}_i - \hat{\mu}_j
+ |d|^j\hat{\sigma}\sqrt{n_i^{-1} + n_j^{-1}}
\right),
& i \neq j.
\end{cases}
$$

It is not known which arm in $S$ is the best, so the final bounds take the widest range over $S$:

$$
L_i = \min_{j \in S} L_{ij},
\qquad
U_i = \max_{j \in S} U_{ij}.
$$

The intervals $[L_i, U_i]$ contain all the $\eta_i$ simultaneously with probability at least $1-\alpha$ (Hsu's Theorem 4.2.6).

## Reading the intervals

- If $U_i < 0$, arm $i$ is not the best. It is between $-U_i$ and $-L_i$ below the best.
- If $U_i = 0$, the data do not show arm $i$ to be any amount below the best. It is at most $-L_i$ below.
- If $S$ contains a single arm $M$, then $L_M = U_M = 0$ and arm $M$ is the best.

Use $S$ (the `candidate_best_set` field), not $U_i = 0$, to see which arms could be the best. Every arm in $S$ has $U_i = 0$, but an arm outside $S$ can have $U_i = 0$ as well: it is ruled out as the best without being shown to be any particular amount below it.

The reported estimate is $\hat{\mu}_i - \max_j \hat{\mu}_j$, which is 0 for the arm with the highest sample mean.

An interval for $\eta_i$ that lies entirely below 0 is also an interval for $\theta_i = \mu_i - \max_{j \neq i}\mu_j$, because the two parameters are equal for an arm that is not the best.
