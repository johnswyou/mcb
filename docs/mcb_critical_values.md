# MCB critical values

`constrained_critical_value` computes $d^i$, the critical value of [constrained MCB](constrained_mcb.md). `unconstrained_edwards_hsu_critical_value` computes $|d|^i$, the critical value of [Edwards-Hsu MCB](edwards_hsu.md). The analysis functions call them once per arm. They can also be called directly:

```python
from math import sqrt

from mcb import constrained_critical_value, unconstrained_edwards_hsu_critical_value

# Balanced design, k = 3 arms, n = 11 per arm: nu = 3 * 10 = 30, and both lambdas are 1 / sqrt(2).
lambdas = (1 / sqrt(2),) * 2
constrained_critical_value(lambdas, nu=30, alpha=0.05)                # 1.989
unconstrained_edwards_hsu_critical_value(lambdas, nu=30, alpha=0.05)  # 2.321
```

These are the values in Hsu's Tables E.2 and E.3. The argument `lambdas` holds $\lambda_{ij}$ for every $j \neq i$, as defined below.

## Setup

Fix an arm $i$ as the control. Under the [one-way normal model](model_notation.md), the standardized estimation errors

$$
Z_m = \frac{\sqrt{n_m}\,(\hat{\mu}_m - \mu_m)}{\sigma},
\qquad m = 1,\ldots,k,
$$

are independent standard normal variables. They are independent of $\hat{\sigma}/\sigma$, whose density is written $\gamma_\nu$; $\nu\hat{\sigma}^2/\sigma^2$ has a chi-square distribution with $\nu$ degrees of freedom.

The sample sizes enter through

$$
\lambda_{ij} =
\left(1 + \frac{n_i}{n_j}\right)^{-1/2}
=
\sqrt{\frac{n_j}{n_i + n_j}},
\qquad j \neq i.
$$

The standardized difference between the errors of arm $i$ and arm $j$ is

$$
T_{ij}
=
\frac{(\hat{\mu}_i - \mu_i) - (\hat{\mu}_j - \mu_j)}
{\sigma \sqrt{n_i^{-1} + n_j^{-1}}}
=
\lambda_{ij} Z_i - \sqrt{1 - \lambda_{ij}^2}\, Z_j.
$$

Given $Z_i = z$ and $\hat{\sigma}/\sigma = s$, the $T_{ij}$ for different $j$ are independent. The probability of a joint event over all $j \neq i$ is therefore a product of normal probabilities, averaged over $z$ and $s$.

## Constrained critical value $d^i$

$d^i$ is the value of $d$ at which the event $T_{ij} > -d\,\hat{\sigma}/\sigma$ for all $j \neq i$ has probability $1 - \alpha$. That probability is

$$
p_i^{\mathrm{con}}(d)
=
\int_0^\infty
\int_{-\infty}^{+\infty}
\prod_{j \neq i}
\Phi\left(
\frac{\lambda_{ij} z + d s}{\sqrt{1 - \lambda_{ij}^2}}
\right)
\, d\Phi(z)\, \gamma_\nu(s)\, ds
$$

(Hsu's equation (4.11)), where $\Phi$ is the standard normal distribution function. The solution of $p_i^{\mathrm{con}}(d^i) = 1 - \alpha$ lies in

$$
t_{\alpha,\nu} \leq d^i \leq t_{\alpha/(k-1),\nu},
$$

where $t_{\beta,\nu}$ is the upper $\beta$ quantile of the $t$ distribution with $\nu$ degrees of freedom. For $k = 2$ the two ends coincide and $d^i = t_{\alpha,\nu}$.

## Edwards-Hsu critical value $|d|^i$

$|d|^i$ is the value of $c$ at which the two-sided event $|T_{ij}| < c\,\hat{\sigma}/\sigma$ for all $j \neq i$ has probability $1 - \alpha$. That probability is

$$
p_i^{\mathrm{EH}}(c)
=
\int_0^\infty
\int_{-\infty}^{+\infty}
\prod_{j \neq i}
\left[
\Phi\left(
\frac{\lambda_{ij} z + c s}{\sqrt{1 - \lambda_{ij}^2}}
\right)
-
\Phi\left(
\frac{\lambda_{ij} z - c s}{\sqrt{1 - \lambda_{ij}^2}}
\right)
\right]
\, d\Phi(z)\, \gamma_\nu(s)\, ds.
$$

The solution of $p_i^{\mathrm{EH}}(|d|^i) = 1 - \alpha$ lies in

$$
t_{\alpha/2,\nu} \leq |d|^i \leq t_{\alpha/[2(k-1)],\nu}.
$$

For $k = 2$, $|d|^i = t_{\alpha/2,\nu}$.

## Numerical method

The computation follows Hsu's Appendix D.2:

- The inner integral over $z$ uses 48-point Gauss-Hermite quadrature.
- The outer integral over $s$ uses 48-point Gauss-Legendre quadrature on the central $1 - 10^{-6}$ probability interval of $\hat{\sigma}/\sigma$.
- The root is found with the Illinois variant of regula falsi, starting from the interval of $t$ quantiles given above. A candidate is accepted when its probability is within $2\times 10^{-5}$ of $1-\alpha$ or when the half-width of the bracketing interval is at most $5\times 10^{-4}$.

When $\nu = \infty$, $\hat{\sigma}$ equals $\sigma$ and only the inner integral remains. The package uses this form when `nu` is infinite or at least $10^6$, where it agrees with the finite-$\nu$ form to within the solver's tolerance. Hsu's programs made the same switch at $\nu > 240$. The package does not, so some of its values differ slightly from the book's (see the [SAT example](sat_example.md#differences-from-the-printed-textbook-values)).

## Balanced designs and the Appendix E tables

In a balanced design $\lambda_{ij} = 1/\sqrt{2}$ for every pair, so the critical values depend only on $\alpha$, $k$ and $\nu$ and can be tabulated. Hsu's Appendix E does this, and the package's test suite checks its results against Tables E.2, E.3, E.5 and E.6.

In an unbalanced design the critical values also depend on the sample-size ratios and differ from arm to arm, so they are computed for each design.
