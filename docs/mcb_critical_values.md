# MCB critical values

This page describes the two critical values exposed by the package:

- `constrained_critical_value` for Hsu's constrained MCB critical value $d^i$
- `unconstrained_edwards_hsu_critical_value` for the direct Edwards-Hsu critical value $|d|^i$

Both are computed arm-by-arm in the unbalanced one-way model, so the value for a focal arm depends on that arm's full pattern of pairwise sample-size ratios.

## Common setup

Under Hsu's model,

$$
Y_{ia} = \mu_i + \varepsilon_{ia},
\qquad
\varepsilon_{ia} \stackrel{iid}{\sim} N(0, \sigma^2),
$$

with pooled variance estimator $\hat{\sigma}^2 = MSE$ and

$$
\nu = \sum_{i=1}^k (n_i - 1).
$$

For a fixed focal arm $i$, define

$$
Z_i = \frac{\sqrt{n_i}(\hat{\mu}_i - \mu_i)}{\sigma},
\qquad
Z_j = \frac{\sqrt{n_j}(\hat{\mu}_j - \mu_j)}{\sigma},
\qquad
S = \frac{\hat{\sigma}}{\sigma}.
$$

Then $Z_i$ and the $Z_j$ values are independent standard normals, and $\nu S^2 \sim \chi^2_\nu$.

The key sample-size-dependent quantity is

$$
\lambda_{ij} =
\left(1 + \frac{n_i}{n_j}\right)^{-1/2}
=
\sqrt{\frac{n_j}{n_i + n_j}},
\qquad j \neq i.
$$

The standardized comparison between the focal arm and arm $j$ is

$$
U_{ij}
=
\frac{(\hat{\mu}_i - \mu_i) - (\hat{\mu}_j - \mu_j)}
{\sigma \sqrt{n_i^{-1} + n_j^{-1}}}
=
\lambda_{ij} Z_i - \sqrt{1 - \lambda_{ij}^2}\, Z_j.
$$

Conditioning on $Z_i = z$ turns the different comparisons into conditionally independent normal events, which is why Hsu's integrals factor over $j \neq i$.

## Constrained critical value $d^i$

For constrained MCB, $d^i$ is the one-sided Dunnett-with-control critical value satisfying

$$
P\left\{
\hat{\mu}_i - \mu_i
>
\hat{\mu}_j - \mu_j - d^i \hat{\sigma} \sqrt{n_i^{-1} + n_j^{-1}}
\text{ for all } j \neq i
\right\}
=
1 - \alpha.
$$

After standardization, the defining probability is

$$
p_i^{\mathrm{con}}(d)
=
\int_0^\infty
\int_{-\infty}^{+\infty}
\prod_{j \neq i}
\Phi\left(
\frac{\lambda_{ij} z + d s}{\sqrt{1 - \lambda_{ij}^2}}
\right)
\, d\Phi(z)\, \gamma_\nu(s)\, ds,
$$

and $d^i$ is the unique solution to

$$
p_i^{\mathrm{con}}(d^i) = 1 - \alpha.
$$

The package brackets the root with

$$
t_{\alpha,\nu} \leq d^i \leq t_{\alpha/(k-1),\nu},
$$

then solves the scalar equation with the Illinois variant of regula falsi.

## Unconstrained Edwards-Hsu critical value $|d|^i$

For the direct unbalanced Edwards-Hsu method, the critical value is two-sided:

$$
P\left\{
\left|
\hat{\mu}_i - \mu_i - (\hat{\mu}_j - \mu_j)
\right|
<
|d|^i \hat{\sigma} \sqrt{n_i^{-1} + n_j^{-1}}
\text{ for all } j \neq i
\right\}
=
1 - \alpha.
$$

The corresponding integral is

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
\, d\Phi(z)\, \gamma_\nu(s)\, ds,
$$

and $|d|^i$ is the unique solution to

$$
p_i^{\mathrm{EH}}(|d|^i) = 1 - \alpha.
$$

The package brackets this root with

$$
t_{\alpha/2,\nu} \leq |d|^i \leq t_{\alpha/[2(k-1)],\nu},
$$

again using the Illinois method for the final solve.

## Numerical implementation

The numerical integration strategy mirrors Hsu's Appendix D.2:

- For finite $\nu$, the inner normal integral is evaluated with 48-point Gauss-Hermite quadrature.
- The outer integral over $S = \hat{\sigma} / \sigma$ is evaluated with 48-point Gauss-Legendre quadrature.
- The outer domain $(0, \infty)$ is truncated to the central $1 - 10^{-6}$ probability interval of $S$.
- When $\nu > 240$, the code switches to the infinite-degrees-of-freedom approximation and treats $S \equiv 1$.

## Practical interpretation

- `constrained_critical_value` is the one-sided ingredient used by `constrained_mcb`.
- `unconstrained_edwards_hsu_critical_value` is the two-sided ingredient used by `unconstrained_mcb_edwards_hsu`.
- In an unbalanced design there is no single universal table indexed only by $(\alpha, k, \nu)$; the sample-size pattern matters through the $\lambda_{ij}$ values.
