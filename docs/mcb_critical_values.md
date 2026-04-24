# MCB Critical Values

This page describes the two critical values exposed by the package:

- `constrained_critical_value` for Hsu's constrained MCB critical value $d^i$
- `unconstrained_edwards_hsu_critical_value` for the direct Edwards-Hsu critical value $|d|^i$

Both are computed arm-by-arm in the unbalanced one-way model, so the value for a focal arm depends on that arm's full pattern of pairwise sample-size ratios. The formulas follow Hsu's Appendix D.2 quadrature strategy for the probabilities in Section 4.1.4 and Section 4.2.4.2.

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

Let $\gamma_\nu(s)$ denote the density of $S = \hat{\sigma}/\sigma$ on $(0,\infty)$.

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

The package brackets the root with the same one-sided Dunnett bounds described in Appendix D.2:

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

The package brackets this root with the same two-sided Dunnett bounds described in Appendix D.2:

$$
t_{\alpha/2,\nu} \leq |d|^i \leq t_{\alpha/[2(k-1)],\nu},
$$

again using the Illinois method for the final solve.

## Numerical Implementation

The numerical integration strategy mirrors Hsu's Appendix D.2:

- For finite $\nu < 10^6$, the inner normal integral is evaluated with 48-point Gauss-Hermite quadrature.
- The outer integral over $S = \hat{\sigma} / \sigma$ is evaluated with 48-point Gauss-Legendre quadrature.
- The outer domain $(0, \infty)$ is truncated to the central $1 - 10^{-6}$ probability interval of $S$.
- When $\nu = \infty$ is supplied explicitly, or when finite $\nu \geq 10^6$, the code
  uses the infinite-degrees-of-freedom approximation and treats $S \equiv 1$ because
  the finite-$\nu$ correction is negligible relative to the solver tolerance; finite
  $\nu$ below this high-df stability switch are not collapsed to Hsu's historical
  `dinfnu = 240` code shortcut.
- The root solver is the Illinois variant of regula falsi.
- A candidate root is accepted when the probability is within $2\times 10^{-5}$ of $1-\alpha$, or when the bracketing half-width is at most $5\times 10^{-4}$.

## Balanced Tables and Unbalanced Designs

Appendix E tabulates balanced-design values for $d$ and $|d|$. In a balanced design,

$$
\lambda_{ij} = \frac{1}{\sqrt{2}}
$$

for every comparison, so the critical value depends only on $\alpha$, $k$, and $\nu$.

In an unbalanced design, there is no single table indexed only by $(\alpha,k,\nu)$. The vector $(\lambda_{ij}: j \neq i)$ changes by focal arm, and so do the critical values.

The package test suite checks the balanced case against Hsu Appendix E Tables E.2, E.3, E.5, and E.6.

## Practical Interpretation

- `constrained_critical_value` is the one-sided ingredient used by `constrained_mcb`.
- `unconstrained_edwards_hsu_critical_value` is the two-sided ingredient used by `unconstrained_mcb_edwards_hsu`.
- For $k=2$, the constrained critical value reduces to the upper $\alpha$ quantile of the $t$ distribution, and the Edwards-Hsu critical value reduces to the upper $\alpha/2$ quantile.
