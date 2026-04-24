# Assumptions and Scope

The simultaneous coverage statements in `mcb` are finite-sample results for Hsu's one-way normal model. The statistical assumptions matter.

## Required Model Conditions

The exact coverage claims assume:

- independent observations within and across arms;
- a common variance $\sigma^2$ across arms;
- normally distributed errors;
- a correctly computed pooled MSE with $\nu = \sum_i(n_i-1)$ error degrees of freedom;
- the larger-is-better direction of comparison.

If your metric is smaller-is-better, transform it before analysis, for example by analyzing `-x`.

## Practical A/B/n Use

The package can be useful for A/B/n workflows with large samples, but exact finite-sample coverage is not guaranteed for heavily skewed, zero-inflated, heteroscedastic, or dependent observations. In those cases, the normal one-way model is an approximation.

For revenue/session and similar metrics, treat the output as most defensible when arm means are approximately normal and the common-variance approximation is credible.

## Ties and the Meaning of "Best"

Hsu's derivations allow ties to be broken arbitrarily for proof purposes. Interpretation depends on your scientific stance:

- If ties are treated as practically impossible, an arm with constrained lower bound 0 can be interpreted as selected as best.
- If ties are possible, a bound at 0 means the arm remains compatible with being one of the best arms.

Either way, a candidate-best set is a confidence-set object, not a probability ranking.

## What Is Not Implemented

`mcb` currently does not implement:

- smaller-is-better formulas directly;
- sample-size planning from Appendix C;
- Tukey or Tukey-Kramer all-pairwise procedures;
- comparison-with-control procedures as standalone APIs;
- nonparametric MCB methods;
- general linear model MCB extensions.

The implemented API is intentionally focused on Hsu's unbalanced one-way constrained MCB and the Edwards-Hsu direct unconstrained MCB method.
