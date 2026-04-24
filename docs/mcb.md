# MCB Guide

Multiple Comparisons with the Best (MCB) is for experiments where the primary question is not every pairwise difference, but which treatments are compatible with being best and how far the others may be from the best.

This package follows the unbalanced one-way normal model in Hsu (1996), Chapter 4. The implementation covers:

- constrained MCB from Section 4.1.4, targeting $\theta_i = \mu_i - \max_{j \neq i} \mu_j$;
- the direct Edwards-Hsu unconstrained method from Section 4.2.4.2, targeting $\eta_i = \mu_i - \max_j \mu_j$;
- critical-value computation based on Appendix D.2, with balanced-case regression tests from Appendix E.

## Reading Path

Read the pages in this order if you want the full statistical story:

1. [Model and notation](model_notation.md): the one-way model, pooled MSE, degrees of freedom, and the two MCB parameterizations.
2. [Constrained MCB](constrained_mcb.md): Hsu's sharper candidate-best procedure.
3. [Edwards-Hsu MCB](edwards_hsu.md): direct intervals for the gap to the global best.
4. [SAT worked example](sat_example.md): a source-grounded reproduction of Hsu's Chapter 4 example.
5. [MCB critical values](mcb_critical_values.md): the numerical integrations and root solves used by the package.
6. [Assumptions and scope](assumptions.md): what the coverage statement relies on and what is not implemented.

## Which Procedure Should I Use?

Use `constrained_mcb` when your main deliverable is a candidate-best set and intervals for each arm against its best competitor:

$$
\theta_i = \mu_i - \max_{j \neq i} \mu_j.
$$

This is the procedure Hsu emphasizes for identifying treatments that are not best while preserving sharper intervals. Its intervals are constrained to include 0.

Use `unconstrained_mcb_edwards_hsu` when you need direct lower and upper bounds for each arm's gap to the global best:

$$
\eta_i = \mu_i - \max_{1 \leq j \leq k} \mu_j \leq 0.
$$

These intervals are not constrained to include 0. If an upper bound is strictly negative, that arm is inferred to be worse than the true best.

Both procedures return `MCBResult`, so downstream code can inspect `labels`, `estimates`, `lower`, `upper`, `critical_values`, and `candidate_best_set` in the same way.
