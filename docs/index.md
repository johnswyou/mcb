# `mcb`

`mcb` implements Hsu's Multiple Comparisons with the Best (MCB) for the one-way normal model, with equal or unequal sample sizes per arm. From the per-arm data of an experiment it reports which arms could be the best and how far each arm may be from the best, as simultaneous confidence intervals.

The package does three things:

- `constrained_mcb` computes simultaneous intervals for $\mu_i - \max_{j \neq i} \mu_j$, arm $i$'s mean minus the best of the other arms' means, and the set of arms that could be the best.
- `unconstrained_mcb_edwards_hsu` computes simultaneous intervals for $\mu_i - \max_j \mu_j$, arm $i$'s mean minus the best mean.
- `constrained_mcb_sample_size` and `unconstrained_mcb_edwards_hsu_sample_size` choose the per-arm sample size of a balanced experiment before it runs.

Source repository: <https://github.com/johnswyou/mcb>

```{toctree}
:maxdepth: 2
:caption: Guide

quickstart
mcb
model_notation
constrained_mcb
edwards_hsu
sat_example
mcb_critical_values
sample_size
assumptions
api
```

## References

- Hsu, J. C. (1996). *Multiple Comparisons: Theory and Methods*. Chapman & Hall. The procedures are in Chapter 4, the sample-size computation in Appendix C, the critical-value computation in Appendix D.2, and the tables used as regression tests in Appendix E.
- Hsu, J. C. (1988). Sample size computation for designing multiple comparison experiments. *Computational Statistics and Data Analysis*, 7:79-91.
- Edwards, D. G. and Hsu, J. C. (1983). Multiple comparisons with the best treatment. *Journal of the American Statistical Association*, 78:965-971.
